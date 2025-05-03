#' Bayesian Lalonde Decision Framework
#'
#' This function implements a Bayesian decision-making framework based on the Lalonde criteria.
#' It calculates posterior distributions for each arm and evaluates posterior inference and decision criteria based on the given endpoint type (continuous or binary).
#'
#' @param endpoint A string. Type of endpoint ("continuous" or "binary").
#' @param data_summary A data frame. Summary of the individual-level data. For continuous endpoints, each row represents the number of patients, the mean of the normal distribution, and the standard deviation of the normal distribution for one arm in one simulation repetition. For binary endpoints, each row represents the number of patients and the number of responses for one arm in one repetition of the simulation.
#' @param prior_params A list. Prior distribution parameters for both the treatment and control arms. It should include delta, w, a, and b values for both arms.
#' @param lrv A scalar. The lower reference value for decision-making.
#' @param tv A scalar. The target value for decision-making.
#' @param fgr A scalar. False go rate for decision-making.
#' @param fsr A scalar. False stop rate for decision-making.
#' @param arm_names A list. Names of the arms in the simulation (e.g., treatment, control, etc.).
#' @param posterior_infer A logical. If TRUE, the function computes posterior inference for credible intervals and posterior probabilities.
#' @param Lalonde_decision A logical. If TRUE, the function makes decisions based on the Lalonde framework.
#' @param EXP_TRANSFORM A logical. If TRUE, the function exponentiates the results, which is useful for log-transformed data.
#'
#' @return A list. The output includes:
#'   - `post_params`: A data frame containing the posterior distribution parameters for each arm.
#'   - `post_est_ci`: A data frame containing evaluating metrics for the posterior distribution.
#'   - `post_inference`: A data frame containing the posterior probabilities and credible interval, and the decisions based on the Lalonde framework.
#' @export
#' @import dplyr
#' @import tidyr
#' @import RBesT
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import data.table
bayesian_lalonde_decision <- function(endpoint, data_summary,
                                      prior_params, lrv = 0, tv = 0, fgr = 0.2, fsr = 0.1, arm_names,
                                      posterior_infer = TRUE, Lalonde_decision = TRUE,
                                      EXP_TRANSFORM = FALSE) {

  if (!"nsim" %in% names(data_summary)) {
    stop("data_summary must contain an 'nsim' column.")
  }
  all_nsims <- sort(unique(data_summary$nsim))
  N_sims <- length(all_nsims)

  if (N_sims == 0) {
    stop("No simulations found in data_summary.")
  }
  cat(sprintf("Starting optimization for %d simulations.\n", N_sims))

  # --- Cluster Setup ---
  is_ci <- Sys.getenv("CI") == "true" ||
    Sys.getenv("TRAVIS") == "true" ||
    Sys.getenv("GITHUB_ACTIONS") == "true" ||
    Sys.getenv("_R_CHECK_LIMIT_CORES_") == "true"

  ncore <- if (is_ci) {
    2
  } else {
    max(1, detectCores() - 1)
  }

  cat(sprintf("Setting up parallel backend with %d cores.\n", ncore))
  cl <- makeCluster(ncore)

  # --- Export Necessary Data and Functions ---
  required_funcs <- c("posterior_distribution", "convert_RBesT_mix",
                      "posterior_inference", "posterior_prob")

  func_envir <- environment()
  global_env <- globalenv()

  funcs_to_export <- Filter(function(f) exists(f, envir = func_envir, inherits = TRUE) || exists(f, envir = global_env, inherits = TRUE), required_funcs)
  if(length(funcs_to_export) != length(required_funcs)){
    warning("Could not find all required functions to export to cluster.")
  }

  required_vars <- c("endpoint", "data_summary", "lrv", "tv", "fgr", "fsr", "arm_names",
                     "prior_params", "EXP_TRANSFORM", "posterior_infer")

  clusterExport(cl, c(required_vars, funcs_to_export), envir = environment())

  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(RBesT)
    library(data.table) # Ensure data.table is available for rbindlist if used
  })

  registerDoParallel(cl)
  cat("Cluster setup complete. Starting parallel computation...\n")

  # --- Combined Parallel Loop ---
  results_list <- foreach(
    nrep = all_nsims,
    .combine = function(a, b) data.table::rbindlist(list(a, b)),
    .errorhandling = 'stop',
    .packages = c("RBesT", "dplyr", "tidyr", "data.table")
  ) %dopar% {

    # 1. Get Data for current simulation
    current_sim_data <- data_summary[data_summary$nsim == nrep, ]

    historical_ctrl_list <- if (any(names(arm_names) == "control_h")) current_sim_data[startsWith(names(current_sim_data), "control_h.")] else NULL
    historical_trt_list <- if (any(names(arm_names) == "treatment_h")) current_sim_data[startsWith(names(current_sim_data), "treatment_h.")] else NULL
    current_ctrl_list <- current_sim_data[startsWith(names(current_sim_data), "control.")]
    current_trt_list  <- current_sim_data[startsWith(names(current_sim_data), "treatment.")]

    if (nrow(current_ctrl_list) == 0 || nrow(current_trt_list) == 0) {
      warning(paste("Missing current control or treatment data for nsim =", nrep))
      return(NULL)
    }

    current_ctrl <- setNames(as.list(current_ctrl_list[1, ]), sub("^control\\.", "", names(current_ctrl_list)))
    current_trt  <- setNames(as.list(current_trt_list[1, ]), sub("^treatment\\.", "", names(current_trt_list)))
    historical_ctrl <- if (!is.null(historical_ctrl_list) && ncol(historical_ctrl_list) > 0) setNames(as.list(historical_ctrl_list[1, ]), sub("^control_h\\.", "", names(historical_ctrl_list))) else NULL
    historical_trt  <- if (!is.null(historical_trt_list) && ncol(historical_trt_list) > 0) setNames(as.list(historical_trt_list[1, ]), sub("^treatment_h\\.", "", names(historical_trt_list))) else NULL

    # 2. Calculate Posterior Parameters
    post_c <- posterior_distribution(endpoint, current = current_ctrl, historical = historical_ctrl,
                                     delta = prior_params[["control.delta"]], w = prior_params[["control.w"]],
                                     a = prior_params[["control.a"]], b = prior_params[["control.b"]],
                                     a0 = prior_params[["control.a0"]], b0 = prior_params[["control.b0"]],
                                     theta0 = prior_params[["control.theta0"]], s0 = prior_params[["control.s0"]],
                                     ess_h = prior_params[["control.ess_h"]])

    post_t <- posterior_distribution(endpoint, current = current_trt, historical = historical_trt,
                                     delta = prior_params[["treatment.delta"]], w = prior_params[["treatment.w"]],
                                     a = prior_params[["treatment.a"]], b = prior_params[["treatment.b"]],
                                     a0 = prior_params[["treatment.a0"]], b0 = prior_params[["treatment.b0"]],
                                     theta0 = prior_params[["treatment.theta0"]], s0 = prior_params[["treatment.s0"]],
                                     ess_h = prior_params[["treatment.ess_h"]])

    post_params_i <- c(unlist(post_t$post), unlist(post_c$post))
    names(post_params_i) <- c(paste0("treatment.", names(post_t$post), "_post"), paste0("control.", names(post_c$post), "_post"))

    prior_ws_i <- c(post_t$w_prior, post_c$w_prior)
    names(prior_ws_i) <- c("treatment.w_prior", "control.w_prior")

    params_results <- c(nsim = nrep, prior_ws_i, post_params_i)

    # 3. Calculate Posterior Inference
    # Convert posterior parameters to RBesT mix objects
    post_t_mix <- convert_RBesT_mix(post = post_t$post, endpoint = endpoint)
    post_c_mix <- convert_RBesT_mix(post = post_c$post, endpoint = endpoint)

    post_est_ci_i <- posterior_inference(post1 = post_t_mix, post2 = post_c_mix,
                                         quantiles = c(0.025, 0.975),
                                         EXP_TRANSFORM = EXP_TRANSFORM)
    names(post_est_ci_i) <- paste0(names(post_est_ci_i), "_95ci")

    # 4. Optional inference for Lalonde
    lalonde_criteria_i <- NULL
    post_prob_i <- NULL
    if (posterior_infer) {
      quantiles <- if (lrv < tv) c(low = fgr, upper = 1 - fsr) else c(low = fsr, upper = 1 - fgr)
      lalonde_criteria_i <- posterior_inference(post1 = post_t_mix, post2 = post_c_mix,
                                                quantiles = quantiles,
                                                EXP_TRANSFORM = EXP_TRANSFORM)
      names(lalonde_criteria_i) <- paste0(names(lalonde_criteria_i), "_lalonde")

      lrv_adj <- if (EXP_TRANSFORM) log(lrv) else lrv
      tv_adj <- if (EXP_TRANSFORM) log(tv) else tv

      pr_m <- posterior_prob(post1 = post_t_mix, value = lrv_adj, post2 = post_c_mix,
                             range_type = ifelse(lrv < tv, "less", "greater"))
      pr_l <- posterior_prob(post1 = post_t_mix, value = ifelse(lrv < tv, lrv_adj, tv_adj),
                             post2 = post_c_mix, range_type = "between",
                             value2 = ifelse(lrv < tv, tv_adj, lrv_adj))
      pr_t <- posterior_prob(post1 = post_t_mix, value = tv_adj, post2 = post_c_mix,
                             range_type = ifelse(lrv < tv, "greater", "less"))

      post_prob_i <- data.frame(pr_m = pr_m, pr_l = pr_l, pr_t = pr_t)
    } else {
      lalonde_placeholder_names <- c("compare_ci_l_lalonde", "compare_ci_u_lalonde", "post1_est_lalonde", "post2_est_lalonde", "diff_est_lalonde", "compare_est_lalonde")
      lalonde_criteria_i <- setNames(data.frame(matrix(NA, ncol = length(lalonde_placeholder_names), nrow = 1)), lalonde_placeholder_names)

      prob_placeholder_names <- c("pr_m", "pr_l", "pr_t")
      post_prob_i <- setNames(data.frame(matrix(NA, ncol = length(prob_placeholder_names), nrow = 1)), prob_placeholder_names)

    }

    result_df_i <- data.table::as.data.table(cbind(
      as.data.frame(t(params_results)),
      as.data.frame(post_est_ci_i),
      as.data.frame(lalonde_criteria_i),
      as.data.frame(post_prob_i)
    ))
    result_df_i
  }

  cat("Parallel computation finished. Stopping cluster.\n")
  stopCluster(cl)
  registerDoSEQ()

  # --- Post-processing ---
  if (!is.data.frame(results_list)) {
    stop("Parallel computation did not return a data frame.")
  }

  numeric_cols <- setdiff(names(results_list), "nsim") # Assuming all others should be numeric
  results_list <- suppressWarnings(results_list %>% mutate(across(all_of(numeric_cols), as.numeric)))

  post_params_df <- results_list %>% select(nsim, ends_with("_post"), ends_with("w_prior"))
  post_est_ci_df <- results_list %>% select(nsim, ends_with("_95ci"))
  post_inference_df <- results_list %>% select(nsim, ends_with("_lalonde"), starts_with("pr_"))

  # 5. Add decision columns if requested
  if (Lalonde_decision && posterior_infer) {
    cat("Start making decisions based on Lalonde framework.\n")

    compare_ci_l_col <- "compare_ci_l_lalonde"
    compare_ci_u_col <- "compare_ci_u_lalonde"

    if(!all(c(compare_ci_l_col, compare_ci_u_col, "pr_t", "pr_l") %in% names(post_inference_df))) {
      warning("Required columns for Lalonde decision not found in results. Skipping decision step.")
      post_inference_df$decision_pr <- NA
      post_inference_df$decision_ci <- NA
    } else {
      post_inference_df <- post_inference_df %>%
        mutate(
          decision_pr = case_when(
            (pr_t > fsr) & ((pr_t + pr_l) > 1 - fgr) ~ "go",
            (pr_t > fsr) & ((pr_t + pr_l) <= 1 - fgr) ~ "consider",
            pr_t <= fsr ~ "no-go",
            TRUE ~ NA_character_
          ),
          decision_ci = if (lrv < tv) {
            case_when(
              !!sym(compare_ci_u_col) > tv & !!sym(compare_ci_l_col) > lrv ~ "go",
              !!sym(compare_ci_u_col) > tv & !!sym(compare_ci_l_col) <= lrv ~ "consider",
              !!sym(compare_ci_u_col) <= tv ~ "no-go",
              TRUE ~ NA_character_
            )
          } else {
            case_when(
              !!sym(compare_ci_l_col) < tv & !!sym(compare_ci_u_col) < lrv ~ "go",
              !!sym(compare_ci_l_col) < tv & !!sym(compare_ci_u_col) >= lrv ~ "consider",
              !!sym(compare_ci_l_col) >= tv ~ "no-go",
              TRUE ~ NA_character_
            )
          }
        )
    }
  } else if (Lalonde_decision && !posterior_infer) {
    post_inference_df$decision_pr <- NA_character_
    post_inference_df$decision_ci <- NA_character_
  }

  cat("Function execution complete.\n")

  return(list(post_params = post_params_df,
              post_est_ci = post_est_ci_df,
              post_inference = post_inference_df))
}


