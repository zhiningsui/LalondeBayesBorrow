#' Bayesian Decision Framework based on Lalonde Criteria
#'
#' @description
#' This function simulates a Bayesian trial across multiple scenarios (simulations),
#' calculating posterior distributions for treatment and control arms (with optional
#' historical borrowing), performing Bayesian inference (estimates, credible
#' intervals, posterior probabilities), and applying decision rules based on
#' the Lalonde framework.
#'
#' @details
#' It leverages `posterior_distribution` to calculate individual arm posteriors,
#' `posterior_inference` to summarize these posteriors and their difference/ratio,
#' and `posterior_prob` to calculate probabilities against thresholds.
#'
#' The function is designed to work with data summarized per arm per simulation
#' repetition, typically the output of a data generation step. Parallel processing
#' is used to speed up calculations across simulations.
#'
#' @param endpoint A character string. The type of endpoint. Must be `'continuous'` or `'binary'`.
#'   This determines the likelihood and prior structure used by `posterior_distribution`.
#' @param data_summary A data frame. For a single analysis, this can be a single-row
#'   data frame without an `nsim` column. For simulations, it must contain an `nsim`
#'   column. Data for each arm and its parameters should be in columns prefixed
#'   by the arm name and a dot (e.g., `control.n`, `treatment.mu_hat`, `control_h.count`).
#'   * For `'continuous'` endpoint, required columns per arm (e.g., `control.`)
#'       are `n` (sample size), `mu_hat` (sample mean), `s` (standard deviation of sample mean).
#'   * For `'binary'` endpoint, required columns per arm (e.g., `control.`)
#'       are `n` (sample size), `count` (number of responses/events).
#'   Historical arm data (e.g., `control_h.n`, `treatment_h.count`) are optional.
#' @param prior_params A named list. Contains prior distribution and borrowing
#'   parameters for each arm (treatment and control, potentially historical).
#'   Names must follow the convention `[arm_name].[parameter]`, matching the
#'   parameters expected by `posterior_distribution`.
#'   * For `'continuous'` arms, includes `[arm_name].theta0`, `[arm_name].s0`.
#'   * For `'binary'` arms, includes `[arm_name].a`, `[arm_name].b`, `[arm_name].a0`, `[arm_name].b0`.
#'   * For arms with potential historical borrowing, includes `[arm_name].delta`, `[arm_name].w`, `[arm_name].ess_h`.
#'   Example names: `treatment.delta`, `control.w`, `treatment.theta0`, `control.s0`,
#'   `treatment.a`, `control.b`, `treatment.a0`, `control.b0`, `control.ess_h`.
#' @param lrv A scalar numeric value. The Lower Reference Value for the comparison
#'   metric (difference or ratio). Used in the Lalonde decision criteria.
#' @param tv A scalar numeric value. The Target Value for the comparison metric.
#'   Used in the Lalonde decision criteria.
#' @param fgr A scalar numeric value between 0 and 1. The False Go Rate. Used in the
#'   Lalonde probability and quantile-based decision criteria.
#' @param fsr A scalar numeric value between 0 and 1. The False Stop Rate. Used in the
#'   Lalonde probability and quantile-based decision criteria.
#' @param arm_names A character vector. A list of arm names present in the simulation.
#'   This parameter is used internally to identify which arm data columns to look for
#'   in `data_summary` and which prior parameters to expect in `prior_params`.
#'   Expected names are typically `"treatment"`, `"control"`, `"treatment_h"`, `"control_h"`.
#' @param posterior_infer A logical value. If `TRUE`, compute posterior inference statistics
#'   (95% credible intervals and posterior probabilities `pr_m`, `pr_l`, `pr_t`). Defaults to `TRUE`.
#' @param Lalonde_decision A logical value. If `TRUE`, apply the Lalonde decision rules
#'   based on the calculated posterior probabilities and quantiles. Defaults to `TRUE`.
#'   Note: Requires `posterior_infer = TRUE` to make decisions.
#' @param EXP_TRANSFORM A logical value. If `TRUE`, results for individual arm means,
#'   credible intervals, and the comparison estimate (`est_compare`) are exponentiated.
#'   This implies that the underlying parameters in `post_t_mix` and `post_c_mix`
#'   (and hence `mu_hat`, `mu` in `data_summary` and `prior_params`) are on a log scale,
#'   and the comparison is interpreted as a ratio. Defaults to `FALSE`.
#'
#' @return A list containing three data frames, each row corresponding to a simulation
#'   repetition (`nsim`):
#'   * `post_params`: Contains the raw posterior parameters (e.g., `w_post`, `mu1_post`, `sigma1_post`, etc.)
#'       for the treatment and control arms, prefixed with `treatment.` or `control.`,
#'       and the prior weights (`treatment.w_prior`, `control.w_prior`).
#'   * `post_est_ci`: Contains standard posterior inference (mean, SD, 95% CI)
#'       for each arm and their comparison, postfixed with `_95ci`. Values are
#'       exponentiated if `EXP_TRANSFORM = TRUE` (except SDs).
#'   * `post_inference`: Contains inference results relevant to the Lalonde
#'       decision framework. Includes posterior probabilities (`pr_m`, `pr_l`, `pr_t`)
#'       and comparison quantiles derived from `fgr`/`fsr` (`_lalonde`). If
#'       `Lalonde_decision = TRUE` and `posterior_infer = TRUE`, includes decision
#'       columns (`decision_pr`, `decision_ci`).
#'
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom data.table as.data.table rbindlist
#' @importFrom dplyr mutate across case_when select everything
#' @importFrom tidyr pivot_longer
#' @importFrom stats setNames
#' @seealso `posterior_distribution`, `posterior_inference`,
#'   `posterior_prob`, `convert_RBesT_mix`,
#'   `RBesT::mix`, `RBesT::mixnorm`, `RBesT::mixbeta`,
#'   `RBesT::pmixdiff`, `RBesT::qmixdiff`
#'
#' @export
bayesian_lalonde_decision <- function(endpoint, data_summary,
                                      prior_params, lrv = 0, tv = 0, fgr = 0.2, fsr = 0.1, arm_names,
                                      posterior_infer = TRUE, Lalonde_decision = TRUE,
                                      EXP_TRANSFORM = FALSE) {

  # --- Input Validation ---
  valid_endpoints <- c("continuous", "binary")
  if (!is.character(endpoint) || length(endpoint) != 1 || !(endpoint %in% valid_endpoints)) {
    stop(paste("Input 'endpoint' must be a single string:", paste(valid_endpoints, collapse = ", ")))
  }
  if (!is.data.frame(data_summary)) {
    stop("Input 'data_summary' must be a data frame.")
  }

  is_simulation <- "nsim" %in% names(data_summary) && length(unique(data_summary$nsim)) > 1

  # If not a simulation, ensure only one row of data is present
  if (!is_simulation && nrow(data_summary) > 1) {
    warning("`data_summary` has more than one row but is not a simulation (no `nsim` column or only one `nsim` value). Only the first row will be used for analysis.")
    data_summary <- data_summary[1, , drop = FALSE]
  }
  # Add nsim column if it's a single analysis
  if (!"nsim" %in% names(data_summary)) {
    data_summary$nsim <- 1
  }
  if (!is.list(prior_params) || is.null(names(prior_params))) {
    stop("Input 'prior_params' must be a named list.")
  }
  if (!is.numeric(lrv) || length(lrv) != 1 || !is.finite(lrv)) {
    stop("Input 'lrv' must be a single finite numeric value.")
  }
  if (!is.numeric(tv) || length(tv) != 1 || !is.finite(tv)) {
    stop("Input 'tv' must be a single finite numeric value.")
  }
  if (!is.numeric(fgr) || length(fgr) != 1 || !is.finite(fgr) || fgr < 0 || fgr > 1) {
    stop("Input 'fgr' must be a single finite numeric value between 0 and 1.")
  }
  if (!is.numeric(fsr) || length(fsr) != 1 || !is.finite(fsr) || fsr < 0 || fsr > 1) {
    stop("Input 'fsr' must be a single finite numeric value between 0 and 1.")
  }
  if (!is.character(arm_names) || length(arm_names) == 0 || !all(nzchar(arm_names))) {
    stop("Input 'arm_names' must be a non-empty character vector with non-empty strings.")
  }
  if (!is.logical(posterior_infer) || length(posterior_infer) != 1) {
    stop("Input 'posterior_infer' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(Lalonde_decision) || length(Lalonde_decision) != 1) {
    stop("Input 'Lalonde_decision' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(EXP_TRANSFORM) || length(EXP_TRANSFORM) != 1) {
    stop("Input 'EXP_TRANSFORM' must be a single logical value (TRUE or FALSE).")
  }

  all_nsims <- sort(unique(data_summary$nsim))
  N_sims <- length(all_nsims)
  cat(sprintf("Processing %d analysis/simulation run(s).\n", N_sims))

  # --- Cluster Setup for Simulations ---
  if (is_simulation) {
    is_restricted <- Sys.getenv("CI") == "true" ||
      Sys.getenv("TRAVIS") == "true" ||
      Sys.getenv("GITHUB_ACTIONS") == "true" ||
      Sys.getenv("_R_CHECK_LIMIT_CORES_") == "true"

    ncore <- if (is_restricted) 2 else max(1, parallel::detectCores() - 1)
    cl <- tryCatch({
      parallel::makeCluster(ncore)
      cat(sprintf("Setting up parallel backend with %d cores.\n", ncore))
    }, error = function(e) {
      warning(paste("Could not create parallel cluster with", ncore, "cores. Falling back to sequential processing:", e$message))
      return(NULL)
    })

    if (!is.null(cl)) {
      required_funcs <- c("posterior_distribution", "convert_RBesT_mix",
                          "posterior_inference", "posterior_prob")
      available_funcs <- sapply(required_funcs, function(f) exists(f, where = environment(NULL), inherits = TRUE))
      if (!all(available_funcs)) {
        missing_funcs <- required_funcs[!available_funcs]
        stop(paste("Required functions for parallel export are missing:", paste(missing_funcs, collapse = ", "),
                   ". Please ensure they are defined in the current session.", sep=""))
      }
      required_vars_for_export <- c("endpoint", "lrv", "tv", "fgr", "fsr", "arm_names",
                                    "prior_params", "EXP_TRANSFORM", "posterior_infer")
      parallel::clusterExport(cl, c(required_vars_for_export, required_funcs), envir = environment(NULL))

      parallel::clusterEvalQ(cl, {
        library(dplyr)
        library(tidyr)
        library(RBesT)
        library(data.table)
      })

      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      cat("Cluster setup complete. Starting parallel computation...\n")
      parallel_enabled <- TRUE
    } else {
      foreach::registerDoSEQ()
      cat("Proceeding with sequential computation.\n")
      parallel_enabled <- FALSE
    }
  } else {
    foreach::registerDoSEQ()
    parallel_enabled <- FALSE
  }

  # --- Combined Parallel Loop ---
  data_summary_dt <- data.table::as.data.table(data_summary)

  results_list <- foreach::foreach(
    nrep_val = all_nsims,
    .combine = function(a, b) data.table::rbindlist(list(a, b)),
    .errorhandling = 'stop',
    .packages = c("RBesT", "dplyr", "tidyr", "data.table")
  ) %dopar% {

    # 1. Get Data for current simulation
    current_sim_data_dt <- data_summary_dt[nsim == nrep_val, ]

    # Extract arm data using prefixes defined by arm_names
    arm_data_list <- list()
    all_params_found <- TRUE # Flag to check if all expected params for current arms are found
    for (arm_nm in c("treatment", "control")) { # Assuming 'treatment' and 'control' current arms always exist
      prefix <- paste0(arm_nm, ".")
      arm_cols <- names(current_sim_data_dt)[startsWith(names(current_sim_data_dt), prefix)]
      if (length(arm_cols) > 0) {
        # Remove prefix for posterior_distribution
        arm_data <- setNames(as.list(current_sim_data_dt[, .SD, .SDcols = arm_cols][1,]), sub(prefix, "", arm_cols))
        arm_data_list[[arm_nm]] <- arm_data

        # Check if required parameters for the endpoint are present in the extracted data
        required_data_params <- if (endpoint == "continuous") c("n", "mu_hat", "s") else c("n", "count")
        if (!all(required_data_params %in% names(arm_data))) {
          warning(paste("Missing required data parameters for arm '", arm_nm, "' for nsim =", nrep_val,
                        ". Expected:", paste(required_data_params, collapse = ", "),
                        ". Found:", paste(names(arm_data), collapse = ", "), sep=""))
          all_params_found <- FALSE # Mark as missing
          break
        }
      } else {
        warning(paste("Missing current data columns for arm '", arm_nm, "' for nsim =", nrep_val))
        all_params_found <- FALSE # Mark as missing
        break
      }
    }

    # Check for historical arms if listed in arm_names
    historical_data_list <- list()
    for (arm_nm in c("treatment_h", "control_h")) { # Assuming these are the only possible historical prefixes
      if (arm_nm %in% arm_names) { # Only look if the historical arm is expected
        prefix <- paste0(arm_nm, ".")
        arm_cols <- names(current_sim_data_dt)[startsWith(names(current_sim_data_dt), prefix)]
        if (length(arm_cols) > 0) {
          historical_data_list[[arm_nm]] <- setNames(as.list(current_sim_data_dt[, .SD, .SDcols = arm_cols][1,]), sub(prefix, "", arm_cols))
        } else {
          cat("check1")
          historical_data_list[[arm_nm]] <- NULL
        }
      } else {
        cat("check2")
        historical_data_list[[arm_nm]] <- NULL
      }
    }

    if (!all_params_found) {
      return(NULL)
    }

    # 2. Calculate Posterior Parameters for Treatment and Control
    posterior_results <- tryCatch({
      post_c_res <- posterior_distribution(endpoint = endpoint,
                                           current = arm_data_list[["control"]],
                                           historical = historical_data_list[["control_h"]],
                                           delta_gate = prior_params[["control.delta_gate"]],
                                           delta_SAM = prior_params[["control.delta_SAM"]],
                                           w = prior_params[["control.w"]],
                                           a = prior_params[["control.a"]], b = prior_params[["control.b"]],
                                           a0 = prior_params[["control.a0"]], b0 = prior_params[["control.b0"]],
                                           theta0 = prior_params[["control.theta0"]], s0 = prior_params[["control.s0"]],
                                           ess_h = prior_params[["control.ess_h"]])

      post_t_res <- posterior_distribution(endpoint = endpoint,
                                           current = arm_data_list[["treatment"]],
                                           historical = historical_data_list[["treatment_h"]],
                                           delta_gate = prior_params[["treatment.delta_gate"]],
                                           delta_SAM = prior_params[["treatment.delta_SAM"]],
                                           w = prior_params[["treatment.w"]],
                                           a = prior_params[["treatment.a"]], b = prior_params[["treatment.b"]],
                                           a0 = prior_params[["treatment.a0"]], b0 = prior_params[["treatment.b0"]],
                                           theta0 = prior_params[["treatment.theta0"]], s0 = prior_params[["treatment.s0"]],
                                           ess_h = prior_params[["treatment.ess_h"]])

      list(post_c = post_c_res, post_t = post_t_res)

    }, error = function(e) {
      warning(paste("Error calculating posterior distribution for nsim =", nrep_val, ":", e$message))
      return(NULL)
    })

    if (is.null(posterior_results)) {
      return(NULL)
    }

    post_c <- posterior_results$post_c
    post_t <- posterior_results$post_t

    # Prepare raw posterior parameters for output df
    post_params_i <- c(unlist(post_t$post), unlist(post_c$post))
    names(post_params_i) <- c(paste0("treatment.", names(post_t$post), "_post"), paste0("control.", names(post_c$post), "_post"))

    prior_ws_i <- c(treatment.w_prior = post_t$w_prior, control.w_prior = post_c$w_prior)
    params_results <- c(nsim = nrep_val, prior_ws_i, post_params_i)

    # 3. Calculate Posterior Inference (95% CI and standard stats)
    post_inference_results <- tryCatch({
      # Convert posterior parameters to RBesT mix objects
      post_t_mix <- convert_RBesT_mix(post = post_t$post, endpoint = endpoint)
      post_c_mix <- convert_RBesT_mix(post = post_c$post, endpoint = endpoint)

      post_est_ci_i <- posterior_inference(post1 = post_t_mix, post2 = post_c_mix,
                                           quantiles = c(0.025, 0.975),
                                           EXP_TRANSFORM = EXP_TRANSFORM)
      names(post_est_ci_i) <- paste0(names(post_est_ci_i), "_95ci")

      list(post_est_ci = post_est_ci_i, post_t_mix = post_t_mix, post_c_mix = post_c_mix)

    }, error = function(e) {
      warning(paste("Error during standard posterior inference (95% CI) for nsim =", nrep_val, ":", e$message))
      placeholder_95ci_names <- c("est1_95ci", "sd1_95ci", "ci1_l_95ci", "ci1_u_95ci",
                                  "est2_95ci", "sd2_95ci", "ci2_l_95ci", "ci2_u_95ci",
                                  "est_compare_95ci", "sd_compare_95ci", "compare_ci_l_95ci", "compare_ci_u_95ci")
      post_est_ci_i <- setNames(data.frame(matrix(NA_real_, ncol = length(placeholder_95ci_names), nrow = 1)), placeholder_95ci_names)
      list(post_est_ci = post_est_ci_i, post_t_mix = NULL, post_c_mix = NULL)

    })

    post_est_ci_i <- post_inference_results$post_est_ci
    post_t_mix_valid <- post_inference_results$post_t_mix
    post_c_mix_valid <- post_inference_results$post_c_mix

    # 4. Optional Lalonde Inference (Quantiles and Probabilities)
    lalonde_criteria_i <- NULL
    post_prob_i <- NULL
    if (posterior_infer && !is.null(post_t_mix_valid) && !is.null(post_c_mix_valid)) {
      lalonde_quantiles <- if (lrv < tv) c(low = fgr, upper = 1 - fsr) else c(low = fsr, upper = 1 - fgr)
      lalonde_quantiles <- sort(lalonde_quantiles)

      lalonde_criteria_i <- tryCatch({
        posterior_inference(post1 = post_t_mix_valid, post2 = post_c_mix_valid,
                            quantiles = lalonde_quantiles,
                            EXP_TRANSFORM = EXP_TRANSFORM)
      }, error = function(e) {
        warning(paste("Error during Lalonde inference (quantiles) for nsim =", nrep_val, ":", e$message))

        placeholder_lalonde_names <- c("est1_lalonde", "sd1_lalonde", "ci1_l_lalonde", "ci1_u_lalonde",
                                       "est2_lalonde", "sd2_lalonde", "ci2_l_lalonde", "ci2_u_lalonde",
                                       "est_compare_lalonde", "sd_compare_lalonde", "compare_ci_l_lalonde", "compare_ci_u_lalonde")
        setNames(data.frame(matrix(NA_real_, ncol = length(placeholder_lalonde_names), nrow = 1)), placeholder_lalonde_names)
      })
      names(lalonde_criteria_i) <- paste0(names(lalonde_criteria_i), "_lalonde")

      # Calculate probabilities for the Lalonde decision rules (P(<lrv), P(lrv to tv), P(>tv))
      post_prob_i <- tryCatch({

        lrv_adj <- if (EXP_TRANSFORM) log(lrv) else lrv
        tv_adj <- if (EXP_TRANSFORM) log(tv) else tv

        if (!is.finite(lrv_adj) || !is.finite(tv_adj)) {
          stop("Adjusted lrv or tv is not finite after EXP_TRANSFORM.")
        }

        pr_m <- posterior_prob(post1 = post_t_mix_valid, value = lrv_adj, post2 = post_c_mix_valid,
                               range_type = ifelse(lrv < tv, "less", "greater"))

        pr_l <- posterior_prob(post1 = post_t_mix_valid, value = ifelse(lrv < tv, lrv_adj, tv_adj),
                               post2 = post_c_mix_valid, range_type = "between",
                               value2 = ifelse(lrv < tv, tv_adj, lrv_adj))

        pr_t <- posterior_prob(post1 = post_t_mix_valid, value = tv_adj, post2 = post_c_mix_valid,
                               range_type = ifelse(lrv < tv, "greater", "less"))

        data.frame(pr_m = pr_m, pr_l = pr_l, pr_t = pr_t)

      }, error = function(e) {
        warning(paste("Error during Lalonde inference (probabilities) for nsim =", nrep_val, ":", e$message))

        prob_placeholder_names <- c("pr_m", "pr_l", "pr_t")
        setNames(data.frame(matrix(NA_real_, ncol = length(prob_placeholder_names), nrow = 1)), prob_placeholder_names)
      })

    } else {
      placeholder_lalonde_names <- c("est1_lalonde", "sd1_lalonde", "ci1_l_lalonde", "ci1_u_lalonde",
                                     "est2_lalonde", "sd2_lalonde", "ci2_l_lalonde", "ci2_u_lalonde",
                                     "est_compare_lalonde", "sd_compare_lalonde", "compare_ci_l_lalonde", "compare_ci_u_lalonde")
      lalonde_criteria_i <- setNames(data.frame(matrix(NA_real_, ncol = length(placeholder_lalonde_names), nrow = 1)), placeholder_lalonde_names)

      prob_placeholder_names <- c("pr_m", "pr_l", "pr_t")
      post_prob_i <- setNames(data.frame(matrix(NA_real_, ncol = length(prob_placeholder_names), nrow = 1)), prob_placeholder_names)
    }

    # 5. Add decision columns if requested
    decision_pr <- NA_character_
    decision_ci <- NA_character_

    if (Lalonde_decision && posterior_infer && !is.null(post_prob_i) && !is.null(lalonde_criteria_i)) {
      if (all(c("pr_t", "pr_l") %in% names(post_prob_i)) &&
          all(c("compare_ci_l_lalonde", "compare_ci_u_lalonde") %in% names(lalonde_criteria_i)) &&
          !anyNA(post_prob_i[, c("pr_t", "pr_l")]) &&
          !anyNA(lalonde_criteria_i[, c("compare_ci_l_lalonde", "compare_ci_u_lalonde")])) {

        cat("Start making decisions based on Lalonde framework.\n")

        # Probability Decision
        sum_probs <- post_prob_i$pr_m + post_prob_i$pr_l + post_prob_i$pr_t
        if (abs(sum_probs - 1) > .Machine$double.eps * 100) {
          warning(paste("Probabilities pr_m, pr_l, pr_t for nsim =", nrep_val, "do not sum to 1 (sum =", round(sum_probs, 3), "). Decision might be unreliable."))
        }

        decision_pr <- case_when(
          (post_prob_i$pr_t > fsr) & (post_prob_i$pr_l + post_prob_i$pr_t > 1 - fgr) ~ "go",
          (post_prob_i$pr_t > fsr) & (post_prob_i$pr_l + post_prob_i$pr_t <= 1 - fgr) ~ "consider",
          post_prob_i$pr_t <= fsr ~ "no-go",
          TRUE ~ NA_character_ # Should not happen if pr_t is not NA
        )

        # CI Decision
        ci_lower_bound <- lalonde_criteria_i$compare_ci_l_lalonde
        ci_upper_bound <- lalonde_criteria_i$compare_ci_u_lalonde

        decision_ci <- if (lrv < tv) { # "Greater is better" implied (lrv < tv)
          case_when(
            ci_upper_bound > tv & ci_lower_bound > lrv ~ "go",
            ci_upper_bound > tv & ci_lower_bound <= lrv ~ "consider",
            ci_upper_bound <= tv ~ "no-go",
            TRUE ~ NA_character_
          )
        } else { # "Less is better" implied (lrv > tv)
          case_when(
            ci_lower_bound < tv & ci_upper_bound < lrv ~ "go",
            ci_lower_bound < tv & ci_upper_bound >= lrv ~ "consider",
            ci_lower_bound >= tv ~ "no-go",
            TRUE ~ NA_character_
          )
        }

      } else if (Lalonde_decision && posterior_infer) {
        warning(paste("Required probability or CI quantile results for Lalonde decision are missing or NA for nsim =", nrep_val, ". Skipping decision for this simulation."))
        decision_pr <- NA_character_
        decision_ci <- NA_character_
      }
    }

    # 6. Combine all results for this simulation into a single row data.table
    # Ensure all components are data.frames or data.tables before cbind/rbindlist
    result_df_i <- data.table::as.data.table(cbind(
      t(as.data.frame(params_results)), # Transpose basic results
      as.data.frame(post_est_ci_i),
      as.data.frame(lalonde_criteria_i),
      as.data.frame(post_prob_i),
      data.frame(decision_pr = decision_pr, decision_ci = decision_ci) # Add decisions
    ))
    result_df_i
  }

  # Stop the cluster if it was created
  if (parallel_enabled && !is.null(cl)) {
    cat("Parallel computation finished. Stopping cluster.\n")
    stopCluster(cl)
    registerDoSEQ() # Register sequential backend after stopping cluster
  } else {
    cat("Sequential computation finished.\n")
  }

  # --- Post-processing ---
  # Check if results_list is empty or contains only NULLs
  if (is.null(results_list) || nrow(results_list) == 0) {
    warning("No successful simulation results were returned by the parallel loop.")

    post_params_df <- data.frame(matrix(ncol = length(unique(c(paste0("treatment.", names(posterior_results$post_t$post), "_post"), paste0("control.", names(posterior_results$post_c$post), "_post"), "treatment.w_prior", "control.w_prior"))), nrow = 0))
    colnames(post_params_df) <- unique(c(paste0("treatment.", names(posterior_results$post_t$post), "_post"), paste0("control.", names(posterior_results$post_c$post), "_post"), "treatment.w_prior", "control.w_prior"))
    post_params_df$nsim <- numeric()

    placeholder_95ci_names <- c("est1_95ci", "sd1_95ci", "ci1_l_95ci", "ci1_u_95ci",
                                "est2_95ci", "sd2_95ci", "ci2_l_95ci", "ci2_u_95ci",
                                "est_compare_95ci", "sd_compare_95ci", "compare_ci_l_95ci", "compare_ci_u_95ci")
    post_est_ci_df <- setNames(data.frame(matrix(NA_real_, ncol = length(placeholder_95ci_names), nrow = 0)), placeholder_95ci_names)
    post_est_ci_df$nsim <- numeric()

    placeholder_lalonde_names <- c("est1_lalonde", "sd1_lalonde", "ci1_l_lalonde", "ci1_u_lalonde",
                                   "est2_lalonde", "sd2_lalonde", "ci2_l_lalonde", "ci2_u_lalonde",
                                   "est_compare_lalonde", "sd_compare_lalonde", "compare_ci_l_lalonde", "compare_ci_u_lalonde")
    prob_placeholder_names <- c("pr_m", "pr_l", "pr_t")
    decision_names <- c("decision_pr", "decision_ci")

    post_inference_df <- setNames(data.frame(matrix(NA_real_, ncol = length(placeholder_lalonde_names) + length(prob_placeholder_names), nrow = 0)), c(placeholder_lalonde_names, prob_placeholder_names))
    post_inference_df[, decision_names] <- NA_character_
    post_inference_df$nsim <- numeric()

    # Order columns consistently (nsim first)
    post_params_df <- post_params_df %>% select(nsim, everything())
    post_est_ci_df <- post_est_ci_df %>% select(nsim, everything())
    post_inference_df <- post_inference_df %>% select(nsim, everything())

    return(list(post_params = post_params_df,
                post_est_ci = post_est_ci_df,
                post_inference = post_inference_df))

  }

  numeric_cols_after_bind <- names(results_list)[sapply(results_list, is.numeric)]

  char_cols <- names(results_list)[sapply(results_list, is.character)]
  decision_cols <- c("decision_pr", "decision_ci")

  cols_to_numeric <- setdiff(numeric_cols_after_bind, "nsim")
  cols_to_character <- intersect(char_cols, decision_cols)

  results_list <- results_list %>%
    mutate(across(all_of(cols_to_numeric), as.numeric)) %>%
    mutate(across(all_of(setdiff(char_cols, decision_cols)), as.numeric))

  post_params_df <- results_list %>%
    select(nsim, ends_with("_post"), ends_with(".w_prior")) %>%
    as.data.frame()

  post_est_ci_df <- results_list %>%
    select(nsim, ends_with("_95ci")) %>%
    as.data.frame()

  post_inference_df <- results_list %>%
    select(nsim, ends_with("_lalonde"), starts_with("pr_"), one_of(decision_cols)) %>%
    as.data.frame()

  cat("Function execution complete.\n")

  return(list(post_params = data.table::as.data.table(post_params_df),
              post_est_ci = data.table::as.data.table(post_est_ci_df),
              post_inference = data.table::as.data.table(post_inference_df)))
}


