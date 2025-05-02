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
#' @importFrom data.table data.table
bayesian_lalonde_decision <- function(endpoint, data_summary,
                                      prior_params, lrv = 0, tv= 0, fgr = 0.2, fsr = 0.1, arm_names,
                                      posterior_infer = TRUE, Lalonde_decision = TRUE,
                                      EXP_TRANSFORM = FALSE) {
  nsims = unique(data_summary$nsim)

  # Detect if running on CRAN, Travis, GitHub Actions, or similar
  is_ci <- Sys.getenv("CI") == "true" ||
    Sys.getenv("TRAVIS") == "true" ||
    Sys.getenv("GITHUB_ACTIONS") == "true" ||
    Sys.getenv("_R_CHECK_LIMIT_CORES_") == "true"

  if (is_ci) {
    ncore <- 2  # safer for CI/CRAN
  } else {
    ncore <- max(1, detectCores() - 1)  # use more cores locally
  }

  cl <- makeCluster(ncore)
  required_vars <- c("endpoint", "data_summary", "lrv", "tv", "fgr", "fsr", "arm_names",
                     "prior_params")
  required_funcs <- c(get_all_functions(environment()), get_all_functions(globalenv()))
  clusterExport(cl, c(required_vars, required_funcs), envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(RBesT)
    library(data.table)
  })
  registerDoParallel(cl)

  post_params <- foreach(nrep = nsims, .combine=bind_rows) %dopar%{
    historical_ctrl <- if (any(names(arm_names) == "control_h")) data_summary[data_summary$nsim==nrep, startsWith(names(data_summary), "control_h.")] else NULL
    historical_trt <- if (any(names(arm_names) == "treatment_h")) data_summary[data_summary$nsim==nrep, startsWith(names(data_summary), "treatment_h.")]  else NULL
    current_ctrl <- data_summary[data_summary$nsim==nrep, startsWith(names(data_summary), "control.")]
    current_trt <- data_summary[data_summary$nsim==nrep, startsWith(names(data_summary), "treatment.")]

    if (!is.null(current_ctrl)) names(current_ctrl) <- sapply(strsplit(names(current_ctrl), "\\."), "[[", 2)
    if (!is.null(historical_ctrl)) names(historical_ctrl) <- sapply(strsplit(names(historical_ctrl), "\\."), "[[", 2)
    if (!is.null(current_trt)) names(current_trt) <- sapply(strsplit(names(current_trt), "\\."), "[[", 2)
    if (!is.null(historical_trt)) names(historical_trt) <- sapply(strsplit(names(historical_trt), "\\."), "[[", 2)

    post_c <-  posterior_distribution(endpoint,
                                      current = current_ctrl,
                                      historical = historical_ctrl,
                                      delta = prior_params[["control.delta"]],
                                      w = prior_params[["control.w"]],
                                      a = prior_params[["control.a"]],
                                      b = prior_params[["control.b"]],
                                      a0 = prior_params[["control.a0"]],
                                      b0 = prior_params[["control.b0"]],
                                      theta0 = prior_params[["control.theta0"]],
                                      s0 = prior_params[["control.s0"]],
                                      ess_h = prior_params[["control.ess_h"]])

    post_t <-  posterior_distribution(endpoint,
                                      current = current_trt,
                                      historical = historical_trt,
                                      delta = prior_params[["treatment.delta"]],
                                      w = prior_params[["treatment.w"]],
                                      a = prior_params[["treatment.a"]],
                                      b = prior_params[["treatment.b"]],
                                      a0 = prior_params[["treatment.a0"]],
                                      b0 = prior_params[["treatment.b0"]],
                                      theta0 = prior_params[["treatment.theta0"]],
                                      s0 = prior_params[["treatment.s0"]],
                                      ess_h = prior_params[["treatment.ess_h"]])

    post_params_i <- c(unlist(post_t$post), unlist(post_c$post))
    names(post_params_i) <- c(paste0("treatment.", names(post_t$post), "_post"), paste0("control.", names(post_c$post), "_post"))

    prior_ws <- c(post_t$w_prior, post_c$w_prior)
    names(prior_ws) <- c("treatment.w_prior", "control.w_prior")
    c(prior_ws, post_params_i)
  }

  stopCluster(cl)

  post_inference <- NULL
  post_prob <- NULL
  post_est_ci <- NULL
  lalonde_criteria <- NULL

  cat("Start evaluating the posterior distribution (bias, standard error, 95% coverage probability).\n")
  if (posterior_infer) {cat("Start obtaining posterior inference (credible intervals and posterior probabilities).\n")}
  for (nrep in nsims) {
    p <- post_params[nrep, ]

    post_t <- p[startsWith(names(p), "treatment.") & endsWith(names(p), "_post")]
    names(post_t) <- sapply(strsplit(gsub("_post", "", names(post_t)), "\\."), "[[", 2)
    post_t_mix <- convert_RBesT_mix(post = post_t, endpoint = endpoint)
    post_c <- p[startsWith(names(p), "control.") & endsWith(names(p), "_post")]
    names(post_c) <- sapply(strsplit(gsub("_post", "", names(post_c)), "\\."), "[[", 2)
    post_c_mix <- convert_RBesT_mix(post = post_c, endpoint = endpoint)
    post_est_ci <- rbind(post_est_ci, posterior_inference(post1 = post_t_mix,
                                                          post2 = post_c_mix,
                                                          quantiles = c(0.025, 0.975),
                                                          EXP_TRANSFORM))

    if (posterior_infer) {
      quantiles <- if (lrv < tv) c(low = fgr, upper = 1 - fsr) else c(low = fsr, upper = 1 - fgr)
      lalonde_criteria <- bind_rows(lalonde_criteria, posterior_inference(post1 = post_t_mix,
                                                                          post2 = post_c_mix,
                                                                          quantiles = quantiles,
                                                                          EXP_TRANSFORM))

      lrv_adj <- if (EXP_TRANSFORM) log(lrv) else lrv
      tv_adj <- if (EXP_TRANSFORM) log(tv) else tv

      pr_m <- posterior_prob(post1 = post_t_mix,
                             value = lrv_adj,
                             post2 = post_c_mix,
                             range_type = ifelse(lrv < tv, "less", "greater"))
      pr_l <- posterior_prob(post1 = post_t_mix,
                             value = ifelse(lrv < tv, lrv_adj, tv_adj),
                             post2 = post_c_mix, range_type = "between",
                             value2 = ifelse(lrv < tv, tv_adj, lrv_adj))
      pr_t <- posterior_prob(post1 = post_t_mix,
                             value = tv_adj,
                             post2 = post_c_mix, range_type = ifelse(lrv < tv, "greater", "less"))

      post_prob <- rbind(post_prob, c(pr_m, pr_l, pr_t))
      colnames(post_prob) <- c("pr_m", "pr_l", "pr_t")
    }

    post_inference <- cbind(lalonde_criteria, post_prob)
  }

  if (Lalonde_decision) {
    cat("Start making decisions based on Lalonde framework.\n")
    post_inference$decision_pr <- with(post_inference, ifelse( (pr_t > fsr) & ((pr_t + pr_l) > 1 - fgr), "go",
                                                              ifelse( (pr_t > fsr) & ((pr_t + pr_l) <= 1 - fgr), "consider",
                                                                     ifelse(pr_t <= fsr, "no-go", NA))))

    post_inference$decision_ci <- if (lrv < tv) {
      with(post_inference, ifelse(compare_ci_u > tv & compare_ci_l > lrv, "go",
                                  ifelse(compare_ci_u > tv & compare_ci_l <= lrv, "consider",
                                         ifelse(compare_ci_u <= tv, "no-go", NA))))

    } else {
      with(post_inference, ifelse(compare_ci_l < tv & compare_ci_u < lrv, "go",
                                  ifelse(compare_ci_l < tv & compare_ci_u >= lrv, "consider",
                                         ifelse(compare_ci_l >= tv, "no-go", NA))))
    }
    rownames(post_inference) <- NULL
  }

  return(list(post_params = post_params, post_est_ci = post_est_ci, post_inference = post_inference))
}



