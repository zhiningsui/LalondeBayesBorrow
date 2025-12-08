#' Calibrate Lalonde Go/No-Go Thresholds (one-time, no-conflict)
#'
#' Runs a one-time calibration of (\eqn{\tau_{TV}}, \eqn{\tau_{LRV}}) across all
#' methods in \code{prior_params_list}, using a *single* no-conflict design.
#' Internally:
#' \enumerate{
#'   \item Simulates data under \eqn{\Delta=\mathrm{TV}} and \eqn{\Delta=\mathrm{LRV}}
#'         for the same (t, c) design and a given historical control.
#'   \item Computes posterior probabilities with \code{bayesian_lalonde_decision(..., posterior_infer=TRUE)}.
#'   \item Calibrates \eqn{\tau_{TV}} (FSR control) and \eqn{\tau_{LRV}} (FGR control) via \code{calibrate_taus()}.
#' }
#'
#' @param prior_params_list A named \code{list} of prior/method specifications.
#'   Each element is passed to \code{bayesian_lalonde_decision()} as \code{prior_params}.
#' @param lrv,tv Numeric. Minimal and target effect margins (on the \code{endpoint="binary"} scale).
#' @param alpha_N,alpha_G Target risks for FSR and FGR, respectively.
#' @param treatment_n,control_n Integers. Sample sizes for treatment and control in the calibration design.
#' @param control_p Numeric in \[0,1\]. Control response rate in the calibration design.
#' @param control_h_n Integer. Historical control sample size.
#' @param control_h_p Numeric in \[0,1\]. Historical control response rate.
#' @param M_cal Integer. Number of Monte Carlo replicates for stable percentiles.
#' @param save_path Optional file path (e.g., \code{"CALIBRATION_NO_CONFLICT.rds"}) to
#'   save the returned object with \code{saveRDS}. If \code{NULL}, nothing is saved.
#' @param verbose Logical. Print progress.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{taus_by_method}: named list of \code{list(tau_TV, tau_LRV)} per method
#'   \item \code{calib_table}: data.frame of calibrated thresholds and diagnostics
#'   \item \code{lrv}, \code{tv}: the margins used
#'   \item \code{note}: a short note about the setup
#' }
#'
#' @details
#' This performs *one* calibration under no prior–data conflict and is intended to be
#' reused across other scenarios (e.g., with varying conflict levels) without re-calibration.
#'
#' @seealso \code{\link{bayesian_lalonde_decision}} (posterior inference) for binary endpoints,
#'   \code{\link{data_gen_binary}} (simulation), \code{\link{create_data_gen_params}} (historical setup).
#'
#' @export
#' @examples
#' \dontrun{
#' res <- lalonde_calibrate_no_conflict(
#'   prior_params_list = my_priors,
#'   lrv = 0.10, tv = 0.20,
#'   treatment_n = 45, control_n = 45, control_p = 0.30,
#'   control_h_n = 180, control_h_p = 0.30,
#'   M_cal = 5e4, save_path = "CALIBRATION_NO_CONFLICT.rds"
#' )
#' res$calib_table
#' }
lalonde_calibrate_no_conflict <- function(prior_params_list,
                                          lrv, tv,
                                          alpha_N = 0.10,
                                          alpha_G = 0.20,
                                          treatment_n = 45,
                                          control_n = 45,
                                          control_p = NULL,
                                          control_h_n = 180,
                                          control_h_p = 0.30,
                                          M_cal = 100000,
                                          save_path = NULL,
                                          verbose = TRUE) {

  if (is.null(control_p)) control_p <- control_h_p # Default = no conflict

  # ---- Historical control parameters ----
  data_gen_params_h <- create_data_gen_params(
    list(ctrl_h_n = control_h_n, ctrl_h_p = control_h_p),
    endpoint = "binary"
  )

  # ---- Simulate under Delta = TV ----
  if (verbose) message("Simulating calibration data under Delta = TV ...")
  data_tv_nc <- data_gen_binary(
    n_arms = 2, arm_names = c("treatment","control"), nsim = M_cal,
    n_list  = c(treatment = treatment_n, control = control_n),
    prob_list = c(treatment = control_p + tv, control = control_p)
  )

  summary_tv_nc <- data_tv_nc$summary %>%
    dplyr::mutate(
      control_h.n = data_gen_params_h$control_h$n,
      control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p
    )

  # ---- Simulate under Delta = LRV ----
  if (verbose) message("Simulating calibration data under Delta = LRV ...")
  data_lrv_nc <- data_gen_binary(
    n_arms = 2, arm_names = c("treatment","control"), nsim = M_cal,
    n_list  = c(treatment = treatment_n, control = control_n),
    prob_list = c(treatment = control_p + lrv, control = control_p)
  )

  summary_lrv_nc <- data_lrv_nc$summary %>%
    dplyr::mutate(
      control_h.n = data_gen_params_h$control_h$n,
      control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p
    )

  # ---- Settings row for bayesian_lalonde_decision() ----
  settings_nc <- data.frame(
    treatment.n = treatment_n,
    control.n  = control_n,
    control_h.n = data_gen_params_h$control_h$n,
    control_h.p = data_gen_params_h$control_h$p,
    # true_value.compare_true = tv, # not used by inference; included for traceability
    stringsAsFactors = FALSE
  )

  # ---- Loop over methods and calibrate ----
  taus_by_method <- list()
  calib_rows <- vector("list", length(prior_params_list))

  if (verbose) {
    message("\n==== Calibrating thresholds ONCE under NO data–prior conflict ====")
  }

  k <- 0
  for (method_name in names(prior_params_list)) {
    k <- k + 1
    prior_params <- prior_params_list[[method_name]]

    # Posterior inference at TV and LRV on the SAME summaries
    post_tv <- bayesian_lalonde_decision(
      endpoint = "binary",
      data_summary = summary_tv_nc,
      settings = settings_nc,
      prior_params = prior_params,
      arm_names = c(treatment="treatment", control="control", control_h="control_h"),
      lrv = lrv, tv = tv,
      posterior_infer = TRUE,
      Lalonde_decision = FALSE,
      verbose = FALSE
    )

    post_lrv <- bayesian_lalonde_decision(
      endpoint = "binary",
      data_summary = summary_lrv_nc,
      settings = settings_nc,
      prior_params = prior_params,
      arm_names = c(treatment="treatment", control="control", control_h="control_h"),
      lrv = lrv, tv = tv,
      posterior_infer = TRUE,
      Lalonde_decision = FALSE,
      verbose = FALSE
    )

    cal <- calibrate_taus(
      post_tv = post_tv, post_lrv = post_lrv,
      alpha_N = alpha_N, alpha_G = alpha_G,
      bisection_tol = 1e-3, max_iter = 40,
      shrink_factor = 0.98, verbose = verbose
    )

    taus_by_method[[method_name]] <- list(tau_TV = cal$tau_TV, tau_LRV = cal$tau_LRV)

    # Flatten prior for reportin
    prior_flat <- tryCatch(
      as.data.frame(t(unlist(prior_params)), stringsAsFactors = FALSE),
      error = function(e) data.frame()
    )

    calib_rows[[k]] <- cbind(
      data.frame(
        method = method_name,
        tau_TV = cal$tau_TV,
        tau_LRV = cal$tau_LRV,
        FSR_hat = cal$FSR_hat, FSR_se = cal$FSR_se, n_tv = cal$n_tv,
        FGR_hat = cal$FGR_hat, FGR_se = cal$FGR_se, n_lrv = cal$n_lrv,
        qT = cal$qT,
        alpha_N = alpha_N, alpha_G = alpha_G,
        stringsAsFactors = FALSE
      ),
      prior_flat
    )
  }

  calib_results <- dplyr::bind_rows(calib_rows)
  out <- list(
    taus_by_method = taus_by_method,
    calib_table = calib_results,
    lrv = lrv, tv = tv,
    note = "Calibrated under NO data–prior conflict; reuse unchanged for all conflict levels."
  )

  if (!is.null(save_path)) {
    saveRDS(out, file = save_path)
    if (verbose) message(sprintf("Saved calibration results to: %s", save_path))
  }

  out
}

#' Calibrator for (\eqn{\tau_{TV}}, \eqn{\tau_{LRV}})
#'
#' Given posterior-inference objects at \eqn{\Delta=TV} and \eqn{\Delta=LRV},
#' computes \eqn{\tau_{TV}} (FSR control) via percentile and the largest admissible
#' \eqn{\tau_{LRV}} (FGR control) via conditional percentile + bisection.
#'
#' Expects \code{post$post_inference} to contain columns for posterior probabilities of
#' target profile \code{pr_t} and minimal profile \code{pr_m}.
#'
#' @param post_tv,post_lrv Output objects from \code{bayesian_lalonde_decision(..., posterior_infer=TRUE)}.
#' @inheritParams lalonde_calibrate_no_conflict
#' @param bisection_tol,max_iter,shrink_factor Tuning knobs for the \eqn{\tau_{LRV}} search.
#' @param verbose Logical; print messages.
#' @return A list with thresholds and diagnostics used upstream by \code{lalonde_calibrate_no_conflict()}.
#' @export
calibrate_taus <- function(post_tv,
                           post_lrv,
                           alpha_N = 0.10,
                           alpha_G = 0.20,
                           bisection_tol = 1e-3,
                           max_iter = 40,
                           shrink_factor = 0.98,
                           verbose = TRUE) {
  # helpers
  pick <- function(df, choices) {
    nm <- intersect(choices, names(df))
    if (!length(nm)) stop("None of ", paste(choices, collapse = ", "), " found in: ",
                          paste(names(df), collapse = ", "))
    df[[nm[1]]]
  }
  mc_se <- function(p, n) sqrt(p * (1 - p) / max(n, 1))

  # extract posterior probabilities
  if (is.null(post_tv$post_inference) || is.null(post_lrv$post_inference)) {
    stop("post_tv$post_inference and post_lrv$post_inference must be present.")
  }
  pti_tv  <- post_tv$post_inference
  pti_lrv <- post_lrv$post_inference

  # P(Delta >= TV) at Delta = TV
  pT_tv <- pick(pti_tv,  c("pr_t", "pT", "pr_T", "post_pr_t"))
  pT_tv <- pT_tv[is.finite(pT_tv)]

  # P(Delta >= TV) and P(Delta <= LRV) at Delta = LRV
  pT_lrv <- pick(pti_lrv, c("pr_t", "pT", "pr_T", "post_pr_t"))
  pM_lrv <- pick(pti_lrv, c("pr_m", "pM", "pr_M", "post_pr_m"))
  ok <- is.finite(pT_lrv) & is.finite(pM_lrv)
  pT_lrv <- pT_lrv[ok]
  pM_lrv <- pM_lrv[ok]

  n_tv <- length(pT_tv)
  n_lrv <- length(pT_lrv)
  if (n_tv == 0 || n_lrv == 0) stop("Empty post_inference after NA filtering.")

  # Step 1: tau_TV by empirical percentile (FSR control)
  tau_TV <- as.numeric(stats::quantile(pT_tv, probs = alpha_N, names = FALSE, na.rm = TRUE))
  FSR_hat <- mean(pT_tv < tau_TV)
  FSR_se  <- mc_se(FSR_hat, n_tv)

  if (verbose) {
    message("-----------------------------------------------------------")
    message(" Step 1: Calibrating tau_TV to control False Stop Risk (FSR) ")
    message("-----------------------------------------------------------")
    message(sprintf(" - Monte Carlo samples: %d", n_tv))
    message(sprintf(" - Percentile rule: tau_TV = Quantile_{alpha_N=%.3f}(pT)", alpha_N))
    message(sprintf(" - Calibrated tau_TV = %.4f", tau_TV))
    message(sprintf(" - Estimated FSR = %.4f (SE = %.4f) [target alpha_N = %.3f]",
                    FSR_hat, FSR_se, alpha_N))
  }

  # Step 2: feasible tau_LRV by conditional percentile (FGR control)
  if (verbose) {
    message("")
    message("-----------------------------------------------------------")
    message(" Step 2: Calibrating tau_LRV to control False Go Risk (FGR) ")
    message("-----------------------------------------------------------")
  }

  T_idx <- (pT_lrv >= tau_TV)        # T = {p_T >= tau_TV} at Delta = LRV
  qT <- mean(T_idx)

  if (qT == 0) {
    tau_LRV_feasible <- as.numeric(stats::quantile(pM_lrv, probs = alpha_G, names = FALSE, na.rm = TRUE))
    FGR_hat_fun <- function(tau) mean((pT_lrv >= tau_TV) & (pM_lrv < tau))
    tau_LRV_star <- tau_LRV_feasible
    FGR_hat <- 0
    FGR_se  <- 0

    if (verbose) {
      message(" - q_T = 0 (no trial meets TV threshold under Delta = LRV).")
      message(" - FGR = 0 for all tau_LRV -> Using unconditional alpha_G-quantile.")
      message(sprintf(" - tau_LRV (unconditional) = %.4f", tau_LRV_feasible))
    }

  } else {
    cond_level <- min(1, alpha_G / qT)
    tau_LRV_feasible <- as.numeric(stats::quantile(pM_lrv[T_idx], probs = cond_level,
                                                   names = FALSE, na.rm = TRUE))
    FGR_hat_fun <- function(tau) mean((pT_lrv >= tau_TV) & (pM_lrv < tau))

    if (verbose) {
      message(sprintf(" - Monte Carlo samples (LRV): %d", n_lrv))
      message(sprintf(" - Estimated q_T = Pr(p_T ≥ tau_TV | Delta=LRV) = %.4f", qT))
      message(sprintf(" - Conditional percentile level = alpha_G / q_T = %.4f", cond_level))
      message(sprintf(" - Feasible tau_LRV (conditional) = %.4f", tau_LRV_feasible))
    }

    # Guard for Monte Carlo overshoot
    tau0 <- tau_LRV_feasible
    overshoot_iter <- 0L
    while (FGR_hat_fun(tau0) > alpha_G && tau0 > 0) {
      overshoot_iter <- overshoot_iter + 1L
      tau0 <- tau0 * shrink_factor

      if (verbose) {
        message(sprintf("    + Iter %02d: shrink tau (%.4f -> %.4f), FGR(prev)=%.4f",
                        overshoot_iter, tau0 / shrink_factor, tau0, FGR_hat_fun(tau0 / shrink_factor)))
      }
      if (overshoot_iter > 20L) {
        warning("  Overshoot correction exceeded 20 iterations. Check Monte Carlo stability.")
        break
      }
    }

    # Refinement: bisection to largest admissible tau_LRV
    if (verbose) message(" - Refining via bisection to largest admissible tau_LRV...")
    low <- tau0
    high <- 1.0
    for (iter in seq_len(max_iter)) {
      if ((high - low) <= bisection_tol) break
      mid <- 0.5 * (low + high)
      if (FGR_hat_fun(mid) <= alpha_G) low <- mid else high <- mid
      if (verbose && (iter %% 5L == 0L)) {
        message(sprintf("    + Iter %02d: [low, high]=[%.4f, %.4f], mid=%.4f, FGR(mid)=%.4f",
                        iter, low, high, mid, FGR_hat_fun(mid)))
      }
    }
    tau_LRV_star <- low
    FGR_hat <- FGR_hat_fun(tau_LRV_star)
    FGR_se  <- mc_se(FGR_hat, n_lrv)

    if (verbose) {
      message(sprintf(" - Final tau_LRV* = %.4f (largest admissible with FGR ≤ alpha_G)", tau_LRV_star))
      message(sprintf(" - Empirical FGR = %.4f (SE = %.4f) [target alpha_G = %.3f]",
                      FGR_hat, FGR_se, alpha_G))
    }
  }

  if (verbose) {
    message("")
    message("===========================================================")
    message(" Calibration Summary")
    message("===========================================================")
    message(sprintf("  tau_TV = %.4f  (controls FSR ≤ alpha_N = %.3f)", tau_TV, alpha_N))
    message(sprintf("  tau_LRV = %.4f (controls FGR ≤ alpha_G = %.3f)", tau_LRV_star, alpha_G))
    message(sprintf("  Achieved FSR = %.4f (SE = %.4f)", FSR_hat, FSR_se))
    message(sprintf("  Achieved FGR = %.4f (SE = %.4f)", FGR_hat, FGR_se))
    message("===========================================================")
  }

  list(
    tau_TV = tau_TV,
    tau_LRV = tau_LRV_star,
    tau_LRV_feasible = if (exists("tau_LRV_feasible")) tau_LRV_feasible else NA_real_,
    alpha_N = alpha_N, alpha_G = alpha_G,
    FSR_hat = FSR_hat, FSR_se = FSR_se, n_tv = n_tv,
    FGR_hat = FGR_hat, FGR_se = FGR_se, n_lrv = n_lrv,
    qT = ifelse(is.finite(qT), qT, NA_real_),
    shrink_factor = shrink_factor
  )
}
