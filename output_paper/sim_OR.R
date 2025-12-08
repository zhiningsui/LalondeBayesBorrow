rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)

# ---------- helpers ---------------------------------------------------
mc_se <- function(p, n) sqrt(p * (1 - p) / max(n, 1))

# ---------- Core calibrator (from posterior outputs) ------------------------
## Calibrate (tau_TV, tau_LRV) for Lalonde framework using post_tv / post_lrv
# Implements Steps 1–3 (FSR, FGR, refinement, verification-ready).
#
# Arguments:
#   post_tv, post_lrv: outputs from bayesian_lalonde_decision(..., posterior_infer=TRUE)
#   alpha_N, alpha_G: target FSR and FGR
#   bisection_tol: numeric tolerance for bisection convergence
#   max_iter: max iterations for bisection
#   shrink_factor: multiplier (<1) used to slightly reduce tau0 if initial FGR > alpha_G
#   verbose: logical; print messages
#
# Returns: list of thresholds and diagnostics
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
  pti_tv <- post_tv$post_inference
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
  tau_TV <- as.numeric(quantile(pT_tv, probs = alpha_N, names = FALSE, na.rm = TRUE))
  FSR_hat <- mean(pT_tv < tau_TV)
  FSR_se <- mc_se(FSR_hat, n_tv)

  if (verbose) {
    message("-----------------------------------------------------------")
    message(" Step 1: Calibrating tau_TV to control False Stop Risk (FSR) ")
    message("-----------------------------------------------------------")
    message(sprintf(" - Monte Carlo samples: %d", n_tv))
    message(sprintf(" - Percentile rule: tau_TV = Quantile_{alpha_N=%.3f}(pT)", alpha_N))
    message(sprintf(" - Calibrated tau_TV = %.4f", tau_TV))
    message(sprintf(" - Estimated FSR = %.4f (SE = %.4f) [target alpha_N = %.3f]",
                    FSR_hat, FSR_se, alpha_N))
    message("  -> FSR is the probability of No-Go when Delta = Delta_TV (target effect).")
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
  if (qT == 0) { # No TV gate passes at LRV -> FGR = 0 for any tau_LRV; use unconditional safe bound
    tau_LRV_feasible <- as.numeric(quantile(pM_lrv, probs = alpha_G, names = FALSE, na.rm = TRUE))
    FGR_hat_fun <- function(tau) mean((pT_lrv >= tau_TV) & (pM_lrv < tau))
    tau_LRV_star <- tau_LRV_feasible
    FGR_hat <- 0
    FGR_se <- 0

    if (verbose) {
      message(" - q_T = 0 (no trial meets TV threshold under Delta = LRV).")
      message(" - FGR = 0 for all tau_LRV -> Using unconditional alpha_G-quantile.")
      message(sprintf(" - tau_LRV (unconditional) = %.4f", tau_LRV_feasible))
    }

  } else {
    cond_level <- min(1, alpha_G / qT)
    tau_LRV_feasible <- as.numeric(quantile(pM_lrv[T_idx], probs = cond_level,
                                            names = FALSE, na.rm = TRUE))
    # Empirical FGR function at Delta = LRV with Go event
    FGR_hat_fun <- function(tau) mean((pT_lrv >= tau_TV) & (pM_lrv < tau))

    if (verbose) {
      message(sprintf(" - Monte Carlo samples: %d", n_lrv))
      message(sprintf(" - Estimated q_T = Pr(p_T ≥ tau_TV | Delta = Delta_LRV) = %.4f", qT))
      message(sprintf(" - Conditional percentile level = alpha_G / q_T = %.4f", cond_level))
      message(sprintf(" - Feasible tau_LRV (conditional) = %.4f", tau_LRV_feasible))
    }

    # Guard for Monte Carlo overshoot
    tau0 <- tau_LRV_feasible
    overshoot_iter <- 0
    while (FGR_hat_fun(tau0) > alpha_G && tau0 > 0) {
      overshoot_iter <- overshoot_iter + 1
      tau0 <- tau0 * shrink_factor

      if (verbose) {
        message(sprintf("    + Iter %02d: FGR(tau=%.4f)=%.4f > alpha_G=%.3f -> shrinking tau x = %.2f -> %.4f",
                        overshoot_iter, tau0 / shrink_factor, FGR_hat_fun(tau0 / shrink_factor),
                        alpha_G, shrink_factor, tau0))
      }

      if (overshoot_iter > 20) {
        warning("  Overshoot correction exceeded 20 iterations. Check Monte Carlo stability.")
        break
      }
    }

    # Refinement: bisection to largest admissible tau_LRV
    if (verbose) message(" - Refining via bisection to largest admissible tau_LRV...")
    low <- tau0
    high <- 1.0
    for (iter in 1:max_iter) {
      if ((high - low) <= bisection_tol) break
      mid <- 0.5 * (low + high)
      if (FGR_hat_fun(mid) <= alpha_G) low <- mid else high <- mid
      if (verbose && iter %% 5 == 0) {
        message(sprintf("    + Iter %02d: interval [%.4f, %.4f] mid=%.4f FGR(mid)=%.4f",
                        iter, low, high, mid, FGR_hat_fun(mid)))
      }
    }
    tau_LRV_star <- low
    FGR_hat <- FGR_hat_fun(tau_LRV_star)
    FGR_se <- mc_se(FGR_hat, n_lrv)

    if (verbose) {
      message(sprintf(" - Final tau_LRV* = %.4f (largest admissible with FGR ≤ alpha_G)", tau_LRV_star))
      message(sprintf(" - Empirical FGR = %.4f (SE = %.4f) [target alpha_G = %.3f]",
                      FGR_hat, FGR_se, alpha_G))
      message("  -> FGR is the probability of Go when Delta = Delta_LRV (minimal effect).")
    }
  }

  # Summary and return
  if (verbose) {
    message("")
    message("===========================================================")
    message(" Calibration Summary")
    message("===========================================================")
    message(sprintf("  tau_TV = %.4f  (controls FSR ≤ alpha_N = %.3f)", tau_TV, alpha_N))
    message(sprintf("  tau_LRV = %.4f  (controls FGR ≤ alpha_G = %.3f)", tau_LRV_star, alpha_G))
    message(sprintf("  Achieved FSR = %.4f (SE = %.4f)", FSR_hat, FSR_se))
    message(sprintf("  Achieved FGR = %.4f (SE = %.4f)", FGR_hat, FGR_se))
    message("===========================================================")
  }

  list(
    tau_TV = tau_TV,
    tau_LRV = tau_LRV_star,
    tau_LRV_feasible = tau_LRV_feasible,
    alpha_N = alpha_N,
    alpha_G = alpha_G,
    FSR_hat = FSR_hat, FSR_se = FSR_se, n_tv = n_tv,
    FGR_hat = FGR_hat, FGR_se = FGR_se, n_lrv = n_lrv,
    qT = ifelse(is.finite(qT), qT, NA_real_),
    shrink_factor = shrink_factor
  )
}

# ---------- Targets and design constants -----------------------------------
alpha_G <- 0.20
alpha_N <- 0.10
lrv <- 0.10
tv <- 0.20

# ---------- Historical control (no prior–data conflict) --------------------
param_grid_h <- expand.grid(ctrl_h_n = 180, ctrl_h_p = 0.3, stringsAsFactors = FALSE)
data_gen_params_h <- create_data_gen_params(as.list(param_grid_h[1,]), endpoint = "binary")

# ---------- Build prior/method list ----------------------------------------

prior_params_list <- list()

# No Borrowing
prior_params_list[["NoBorrow"]] <- list(
  method = "SAM",
  treatment.a0 = 1, treatment.b0 = 1, treatment.a = 1, treatment.b = 1, treatment.w = 0,
  control.a0  = 1, control.b0  = 1, control.a  = 1, control.b  = 1, control.w  = 0
)

# Modified SAM (gated)
prior_grid <- expand.grid(
  method = "SAM",
  treatment.a0 = 1, treatment.b0 = 1, treatment.a = 1, treatment.b = 1, treatment.w = 0,
  control.a0 = 1, control.b0 = 1, control.a = 1, control.b = 1,
  control.w = NA,                      # allow SAM prior
  control.delta_gate = c(0.1, 0.15),
  control.ess_h = c(45, 90, 180),
  stringsAsFactors = FALSE
)
prior_grid$control.delta_SAM <- prior_grid$control.delta_gate
prior_params_list2 <- lapply(rows_to_list(prior_grid), na_to_null)
names(prior_params_list2) <- paste0("gated_SAM.", seq_along(prior_params_list2))

# Naive SAM
prior_grid <- expand.grid(
  method = "SAM",
  treatment.a0 = 1, treatment.b0 = 1, treatment.a = 1, treatment.b = 1, treatment.w = 0,
  control.a0 = 1, control.b0 = 1, control.a = 1, control.b = 1,
  control.w = NA,                      # allow SAM prior
  control.delta_SAM = c(0.1, 0.15),
  control.delta_gate = NA,
  control.ess_h = NA,
  stringsAsFactors = FALSE
)
prior_params_list3 <- lapply(rows_to_list(prior_grid), na_to_null)
names(prior_params_list3) <- paste0("naive_SAM.", seq_along(prior_params_list3))

# DPP (Empirical Bayes)
prior_grid <- expand.grid(
  method = "DPP",
  DPP.method = "Empirical Bayes",
  treatment.a = 1, treatment.b = 1,
  control.a = 1, control.b = 1,
  control.delta_gate = c(0.1, 0.15),
  control.ess_h = c(45, 90, 180),
  DPP.theta = NA, DPP.eta = NA,
  stringsAsFactors = FALSE
)
prior_params_list4 <- lapply(rows_to_list(prior_grid), na_to_null)
names(prior_params_list4) <- paste0("DPP.", seq_along(prior_params_list4))

prior_params_list <- c(prior_params_list, prior_params_list2, prior_params_list3, prior_params_list4)


# ---------- ONE-TIME CALIBRATION (NO conflict) -----------------------------
# Choose a single "no-conflict" study design for calibration:
n_t_cal <- 45
n_c_cal <- 45
p_c_cal <- 0.30
M_cal <- 50000  # large M for stable percentiles

## Under target case (Delta = Delta_TV)
data_tv_nc <- data_gen_binary(
  n_arms = 2, arm_names = c("treatment","control"), nsim = M_cal,
  n_list  = c(treatment = n_t_cal, control = n_c_cal),
  prob_list = c(treatment = p_c_cal + tv, control = p_c_cal),
  seed = 1001
)
summary_tv_nc <- data_tv_nc$summary %>%
  mutate(control_h.n = data_gen_params_h$control_h$n,
         control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p)

data_tv_nc2 <- data_gen_binary(
  n_arms = 2, arm_names = c("treatment","control"), nsim = 10000,
  n_list  = c(treatment = n_t_cal, control = n_c_cal),
  prob_list = c(treatment = p_c_cal + tv, control = p_c_cal),
  seed = 1001
)
summary_tv_nc2 <- data_tv_nc2$summary %>%
  mutate(control_h.n = data_gen_params_h$control_h$n,
         control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p)

summary_tv_nc1 <- summary_tv_nc %>%
  filter(nsim %in% 1:10000)


all.equal(summary_tv_nc1, summary_tv_nc2)


# Rows in summary_tv_nc1 but not in summary_tv_nc2
diff_1 <- anti_join(summary_tv_nc, summary_tv_nc2)

# Rows in df2 but not in df1
diff_2 <- anti_join(summary_tv_nc2, summary_tv_nc)

# All differing rows (from either side)
diff_all <- bind_rows(diff_1, diff_2)



## Under minimal case (Delta = Delta_LRV)
data_lrv_nc <- data_gen_binary(
  n_arms = 2, arm_names = c("treatment","control"), nsim = M_cal,
  n_list = c(treatment = n_t_cal, control = n_c_cal),
  prob_list = c(treatment = p_c_cal + lrv, control = p_c_cal),
  seed = 1002
)
summary_lrv_nc <- data_lrv_nc$summary %>%
  mutate(control_h.n = data_gen_params_h$control_h$n,
         control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p)

settings_nc <- data.frame(
  treatment.n = n_t_cal, control.n = n_c_cal,
  control_h.n = data_gen_params_h$control_h$n,
  control_h.p = data_gen_params_h$control_h$p,
  true_value  = data_tv_nc$true_value
)

# --- Calibrate tau's for EVERY method on SAME summaries (no conflict) ------

taus_by_method <- list()
calib_rows <- list()

cat("\n==== Calibrating thresholds ONCE under NO data–prior conflict ====\n")
for (k in seq_along(prior_params_list)) {
  prior_params <- prior_params_list[[k]]
  method_name <- names(prior_params_list)[k]

  # Posterior inference at Delta=TV and Delta=LRV on the SAME summaries

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

  cal <- calibrate_taus(post_tv = post_tv, post_lrv = post_lrv,
                        alpha_N = alpha_N, alpha_G = alpha_G,
                        bisection_tol = 1e-3, max_iter = 40,
                        shrink_factor = 0.98, verbose = TRUE)

  # flatten prior for reporting
  prior_flat <- tryCatch(as.data.frame(t(unlist(prior_params)), stringsAsFactors = FALSE),
                         error = function(e) data.frame())

  taus_by_method[[method_name]] <- list(tau_TV = cal$tau_TV, tau_LRV = cal$tau_LRV)

  calib_rows[[k]] <- cbind(
    data.frame(
      method = method_name,
      tau_TV = cal$tau_TV,
      tau_LRV = cal$tau_LRV,
      FSR_hat = cal$FSR_hat,
      FSR_se = cal$FSR_se,
      n_tv = cal$n_tv,
      FGR_hat = cal$FGR_hat,
      FGR_se = cal$FGR_se,
      n_lrv = cal$n_lrv,
      qT = cal$qT,
      alpha_N = alpha_N,
      alpha_G = alpha_G,
      stringsAsFactors = FALSE
    ),
    prior_flat
  )
}

calib_results <- bind_rows(calib_rows)
saveRDS(list(taus_by_method = taus_by_method,
             calib_table = calib_results,
             lrv = lrv, tv = tv,
             note = "Calibrated under NO data–prior conflict; reuse unchanged for all conflict levels."),
        file = "CALIBRATION_NO_CONFLICT.rds")
print(calib_results)




# RUN SIMULATIONS ---------------------------------------------------------

# ---- helpers ---------------------------------------------------------------
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}
make_key <- function(i, meth) sprintf("design%03d__%s", i, safe_name(meth))


prior_params_list <- list()

# No Borrowing
prior_params_list[["NoBorrow"]] <- list(
  method = "SAM",
  treatment.a0 = 1, treatment.b0 = 1, treatment.a = 1, treatment.b = 1, treatment.w = 0,
  control.a0  = 1, control.b0  = 1, control.a  = 1, control.b  = 1, control.w  = 0
)

# Modified SAM (gated)
prior_grid <- expand.grid(
  method = "SAM",
  treatment.a0 = 1, treatment.b0 = 1, treatment.a = 1, treatment.b = 1, treatment.w = 0,
  control.a0 = 1, control.b0 = 1, control.a = 1, control.b = 1,
  control.w = NA,                      # allow SAM prior
  control.delta_gate = c(0.1, 0.15),
  control.ess_h = c(45, 90, 180),
  stringsAsFactors = FALSE
)
prior_grid$control.delta_SAM <- prior_grid$control.delta_gate
prior_params_list2 <- lapply(rows_to_list(prior_grid), na_to_null)
names(prior_params_list2) <- paste0("gated_SAM.", seq_along(prior_params_list2))

# Naive SAM
prior_grid <- expand.grid(
  method = "SAM",
  treatment.a0 = 1, treatment.b0 = 1, treatment.a = 1, treatment.b = 1, treatment.w = 0,
  control.a0 = 1, control.b0 = 1, control.a = 1, control.b = 1,
  control.w = NA,                      # allow SAM prior
  control.delta_SAM = c(0.1, 0.15),
  control.delta_gate = NA,
  control.ess_h = NA,
  stringsAsFactors = FALSE
)
prior_params_list3 <- lapply(rows_to_list(prior_grid), na_to_null)
names(prior_params_list3) <- paste0("naive_SAM.", seq_along(prior_params_list3))

# DPP (Empirical Bayes)
prior_grid <- expand.grid(
  method = "DPP",
  DPP.method = "Empirical Bayes",
  treatment.a = 1, treatment.b = 1,
  control.a = 1, control.b = 1,
  control.delta_gate = c(0.1, 0.15),
  control.ess_h = c(45, 90, 180),
  DPP.theta = NA, DPP.eta = NA,
  stringsAsFactors = FALSE
)
prior_params_list4 <- lapply(rows_to_list(prior_grid), na_to_null)
names(prior_params_list4) <- paste0("DPP.", seq_along(prior_params_list4))

prior_params_list <- c(prior_params_list, prior_params_list2, prior_params_list3, prior_params_list4)

# ---------- Simulation design grid -----------------------------------------

param_grid_h <- expand.grid(ctrl_h_n = 180, ctrl_h_p = 0.3, stringsAsFactors = FALSE)
data_gen_params_h <- create_data_gen_params(as.list(param_grid_h[1,]), endpoint = "binary")

param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.1, 0.5, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_p = lrv + ctrl_p)

param_grid2 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.1, 0.5, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_p = tv + ctrl_p)


param_grid <- bind_rows(param_grid1, param_grid2)
data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "binary")

# ---------- Evaluate ALL scenarios with fixed taus -------------------------

CALIBRATION_NO_CONFLICT <- readRDS("CALIBRATION_NO_CONFLICT.rds")
list2env(CALIBRATION_NO_CONFLICT, .GlobalEnv)

settings_cal_tv <- data.frame(
  treatment.n = n_t_cal,
  control.n = n_c_cal,
  control_h.n = data_gen_params_h$control_h$n,
  treatment.p = p_c_cal + tv,
  control.p = p_c_cal,
  control_h.p = data_gen_params_h$control_h$p,
  seed = 1001
)

settings_cal_lrv <- data.frame(
  treatment.n = n_t_cal,
  control.n = n_c_cal,
  control_h.n = data_gen_params_h$control_h$n,
  treatment.p = p_c_cal + lrv,
  control.p = p_c_cal,
  control_h.p = data_gen_params_h$control_h$p,
  seed = 1002
)

nsim_eval <- 10000
bayes_results <- list()

cat("\n==== Applying calibrated taus to ALL conflict scenarios ====\n")
for (i in seq_along(data_gen_params_list)) {
  start_time_i <- Sys.time()

  data_gen_params <- data_gen_params_list[[i]]

  cat("\n===========================================================\n")
  cat("Design set", i, "\n")
  for (arm in names(data_gen_params)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params[[arm]]$n, "\n")
    cat("\t  p:", data_gen_params[[arm]]$p, "\n")
  }
  for (arm in names(data_gen_params_h)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params_h[[arm]]$n, "\n")
    cat("\t  p:", data_gen_params_h[[arm]]$p, "\n")
  }

  settings_eval1 <- data.frame(
    treatment.n = data_gen_params[["treatment"]]$n,
    control.n = data_gen_params[["control"]]$n,
    control_h.n = data_gen_params_h[["control_h"]]$n,
    treatment.p = data_gen_params[["treatment"]]$p,
    control.p = data_gen_params[["control"]]$p,
    control_h.p = data_gen_params_h[["control_h"]]$p
  )

  if (isTRUE(all.equal(settings_eval1,
                       settings_cal_tv[names(settings_cal_tv) != "seed"]))) {
    s1 = settings_cal_tv$seed
    cat("\n===========================================================\n",
        "Same design set as the calibration under target case. Use seed =", s1, "\n")
  } else if (isTRUE(all.equal(settings_eval1,
                              settings_cal_lrv[names(settings_cal_lrv) != "seed"]))) {
    s1 = settings_cal_lrv$seed
    cat("\n===========================================================\n",
        "Same design set as the calibration under minimal case. Use seed =", s1, "\n")
  } else {
    s1 = 1000
    cat("\n===========================================================\n",
        "Different design set as the calibration. Use seed =", s1, "\n")
  }

  # --- Generate an evaluation dataset -----------
  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- sapply(data_gen_params, function(x) x$n)
  prob_list <- sapply(data_gen_params, function(x) x$p)

  data_eval <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim_eval,
                               n_list = n_list, prob_list = prob_list, seed = s1)

  # Obtain data summary for bayesian_lalonde_decision
  summary_eval <- data_eval$summary %>%
    mutate(control_h.n = data_gen_params_h$control_h$n,
           control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p)

  cat("Start Bayesian evaluation analysis on this design...\n")
  for (k in seq_along(prior_params_list)) {
    prior <- prior_params_list[[k]]
    meth <- names(prior_params_list)[k]
    key <- make_key(i, meth)

    if (!is.null(names(bayes_results)) && key %in% names(bayes_results)) {
      cat("Skipping (already finished): ", key, "\n")
      next
    }

    start_time_k <- Sys.time()

    taus <- taus_by_method[[meth]]

    settings_eval <- c(list(true_value = data_eval$true_value),
                       data_gen_params, data_gen_params_h, taus, prior)
    settings_eval <- as.data.frame(t(unlist(settings_eval)), stringsAsFactors = FALSE)
    settings_num <- unlist(apply(settings_eval, 2, num_or_null))
    settings_eval <- as.data.frame(c(settings_eval[!names(settings_eval) %in% names(settings_num)], settings_num))

    post_eval <- bayesian_lalonde_decision(
      endpoint = "binary",
      data_summary = summary_eval,
      settings = settings_eval,
      prior_params = prior,
      arm_names = c(treatment="treatment", control="control", control_h="control_h"),
      lrv = lrv, tv = tv,
      fgr = taus$tau_LRV,
      fsr = taus$tau_TV,
      posterior_infer = TRUE,
      Lalonde_decision = TRUE,
      verbose = FALSE
    )

    post_eval$metrics_post_dist <- cbind(
      settings_eval,
      calc_post_dist_metrics(
        endpoint = "OR",
        true_value = unique(post_eval$post_est_ci$true_value.compare_true),
        post_est_ci = post_eval$post_est_ci
      )
    )
    # Save elapsed time for this prior ---
    end_time_k <- Sys.time()
    post_eval$elapsed_time_secs <- as.numeric(difftime(end_time_k, start_time_k, units = "secs"))

    bayes_results[[key]] <- post_eval
    cat("Completed Bayesian analysis: ", key, "\n\n")
  }

  end_time_i <- Sys.time()
  cat("Total time for design set", i, "=", round(as.numeric(difftime(end_time_i, start_time_i, units="secs")), 2), "seconds\n")
}

# ---------- Save global outputs --------------------------------------------
# saveRDS(ALL_CALIB_RESULTS, file = "ALL_CALIB_RESULTS.rds")
saveRDS(bayes_results, file = "sim_OR_bayes_results_1015.rds")

cat("\nAll done. Calibrated τ’s saved to ALL_CALIB_RESULTS.rds.\n")
cat("Remember: these τ’s were calibrated under NO data–prior conflict and should be APPLIED as-is across all conflict levels.\n")


# saveRDS(bayes_results, file = "sim_OR_bayes_results_0824.rds")

bayes_results1 <- readRDS("sim_OR_bayes_results_1010.rds")
bayes_results2 <- readRDS("sim_OR_bayes_results_1012.rds")

names(bayes_results2) %in% names(bayes_results1)

bayes_results <- bayes_results1
bayes_results[names(bayes_results2)] <- bayes_results2

saveRDS(bayes_results, file = "sim_OR_bayes_results_new.rds")




# # Calibration -------------------------------------------------------------
#
# nsim_cal <- 50000
#
# # Historical control (no prior-data conflict)
# param_h <- list(ctrl_h_n = 180, ctrl_h_p = 0.3)
# data_gen_params_h <- create_data_gen_params(param_h, endpoint = "binary")
#
# # Current trial design for calibration
# param_cal_tv <- list(trt_n = 45, ctrl_n = 45, ctrl_p = 0.3, trt_p = 0.3 + tv)
# data_gen_params_cal_tv <- create_data_gen_params(param_cal_tv, endpoint = "binary")
#
# param_cal_lrv <- list(trt_n = 45, ctrl_n = 45, ctrl_p = 0.3, trt_p = 0.3 + lrv)
# data_gen_params_cal_lrv <- create_data_gen_params(param_cal_lrv, endpoint = "binary")
#
# # Generate data summaries under calibration scenarios
# cat("Generating calibration datasets...\n")
#
# ## Under target case (Delta = Delta_TV)
# data_tv <- data_gen_binary(
#   n_arms = 2, arm_names = c("treatment","control"), nsim = nsim_cal,
#   n_list = c(treatment = data_gen_params_cal_tv$treatment$n, control = data_gen_params_cal_tv$control$n),
#   prob_list = c(treatment = data_gen_params_cal_tv$treatment$p, control = data_gen_params_cal_tv$control$p)
# )
# summary_tv <- data_tv$summary %>%
#   mutate(control_h.n = data_gen_params_h$control_h$n,
#          control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p)
#
# ## Under minimal case (Delta = Delta_LRV)
# data_lrv <- data_gen_binary(
#   n_arms = 2, arm_names = c("treatment","control"), nsim = nsim_cal,
#   n_list = c(treatment = data_gen_params_cal_lrv$treatment$n, control = data_gen_params_cal_lrv$control$n),
#   prob_list = c(treatment = data_gen_params_cal_lrv$treatment$p, control = data_gen_params_cal_lrv$control$p)
# )
# summary_lrv <- data_lrv$summary %>%
#   mutate(control_h.n = data_gen_params_h$control_h$n,
#          control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p)
#
# save(data_tv, summary_tv, data_lrv, summary_lrv, file = "sim_OR_calibration_datasets.RData")
#
#
# # --- Choose a representative prior (any method can be calibrated separately) ---
#
# param_grid_h <- expand.grid(
#   ctrl_h_n = 180,
#   ctrl_h_p = 0.3,
#   stringsAsFactors = FALSE
# )
#
# param_grid1 <- expand.grid(
#   trt_n = c(45),
#   ctrl_p = seq(0.1, 0.5, 0.05),
#   stringsAsFactors = FALSE
# ) %>%
#   mutate(ctrl_n = trt_n,
#          trt_p = 0.1 + ctrl_p)
#
# param_grid2 <- expand.grid(
#   trt_n = c(45),
#   ctrl_p = seq(0.1, 0.5, 0.05),
#   stringsAsFactors = FALSE
# ) %>%
#   mutate(ctrl_n = trt_n,
#          trt_p = 0.2 + ctrl_p)
#
#
# param_grid <- bind_rows(param_grid1, param_grid2)
#
# data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
#                                  create_data_gen_params, endpoint = "binary")
#
# data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
#                                create_data_gen_params, endpoint = "binary")
#
#
# alpha_G <- 0.2  # target FGR
# alpha_N <- 0.1  # target FSR
# lrv <- 0.10
# tv <- 0.20
#
# # Build a named list of priors/methods
# prior_list <- list(
#   "NoBorrow" = list(
#     method="SAM",
#     treatment.a0=1, treatment.b0=1, treatment.a=1, treatment.b=1, treatment.w=0,
#     control.a0=1,   control.b0=1,   control.a=1,   control.b=1,   control.w=0
#   ),
#   "SAM_gated_ess90" = list(
#     method="SAM",
#     treatment.a0=1, treatment.b0=1, treatment.a=1, treatment.b=1, treatment.w=0,
#     control.a0=1, control.b0=1, control.a=1, control.b=1,
#     control.w=NA, control.delta_gate=0.10, control.ess_h=90
#   ),
#   "DPP_EB_ess90" = list(
#     method="DPP",
#     DPP.method="Empirical Bayes",
#     treatment.a=1, treatment.b=1,
#     control.a=1,   control.b=1,
#     control.delta_gate=0.10,
#     control.ess_h=90
#   )
# )
#
#
#
#
#
#
# #
# # prior_params <- list(
# #   method = "SAM",
# #   treatment.a0 = 1, treatment.b0 = 1, treatment.a = 1, treatment.b = 1, treatment.w = 0,
# #   control.a0 = 1, control.b0 = 1, control.a = 1, control.b = 1,
# #   control.w = NA, control.delta_gate = 0.1, control.ess_h = 90
# # )
# #
# # settings <- data.frame(
# #   treatment.n = 45, control.n = 45,
# #   control_h.n = 180, control_h.p = 0.3,
# #   true_value = data_tv$true_value
# # )
# #
# # ## Step 1: Compute posterior probabilities at target case
# # post_tv <- bayesian_lalonde_decision(
# #   endpoint = "binary",
# #   data_summary = summary_tv,
# #   settings = settings,
# #   prior_params = prior_params,
# #   arm_names = c(treatment = "treatment",
# #                 control = "control",
# #                 control_h = "control_h"),
# #   lrv = lrv, tv = tv,
# #   posterior_infer = TRUE,
# #   Lalonde_decision = FALSE,
# #   verbose = FALSE
# # )
# #
# # post_lrv <- bayesian_lalonde_decision(
# #   endpoint = "binary",
# #   data_summary = summary_lrv,
# #   settings = settings,
# #   prior_params = prior_params,
# #   arm_names = c(treatment = "treatment",
# #                 control  = "control",
# #                 control_h = "control_h"),
# #   lrv = lrv, tv = tv,
# #   posterior_infer = TRUE,
# #   Lalonde_decision = FALSE,
# #   verbose = FALSE
# # )
# #
# # cal <- calibrate_taus(
# #   post_tv = post_tv,
# #   post_lrv = post_lrv,
# #   alpha_N = alpha_N,
# #   alpha_G = alpha_G,
# #   verbose = TRUE
# # )
#



