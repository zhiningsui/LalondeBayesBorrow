rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)

# Simulation ------------------------------------------------------------

# Define historical control configurations
param_grid_h <- expand.grid(
  ctrl_h_n = 180,
  ctrl_h_mu = -0.2,
  ctrl_h_sigma = 0.4,
  stringsAsFactors = FALSE
)

# Define current trial configurations
param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_mu = seq(-0.4, 0, 0.05),
  trt_sigma = 0.4,
  ctrl_sigma = 0.4,
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_mu = ctrl_mu - 0.1)

param_grid2 <- expand.grid(
  trt_n = c(45),
  ctrl_mu = seq(-0.4, 0, 0.05),
  trt_sigma = 0.4,
  ctrl_sigma = 0.4,
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_mu = ctrl_mu - 0.2)

param_grid <- bind_rows(param_grid1, param_grid2)

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "continuous")

data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "continuous")


nsim = 10000
bayes_results <- list()
data_gen_params_h <- data_gen_params_list_h[[1]]

for (i in seq_along(data_gen_params_list)) {
  start_time_i <- Sys.time()

  data_gen_params <- data_gen_params_list[[i]]

  cat("Running", nsim, "simulations for data_gen_params set", i, "\n\n")

  for (arm in names(data_gen_params)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params[[arm]]$n, "\n")
    cat("\t  mu:", data_gen_params[[arm]]$mu, "\n")
    cat("\t  sigma:", data_gen_params[[arm]]$sigma, "\n\n")
  }
  for (arm in names(data_gen_params_h)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params_h[[arm]]$n, "\n")
    cat("\t  mu:", data_gen_params_h[[arm]]$mu, "\n")
    cat("\t  sigma:", data_gen_params_h[[arm]]$sigma, "\n\n")
  }

  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- lapply(data_gen_params, function(x) x$n)
  mu_list <- lapply(data_gen_params, function(x) x$mu)
  sigma_list <- lapply(data_gen_params, function(x) x$sigma)

  # --- Generate simulation dataset ---
  data <- data_gen_continuous(arm_names = arm_names,
                              nsim = nsim, n_list = n_list,
                              mu_list = mu_list, sigma_list = sigma_list,
                              seed = 123)

  # --- Obtain data summary for bayesian_lalonde_decision ---
  summary <- data$summary
  summary_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                          control_h.mu_hat = data_gen_params_h$control_h$mu,
                          control_h.s = data_gen_params_h$control_h$sigma / sqrt(data_gen_params_h$control_h$n))
  summary <- cbind(summary, summary_h)

  cat("proportion of |conflict| < 0.1 =", sum(abs(summary$control.mu_hat-summary$control_h.mu_hat) < 0.1)/nsim, "\n")
  cat("proportion of |conflict| < 0.15 =", sum(abs(summary$control.mu_hat-summary$control_h.mu_hat) < 0.15)/nsim, "\n")

  # --- Define prior parameters ---
  prior_params_list <- list()

  # No Borrowing ------------------------------------------------------------
  prior_params_list[["No Borrow"]] <- list(
    method = "SAM",
    treatment.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    treatment.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    treatment.w = 0,
    control.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    control.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    control.w = 0
  )

  # Modified SAM ------------------------------------------------------------
  prior_grid <- expand.grid(
    method = "SAM",
    treatment.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    treatment.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    treatment.w = 0,
    control.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    control.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    control.w = NA,
    control.delta_gate = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  ) %>%
    mutate(control.delta_SAM = control.delta_gate)

  prior_params_list2 <- lapply(rows_to_list(prior_grid), na_to_null)
  names(prior_params_list2) <- paste0("gated_SAM.", seq_along(prior_params_list2))

  # Naive SAM ---------------------------------------------------------------
  prior_grid <- expand.grid(
    method = "SAM",
    treatment.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    treatment.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    treatment.w = 0,
    control.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    control.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    control.w = NA,                      # allow SAM prior
    control.delta_SAM = c(0.1, 0.15),
    control.delta_gate = NA,
    control.ess_h = NA,
    stringsAsFactors = FALSE
  )
  prior_params_list3 <- lapply(rows_to_list(prior_grid), na_to_null)
  names(prior_params_list3) <- paste0("naive_SAM.", seq_along(prior_params_list3))


  prior_params_list <- c(prior_params_list, prior_params_list2, prior_params_list3)

  cat("Start Bayesian analysis.\n")
  for (k in seq_along(prior_params_list)) {
    start_time_k <- Sys.time()

    prior_params <- prior_params_list[[k]]

    settings <- c(list(true_value = data$true_value),
                  data_gen_params, data_gen_params_h,
                  prior_params)

    settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE)

    settings_num <- unlist(apply(settings, 2, num_or_null))
    settings <- as.data.frame(c(settings[!names(settings) %in% names(settings_num)], settings_num))

    post <- bayesian_lalonde_decision(endpoint = "continuous",
                                      data_summary = summary,
                                      settings = settings,
                                      prior_params = prior_params,
                                      arm_names = c(treatment = "treatment",
                                                    control = "control",
                                                    control_h = "control_h"),
                                      lrv = -0.1, tv = -0.2,
                                      fgr = 0.2, fsr = 0.1,
                                      posterior_infer = T,
                                      Lalonde_decision = T,
                                      verbose = F)

    post$metrics_post_dist <- cbind(settings,
                                    calc_post_dist_metrics(
                                      endpoint = "OR",
                                      true_value = data$true_value$compare_true,
                                      post_est_ci = post$post_est_ci))

    # --- Save elapsed time for this prior ---
    end_time_k <- Sys.time()
    elapsed_k  <- as.numeric(difftime(end_time_k, start_time_k, units = "secs"))
    post$elapsed_time_secs <- elapsed_k

    bayes_results[[length(bayes_results) + 1]] <- post

    cat("Time for Bayesian analysis with prior parameter list", k, "for data_gen_params set", i, "=",
        round(elapsed_k, 2), "seconds\n\n")
  }

  end_time_i <- Sys.time()
  elapsed_i  <- as.numeric(difftime(end_time_i, start_time_i, units = "secs"))
  cat("Total time for data_gen_params set", i, "=",
      round(elapsed_i, 2), "seconds\n\n")
}

saveRDS(bayes_results, file = "sim_DpR_bayes_results_0825.rds")
