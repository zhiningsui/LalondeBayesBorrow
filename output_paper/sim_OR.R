rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)

# Simulations -------------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 180,
  ctrl_h_p = 0.3,
  stringsAsFactors = FALSE
)

param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.1, 0.5, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_p = 0.1 + ctrl_p)

param_grid2 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.1, 0.5, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_p = 0.2 + ctrl_p)


param_grid <- bind_rows(param_grid1, param_grid2)

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "binary")

data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "binary")

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
    cat("\t  p:", data_gen_params[[arm]]$p, "\n")
  }
  for (arm in names(data_gen_params_h)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params_h[[arm]]$n, "\n")
    cat("\t  p:", data_gen_params_h[[arm]]$p, "\n")
  }

  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- sapply(data_gen_params, function(x) x$n)
  prob_list <- sapply(data_gen_params, function(x) x$p)

  # --- Generate simulation dataset ---
  data <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list, seed = 123)

  # --- Obtain data summary for bayesian_lalonde_decision ---
  summary <- data$summary
  summary_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                          control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p)
  summary <- cbind(summary, summary_h)

  # --- Define prior parameters ---
  prior_grid <- expand.grid(
    treatment.a0 = 1,
    treatment.b0 = 1,
    treatment.a = 1,
    treatment.b = 1,
    treatment.w = 0,
    control.a0 = 1,
    control.b0 = 1,
    control.a = 1,
    control.b = 1,
    control.w = NA, # allow SAM prior
    control.delta_gate = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  ) %>%
    mutate(control.delta_SAM = control.delta_gate)

  # Convert to list and replace NA with NULL
  prior_params_list <- apply(prior_grid, 1, function(row) {
    out <- as.list(row)
    if (is.na(out$control.w)) out$control.w <- NULL
    if (is.na(out$control.delta_gate)) out$control.delta_gate <- NULL
    out
  })

  # Add no borrowing prior
  prior_params_list[[length(prior_params_list)+1]] <- list(
    treatment.a0 = 1,
    treatment.b0 = 1,
    treatment.a = 1,
    treatment.b = 1,
    treatment.w = 0,
    control.a0 = 1,
    control.b0 = 1,
    control.a = 1,
    control.b = 1,
    control.w = 0
  )

  # --- Run analysis ---
  cat("Start Bayesian analysis.\n")
  for (k in seq_along(prior_params_list)) {
    start_time_k <- Sys.time()

    prior_params <- prior_params_list[[k]]

    settings <- c(list(true_value = data$true_value),
                  data_gen_params, data_gen_params_h,
                  prior_params)

    settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
      mutate(across(
        .cols = -c(ends_with(".name")),
        .fns = as.numeric
      ))

    post <- bayesian_lalonde_decision(endpoint = "binary",
                                      data_summary = summary,
                                      settings = settings,
                                      prior_params = prior_params,
                                      arm_names = c(treatment = "treatment",
                                                    control = "control",
                                                    control_h = "control_h"),
                                      lrv = 0.1, tv = 0.2,
                                      fgr = 0.2, fsr = 0.1,
                                      posterior_infer = T,
                                      Lalonde_decision = T,
                                      verbose = F)

    post$metrics_post_dist <- cbind(settings,
                                    calc_post_dist_metrics(endpoint = "OR",
                                                           true_value = post$post_est_ci$true_value.compare_true,
                                                           post_est_ci = post$post_est_ci))

    bayes_results[[length(bayes_results)+1]] <- post

    end_time_k <- Sys.time()
    cat("Time for Bayesian analysis with prior parameter list", k, "=", round(difftime(end_time_k, start_time_k, units = "secs"), 2), "seconds\n\n")
  }
  end_time_i <- Sys.time()
  cat("Total time for data_gen_params set", i, "=", round(difftime(end_time_i, start_time_i, units = "secs"), 2), "seconds\n\n")
}

saveRDS(bayes_results, file = "output_paper/sim_OR_bayes_results.rds")


