rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)


# Simulations I -------------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 180,
  ctrl_h_p = 0.3,
  stringsAsFactors = FALSE
)

param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.15, 0.45, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_p = 0.1 + ctrl_p)

param_grid2 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.15, 0.45, 0.05),
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

  # Generate simulation dataset
  data <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list, seed = 123)

  summary <- data$summary
  summary_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                          control_h.count = data_gen_params_h$control_h$n*data_gen_params_h$control_h$p)
  summary <- cbind(summary, summary_h)

  prior_grid1 <- expand.grid(
    treatment.a0 = 1,
    treatment.b0 = 1,
    treatment.a = 1,
    treatment.b = 1,
    treatment.w = 0,
    control.a0 = 1,
    control.b0 = 1,
    control.a = 1,
    control.b = 1,
    control.w = c(NA, 0),
    control.delta_SAM = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  )

  prior_grid2 <- expand.grid(
    treatment.a0 = 1,
    treatment.b0 = 1,
    treatment.a = 1,
    treatment.b = 1,
    treatment.w = 0,
    control.a0 = 1,
    control.b0 = 1,
    control.a = 1,
    control.b = 1,
    control.w = c(NA, 0),
    control.delta_gate = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  ) %>%
    mutate(control.delta_SAM = control.delta_gate)

  prior_grid <- bind_rows(prior_grid1, prior_grid2)

  # Convert to list and replace NA with NULL
  prior_params_list <- apply(prior_grid, 1, function(row) {
    out <- as.list(row)
    if (is.na(out$control.w)) out$control.w <- NULL
    if (is.na(out$control.delta_gate)) out$control.delta_gate <- NULL
    out
  })

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
                                      prior_params = prior_params,
                                      arm_names = c(treatment = "treatment",
                                                       control = "control",
                                                       control_h = "control_h"),
                                      lrv = 0.1, tv = 0.2,
                                      fgr = 0.2, fsr = 0.1,
                                      posterior_infer = T,
                                      Lalonde_decision = T)

    post$metrics_post_dist <- calc_post_dist_metrics(endpoint = "OR",
                                                     true_value = data$true_value$compare_true,
                                                     post$post_est_ci)
    post$settings <- settings
    post$data_summary <- cbind(i, summary)

    bayes_results[[length(bayes_results)+1]] <- post

    end_time_k <- Sys.time()
    cat("Time for Bayesian analysis with prior parameter list", k, "=", round(difftime(end_time_k, start_time_k, units = "secs"), 2), "seconds\n\n")
  }
  end_time_i <- Sys.time()
  cat("Total time for data_gen_params set", i, "=", round(difftime(end_time_i, start_time_i, units = "secs"), 2), "seconds\n\n")
}

# old <- readRDS("simulations/sim_OR_bayes_results_5_19.rds")
# bayes_results <- c(old, bayes_results)
saveRDS(bayes_results, file = "simulations/sim_OR_bayes_results_5_19.rds")


post_params_all <- data.frame()
post_est_ci_all <- data.frame()
post_inference_all <- data.frame()
metrics_post_dist_all <- data.frame()
oc_all <- data.frame()
data_summary_all <- data.frame()
for (i in seq_along(bayes_results)) {
  res <- bayes_results[[i]]
  setting <- res$settings
  post_params <- cbind(setting, res$post_params)
  post_est_ci <- cbind(setting, res$post_est_ci)
  post_inference <- cbind(setting, res$post_inference)
  metrics_post_dist <- cbind(setting, res$metrics_post_dist)

  post_inference$decision_pr <- ifelse(is.na(post_inference$decision_pr), post_inference$decision_ci, post_inference$decision_pr)
  post_inference$decision_ci <- ifelse(is.na(post_inference$decision_ci), post_inference$decision_pr, post_inference$decision_ci)

  oc <- cbind(setting, obtain_oc(post_inference))

  data_summary_all <- bind_rows(data_summary_all, res$data_summary)
  post_params_all <- bind_rows(post_params_all, cbind(i, post_params))
  post_est_ci_all <- bind_rows(post_est_ci_all, cbind(i, post_est_ci))
  post_inference_all <- bind_rows(post_inference_all, cbind(i, post_inference))
  metrics_post_dist_all <- bind_rows(metrics_post_dist_all, cbind(i, metrics_post_dist))
  oc_all <- bind_rows(oc_all, cbind(i, oc))
}

bayes_results_all <- list(
  data_summary_all = data_summary_all,
  metrics_post_dist_all = metrics_post_dist_all,
  post_params_all = post_params_all,
  post_est_ci_all = post_est_ci_all,
  post_inference_all = post_inference_all,
  oc_all = oc_all
)

saveRDS(bayes_results_all, "simulations/sim_OR_bayes_results_df_5_19.rds")



# Simulations II -------------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 180,
  ctrl_h_p = 0.3,
  stringsAsFactors = FALSE
)

param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.15, 0.45, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_p = 0.1 + ctrl_p)

param_grid2 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.15, 0.45, 0.05),
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

  # Generate simulation dataset
  data <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list, seed = 123)

  summary <- data$summary
  summary_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                          control_h.count = data_gen_params_h$control_h$n*data_gen_params_h$control_h$p)
  summary <- cbind(summary, summary_h)

  prior_grid1 <- expand.grid(
    treatment.a0 = 0.1,
    treatment.b0 = 0.1,
    treatment.a = 0.1,
    treatment.b = 0.1,
    treatment.w = 0,
    control.a0 = 0.1,
    control.b0 = 0.1,
    control.a = 0.1,
    control.b = 0.1,
    control.w = c(NA, 0),
    control.delta_SAM = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  )

  prior_grid2 <- expand.grid(
    treatment.a0 = 0.1,
    treatment.b0 = 0.1,
    treatment.a = 0.1,
    treatment.b = 0.1,
    treatment.w = 0,
    control.a0 = 0.1,
    control.b0 = 0.1,
    control.a = 0.1,
    control.b = 0.1,
    control.w = c(NA, 0),
    control.delta_gate = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  ) %>%
    mutate(control.delta_SAM = control.delta_gate)

  prior_grid <- bind_rows(prior_grid1, prior_grid2)

  # Convert to list and replace NA with NULL
  prior_params_list <- apply(prior_grid, 1, function(row) {
    out <- as.list(row)
    if (is.na(out$control.w)) out$control.w <- NULL
    if (is.na(out$control.delta_gate)) out$control.delta_gate <- NULL
    out
  })

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
                                      prior_params = prior_params,
                                      arm_names = c(treatment = "treatment",
                                                    control = "control",
                                                    control_h = "control_h"),
                                      lrv = 0.1, tv = 0.2,
                                      fgr = 0.2, fsr = 0.1,
                                      posterior_infer = T,
                                      Lalonde_decision = T)

    post$metrics_post_dist <- calc_post_dist_metrics(endpoint = "OR",
                                                     true_value = data$true_value$compare_true,
                                                     post$post_est_ci)
    post$settings <- settings
    post$data_summary <- cbind(i, summary)

    bayes_results[[length(bayes_results)+1]] <- post

    end_time_k <- Sys.time()
    cat("Time for Bayesian analysis with prior parameter list", k, "=", round(difftime(end_time_k, start_time_k, units = "secs"), 2), "seconds\n\n")
  }
  end_time_i <- Sys.time()
  cat("Total time for data_gen_params set", i, "=", round(difftime(end_time_i, start_time_i, units = "secs"), 2), "seconds\n\n")
}

# old <- readRDS("simulations/sim_OR_bayes_results_5_19.rds")
# bayes_results <- c(old, bayes_results)
saveRDS(bayes_results, file = "simulations/sim_OR_bayes_results_5_19_beta01.rds")


post_params_all <- data.frame()
post_est_ci_all <- data.frame()
post_inference_all <- data.frame()
metrics_post_dist_all <- data.frame()
oc_all <- data.frame()
data_summary_all <- data.frame()
for (i in seq_along(bayes_results)) {
  res <- bayes_results[[i]]
  setting <- res$settings
  post_params <- cbind(setting, res$post_params)
  post_est_ci <- cbind(setting, res$post_est_ci)
  post_inference <- cbind(setting, res$post_inference)
  metrics_post_dist <- cbind(setting, res$metrics_post_dist)

  post_inference$decision_pr <- ifelse(is.na(post_inference$decision_pr), post_inference$decision_ci, post_inference$decision_pr)
  post_inference$decision_ci <- ifelse(is.na(post_inference$decision_ci), post_inference$decision_pr, post_inference$decision_ci)

  oc <- cbind(setting, obtain_oc(post_inference))

  data_summary_all <- bind_rows(data_summary_all, res$data_summary)
  post_params_all <- bind_rows(post_params_all, cbind(i, post_params))
  post_est_ci_all <- bind_rows(post_est_ci_all, cbind(i, post_est_ci))
  post_inference_all <- bind_rows(post_inference_all, cbind(i, post_inference))
  metrics_post_dist_all <- bind_rows(metrics_post_dist_all, cbind(i, metrics_post_dist))
  oc_all <- bind_rows(oc_all, cbind(i, oc))
}

bayes_results_all <- list(
  data_summary_all = data_summary_all,
  metrics_post_dist_all = metrics_post_dist_all,
  post_params_all = post_params_all,
  post_est_ci_all = post_est_ci_all,
  post_inference_all = post_inference_all,
  oc_all = oc_all
)

saveRDS(bayes_results_all, "simulations/sim_OR_bayes_results_df_5_19_beta01.rds")



# Simulations III -------------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 180,
  ctrl_h_p = 0.3,
  stringsAsFactors = FALSE
)

param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.15, 0.45, 0.05),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_p = 0.1 + ctrl_p)

param_grid2 <- expand.grid(
  trt_n = c(45),
  ctrl_p = seq(0.15, 0.45, 0.05),
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

  # Generate simulation dataset
  data <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list, seed = 123)

  summary <- data$summary
  summary_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                          control_h.count = data_gen_params_h$control_h$n*data_gen_params_h$control_h$p)
  summary <- cbind(summary, summary_h)

  prior_grid1 <- expand.grid(
    treatment.a0 = 0.01,
    treatment.b0 = 0.01,
    treatment.a = 0.01,
    treatment.b = 0.01,
    treatment.w = 0,
    control.a0 = 0.01,
    control.b0 = 0.01,
    control.a = 0.01,
    control.b = 0.01,
    control.w = c(NA, 0),
    control.delta_SAM = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  )

  prior_grid2 <- expand.grid(
    treatment.a0 = 0.01,
    treatment.b0 = 0.01,
    treatment.a = 0.01,
    treatment.b = 0.01,
    treatment.w = 0,
    control.a0 = 0.01,
    control.b0 = 0.01,
    control.a = 0.01,
    control.b = 0.01,
    control.w = c(NA, 0),
    control.delta_gate = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  ) %>%
    mutate(control.delta_SAM = control.delta_gate)

  prior_grid <- bind_rows(prior_grid1, prior_grid2)

  # Convert to list and replace NA with NULL
  prior_params_list <- apply(prior_grid, 1, function(row) {
    out <- as.list(row)
    if (is.na(out$control.w)) out$control.w <- NULL
    if (is.na(out$control.delta_gate)) out$control.delta_gate <- NULL
    out
  })

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
                                      prior_params = prior_params,
                                      arm_names = c(treatment = "treatment",
                                                    control = "control",
                                                    control_h = "control_h"),
                                      lrv = 0.1, tv = 0.2,
                                      fgr = 0.2, fsr = 0.1,
                                      posterior_infer = T,
                                      Lalonde_decision = T)

    post$metrics_post_dist <- calc_post_dist_metrics(endpoint = "OR",
                                                     true_value = data$true_value$compare_true,
                                                     post$post_est_ci)
    post$settings <- settings
    post$data_summary <- cbind(i, summary)

    bayes_results[[length(bayes_results)+1]] <- post

    end_time_k <- Sys.time()
    cat("Time for Bayesian analysis with prior parameter list", k, "=", round(difftime(end_time_k, start_time_k, units = "secs"), 2), "seconds\n\n")
  }
  end_time_i <- Sys.time()
  cat("Total time for data_gen_params set", i, "=", round(difftime(end_time_i, start_time_i, units = "secs"), 2), "seconds\n\n")
}

# old <- readRDS("simulations/sim_OR_bayes_results_5_19.rds")
# bayes_results <- c(old, bayes_results)
saveRDS(bayes_results, file = "simulations/sim_OR_bayes_results_5_19_beta001.rds")


post_params_all <- data.frame()
post_est_ci_all <- data.frame()
post_inference_all <- data.frame()
metrics_post_dist_all <- data.frame()
oc_all <- data.frame()
data_summary_all <- data.frame()
for (i in seq_along(bayes_results)) {
  res <- bayes_results[[i]]
  setting <- res$settings
  post_params <- cbind(setting, res$post_params)
  post_est_ci <- cbind(setting, res$post_est_ci)
  post_inference <- cbind(setting, res$post_inference)
  metrics_post_dist <- cbind(setting, res$metrics_post_dist)

  post_inference$decision_pr <- ifelse(is.na(post_inference$decision_pr), post_inference$decision_ci, post_inference$decision_pr)
  post_inference$decision_ci <- ifelse(is.na(post_inference$decision_ci), post_inference$decision_pr, post_inference$decision_ci)

  oc <- cbind(setting, obtain_oc(post_inference))

  data_summary_all <- bind_rows(data_summary_all, res$data_summary)
  post_params_all <- bind_rows(post_params_all, cbind(i, post_params))
  post_est_ci_all <- bind_rows(post_est_ci_all, cbind(i, post_est_ci))
  post_inference_all <- bind_rows(post_inference_all, cbind(i, post_inference))
  metrics_post_dist_all <- bind_rows(metrics_post_dist_all, cbind(i, metrics_post_dist))
  oc_all <- bind_rows(oc_all, cbind(i, oc))
}

bayes_results_all <- list(
  data_summary_all = data_summary_all,
  metrics_post_dist_all = metrics_post_dist_all,
  post_params_all = post_params_all,
  post_est_ci_all = post_est_ci_all,
  post_inference_all = post_inference_all,
  oc_all = oc_all
)

saveRDS(bayes_results_all, "simulations/sim_OR_bayes_results_df_5_19_beta001.rds")



# Analysis I ----------------------------------------------------------------

bayes_results_all <- readRDS("simulations/sim_OR_bayes_results_df_5_19.rds")

# + Check posterior weights -----------------------------------------------

post_params_all <- bayes_results_all$post_params_all %>%
  mutate(across(
    .cols = -c(ends_with(".name")),
    .fns = as.numeric
  )) %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  filter(borrowing == "Borrowing: Yes") %>%
  select("nsim", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p.diff",
         "control.w_prior", "control.w_post") %>%
  mutate(control.p.diff = factor(control.p.diff,
                                 levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15")),
         control.delta_SAM = factor(control.delta_SAM,
                                     levels = c("0.1", "0.15"),
                                     labels = c(expression(paste(delta, " = 0.1")),
                                                expression(paste(delta, " = 0.15")))),
         control.ess_h = factor(control.ess_h,
                                levels = c("45", "90", "180"),
                                labels = c(expression(n[ch*","*e] == 45),
                                           expression(n[ch*","*e] == 90),
                                           expression(n[ch*","*e] == 180))))

shade_df <- data.frame(
  control.delta_SAM = factor(c(0.1, 0.15), levels = c("0.1", "0.15"),
                              labels = c(expression(paste(delta, " = 0.1")),
                                         expression(paste(delta, " = 0.15")))),
  xmin = c(2, 1),
  xmax = c(6, 7),
  ymin = -Inf,
  ymax = Inf)

p1 <- ggplot(post_params_all, aes(x = control.p.diff, y = control.w_prior)) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  facet_grid(control.delta_SAM ~ control.ess_h + gate,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(
    x = expression(paste(theta[c], " - ", theta[ch])),
    y = "Prior Weight",
    title = "Prior Weight of Informative Component of SAM Prior"
  ) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10))

p2 <- ggplot(post_params_all, aes(x = control.p.diff, y = control.w_post)) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  facet_grid(control.delta_SAM ~ control.ess_h + gate,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(x = expression(paste(theta[c], " - ", theta[ch])),
       y = "Posterior Weight",
       title = "Posterior Weight of Informative Component of SAM Prior") +
  ylim(c(0,1)) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10))

library(patchwork)
p1 / p2
ggsave(filename = "simulations/sim_OR_SAM_weight_5_19.jpg", p1 / p2,
       width = 13, height = 8)


# + Generate summary table for paper --------------------------------------

post_inference_all <-  bayes_results_all$post_inference_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  select("nsim", "borrowing", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p",
         "est2_lalonde") %>%
  pivot_wider(names_from = borrowing, values_from = est2_lalonde) %>%
  mutate(pmd = `Borrowing: Yes` - `Borrowing: No`) %>%
  select(-c(`Borrowing: Yes`, `Borrowing: No`, true_value.compare_true))

# Summary statistics of PMD on control
metrics_tab <- post_inference_all %>%
  group_by(gate, control.delta_SAM, control.ess_h, control.p) %>%
  summarize(mean_pmd = mean(pmd),
            sd_pmd = sd(pmd),
            .groups = "drop")

# risks
risk_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  filter(true_value.compare_true %in% c(0.1, 0.2),
         decision_pr %in% c("go", "no-go")) %>%
  mutate(
    CorrectLabel = case_when(
      true_value.compare_true == 0.2 & decision_pr == "go" ~ "CGR",
      true_value.compare_true == 0.1 & decision_pr == "no-go" ~ "CSR",
      TRUE ~ NA_character_
    ),
    RiskLabel = case_when(
      true_value.compare_true == 0.1 & decision_pr == "go" ~ "FGR",
      true_value.compare_true == 0.2 & decision_pr == "no-go" ~ "FSR",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(cols = c(CorrectLabel, RiskLabel), names_to = "type", values_to = "label") %>%
  filter(!is.na(label)) %>%
  group_by(gate, borrowing, control.delta_SAM, control.ess_h, control.p, label) %>%
  summarize(prop = mean(proportion_pr), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = prop)


metrics_tab <- merge(metrics_tab, risk_df,
                     by = c("gate", "control.delta_SAM", "control.ess_h","control.p"),
                     all = TRUE)

metrics_tab_1 <- metrics_tab %>%
  mutate(mean_pmd = ifelse(borrowing == "Borrowing: No", NA, mean_pmd),
         sd_pmd = ifelse(borrowing == "Borrowing: No", NA, sd_pmd)) %>%
  pivot_longer(c(FGR, FSR, CGR, CSR, mean_pmd, sd_pmd)) %>%
  pivot_wider(names_from = c(borrowing, control.delta_SAM, control.ess_h),
              values_from = value, names_sort = T) %>%
  mutate(name = factor(name, levels = c("FGR", "FSR", "CGR", "CSR", "mean_pmd", "sd_pmd"))) %>%
  select(c(1,3,2,4,10:15)) %>%
  arrange(gate, name, control.p) %>%
  select(c(1:3, 5:10, 4))

# kbl(metrics_tab_1[, -1],
#     escape = F,
#     # format = "latex",
#     row.names = FALSE,
#     digits = 4,
#     col.names = c("", "control.p",
#                   rep(c("45", "90", "180"), 2),
#                   "Borrowing")) %>%
#   kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
#   kable_classic(full_width = FALSE) %>%
#   add_header_above(c(" " = 1, " " = 1,
#                      "Borrowing with delta = 0.1" = 3,
#                      "Borrowing with delta = 0.15" = 3,
#                      "No" = 1)) %>%
#   pack_rows("Gated Control: No", 1, 42) %>%
#   pack_rows("Gated Control: Yes", 43, 84) %>%
#   collapse_rows(columns = 1, valign = "middle")

metrics_tab_1 <- metrics_tab_1 %>%
  pivot_wider(names_from = gate, values_from = 4:10, names_vary = "slowest")

kbl(metrics_tab_1[, -1],
    escape = FALSE,
    row.names = FALSE,
    digits = 4,
    col.names = c(
      "control.p",
      rep(c("45", "90", "180"), 2), "",
      rep(c("45", "90", "180"), 2), ""
    )) %>%
  kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(
    " " = 1,
    "Borrowing with delta = 0.1" = 3,
    "Borrowing with delta = 0.15" = 3,
    "No Borrowing" = 1,
    "Borrowing with delta = 0.1" = 3,
    "Borrowing with delta = 0.15" = 3,
    "No Borrowing" = 1
  )) %>%
  add_header_above(c(
    " " = 1,
    "Gated Control: No" = 7,
    "Gated Control: Yes" = 7
  )) %>%
  pack_rows("FGR", 1, 7) %>%
  pack_rows("FSR", 8, 14) %>%
  pack_rows("CGR", 14, 21) %>%
  pack_rows("CSR", 22, 28) %>%
  pack_rows("Mean Control PMD", 29, 35) %>%
  pack_rows("SD Control PMD", 36, 42)


# + Visualize Risks -------------------------------------------------------

# Visualize risks
risk_df <- risk_df %>%
  pivot_longer(c(FSR, FGR)) %>%
  mutate(control.p = factor(control.p,
                            levels = c("0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45")),
         control.delta_SAM = factor(control.delta_SAM,
                                     levels = c("0.1", "0.15"),
                                     labels = c(expression(paste(delta, " = 0.1")),
                                                expression(paste(delta, " = 0.15")))),
         control.ess_h = factor(control.ess_h,
                                levels = c("45", "90", "180"),
                                labels = c(expression(n[ch*","*e] == 45),
                                           expression(n[ch*","*e] == 90),
                                           expression(n[ch*","*e] == 180))))

shade_df <- data.frame(
  control.delta_SAM = factor(c(0.1, 0.15), levels = c("0.1", "0.15"),
                              labels = c(expression(paste(delta, " = 0.1")),
                                         expression(paste(delta, " = 0.15")))),
  xmin = c(2, 1),
  xmax = c(6, 7),
  ymin = -Inf,
  ymax = Inf)

p_risk <- ggplot(risk_df %>% filter(borrowing == "Borrowing: Yes"),
                 aes(x = control.p, y = value, color = name, shape = gate)) +
  geom_point(size = 2) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  facet_grid(control.delta_SAM ~ control.ess_h,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(y = "Risk Value",
       color = "", shape = "") +
  xlab(~ paste(theta[c], " - ", theta[ch])) +
  geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
  )

p_risk

ggsave( "simulations/sim_OR_conflict_vs_risk_5_19.jpg", p_risk, width = 10, height = 5)



# + Visualize OC ----------------------------------------------------------

oc_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  select("i", "gate", "borrowing",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p.diff",
         "decision_pr", "proportion_pr", "proportion_ci") %>%
  mutate(decision_pr = factor(decision_pr, levels = c("no-go", "consider", "go"),
                              labels = c("No-Go", "Consider", "Go")),
         true_value.compare_true = ifelse(true_value.compare_true == 0.2, "0.2 (TV)", "0.1 (LRV)"))

oc_new <- data.frame()
oc_tmp <- na.omit(oc_df)
for (i in seq(1, nrow(oc_tmp), 3)) {
  tmp <- oc_tmp[i:(i+2),]
  tmp[tmp$decision_pr == "Go", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"]
  tmp[tmp$decision_pr == "Consider", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] + tmp[tmp$decision_pr == "Consider", "proportion_pr"]
  tmp[tmp$decision_pr == "No-Go", "label_ypos"] <- 1
  oc_new <- rbind(oc_new, tmp)
}

oc_new <- oc_new[, c("true_value.compare_true", "control.delta_SAM", "control.ess_h",
                     "gate", "borrowing",
                     "control.p.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]

oc_new$control.p.diff <- factor(oc_new$control.p.diff)
myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)

oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
oc_new$control.p.diff <- factor(oc_new$control.p.diff,
                                levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15"))
oc_new$control.delta_SAM <- factor(oc_new$control.delta_SAM,
                                    levels = c("0.1", "0.15"),
                                    labels = c(expression(paste(delta, " = 0.1")),
                                               expression(paste(delta, " = 0.15"))))
oc_new$control.ess_h <- factor(oc_new$control.ess_h,
                               levels = c("45", "90", "180"),
                               labels = c(expression(n[ch*","*e] == 45),
                                          expression(n[ch*","*e] == 90),
                                          expression(n[ch*","*e] == 180)))


p_list <- list()
for (i in levels(oc_new$true_value.compare_true)) {
  oc_tmp <- oc_new[oc_new$true_value.compare_true == i, ]

  p_list[[i]] <- ggplot(oc_tmp,
                        aes(x = control.p.diff, y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
              vjust = 1.6,
              fontface = "bold",
              size = 2) +
    scale_fill_manual(name = "Decision", values = myColors)+
    facet_grid(control.ess_h+gate~control.delta_SAM+borrowing,
               labeller = labeller(control.ess_h = label_parsed,
                                   control.delta_SAM = label_parsed,
                                   borrowing = label_parsed)) +
    ggtitle(bquote(Delta == .(i))) +
    labs(x =expression(paste(Delta[c], " = ", theta[c], " - ", theta[ch])), y = "Probability") +
    theme(legend.position = "bottom",
          title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 12)) +
    theme_bw()
}

library(patchwork)
p_oc <- p_list[[1]] + p_list[[2]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        title = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12))
p_oc
ggsave("simulations/sim_OR_zone_size_5_19.jpg", p_oc, width = 20, height = 15)




# Analysis II ----------------------------------------------------------------

bayes_results_all <- readRDS("simulations/sim_OR_bayes_results_df_5_19_betea01.rds")

# + Check posterior weights -----------------------------------------------

post_params_all <- bayes_results_all$post_params_all %>%
  mutate(across(
    .cols = -c(ends_with(".name")),
    .fns = as.numeric
  )) %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  filter(borrowing == "Borrowing: Yes") %>%
  select("nsim", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p.diff",
         "control.w_prior", "control.w_post") %>%
  mutate(control.p.diff = factor(control.p.diff,
                                 levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15")),
         control.delta_SAM = factor(control.delta_SAM,
                                    levels = c("0.1", "0.15"),
                                    labels = c(expression(paste(delta, " = 0.1")),
                                               expression(paste(delta, " = 0.15")))),
         control.ess_h = factor(control.ess_h,
                                levels = c("45", "90", "180"),
                                labels = c(expression(n[ch*","*e] == 45),
                                           expression(n[ch*","*e] == 90),
                                           expression(n[ch*","*e] == 180))))

shade_df <- data.frame(
  control.delta_SAM = factor(c(0.1, 0.15), levels = c("0.1", "0.15"),
                             labels = c(expression(paste(delta, " = 0.1")),
                                        expression(paste(delta, " = 0.15")))),
  xmin = c(2, 1),
  xmax = c(6, 7),
  ymin = -Inf,
  ymax = Inf)

p1 <- ggplot(post_params_all, aes(x = control.p.diff, y = control.w_prior)) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  facet_grid(control.delta_SAM ~ control.ess_h + gate,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(
    x = expression(paste(theta[c], " - ", theta[ch])),
    y = "Prior Weight",
    title = "Prior Weight of Informative Component of SAM Prior"
  ) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10))

p2 <- ggplot(post_params_all, aes(x = control.p.diff, y = control.w_post)) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  facet_grid(control.delta_SAM ~ control.ess_h + gate,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(x = expression(paste(theta[c], " - ", theta[ch])),
       y = "Posterior Weight",
       title = "Posterior Weight of Informative Component of SAM Prior") +
  ylim(c(0,1)) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10))

library(patchwork)
p1 / p2
ggsave(filename = "simulations/sim_OR_SAM_weight_5_19_beta01.jpg", p1 / p2,
       width = 13, height = 8)


# + Generate summary table for paper --------------------------------------

post_inference_all <-  bayes_results_all$post_inference_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  select("nsim", "borrowing", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p",
         "est2_lalonde") %>%
  pivot_wider(names_from = borrowing, values_from = est2_lalonde) %>%
  mutate(pmd = `Borrowing: Yes` - `Borrowing: No`) %>%
  select(-c(`Borrowing: Yes`, `Borrowing: No`, true_value.compare_true))

# Summary statistics of PMD on control
metrics_tab <- post_inference_all %>%
  group_by(gate, control.delta_SAM, control.ess_h, control.p) %>%
  summarize(mean_pmd = mean(pmd),
            sd_pmd = sd(pmd),
            .groups = "drop")

# risks
risk_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  filter(true_value.compare_true %in% c(0.1, 0.2),
         decision_pr %in% c("go", "no-go")) %>%
  mutate(
    CorrectLabel = case_when(
      true_value.compare_true == 0.2 & decision_pr == "go" ~ "CGR",
      true_value.compare_true == 0.1 & decision_pr == "no-go" ~ "CSR",
      TRUE ~ NA_character_
    ),
    RiskLabel = case_when(
      true_value.compare_true == 0.1 & decision_pr == "go" ~ "FGR",
      true_value.compare_true == 0.2 & decision_pr == "no-go" ~ "FSR",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(cols = c(CorrectLabel, RiskLabel), names_to = "type", values_to = "label") %>%
  filter(!is.na(label)) %>%
  group_by(gate, borrowing, control.delta_SAM, control.ess_h, control.p, label) %>%
  summarize(prop = mean(proportion_pr), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = prop)


metrics_tab <- merge(metrics_tab, risk_df,
                     by = c("gate", "control.delta_SAM", "control.ess_h","control.p"),
                     all = TRUE)

metrics_tab_1 <- metrics_tab %>%
  mutate(mean_pmd = ifelse(borrowing == "Borrowing: No", NA, mean_pmd),
         sd_pmd = ifelse(borrowing == "Borrowing: No", NA, sd_pmd)) %>%
  pivot_longer(c(FGR, FSR, CGR, CSR, mean_pmd, sd_pmd)) %>%
  pivot_wider(names_from = c(borrowing, control.delta_SAM, control.ess_h),
              values_from = value, names_sort = T) %>%
  mutate(name = factor(name, levels = c("FGR", "FSR", "CGR", "CSR", "mean_pmd", "sd_pmd"))) %>%
  select(c(1,3,2,4,10:15)) %>%
  arrange(gate, name, control.p) %>%
  select(c(1:3, 5:10, 4))

# kbl(metrics_tab_1[, -1],
#     escape = F,
#     # format = "latex",
#     row.names = FALSE,
#     digits = 4,
#     col.names = c("", "control.p",
#                   rep(c("45", "90", "180"), 2),
#                   "Borrowing")) %>%
#   kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
#   kable_classic(full_width = FALSE) %>%
#   add_header_above(c(" " = 1, " " = 1,
#                      "Borrowing with delta = 0.1" = 3,
#                      "Borrowing with delta = 0.15" = 3,
#                      "No" = 1)) %>%
#   pack_rows("Gated Control: No", 1, 42) %>%
#   pack_rows("Gated Control: Yes", 43, 84) %>%
#   collapse_rows(columns = 1, valign = "middle")

metrics_tab_1 <- metrics_tab_1 %>%
  pivot_wider(names_from = gate, values_from = 4:10, names_vary = "slowest")

kbl(metrics_tab_1[, -1],
    escape = FALSE,
    row.names = FALSE,
    digits = 4,
    col.names = c(
      "control.p",
      rep(c("45", "90", "180"), 2), "",
      rep(c("45", "90", "180"), 2), ""
    )) %>%
  kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(
    " " = 1,
    "Borrowing with delta = 0.1" = 3,
    "Borrowing with delta = 0.15" = 3,
    "No Borrowing" = 1,
    "Borrowing with delta = 0.1" = 3,
    "Borrowing with delta = 0.15" = 3,
    "No Borrowing" = 1
  )) %>%
  add_header_above(c(
    " " = 1,
    "Gated Control: No" = 7,
    "Gated Control: Yes" = 7
  )) %>%
  pack_rows("FGR", 1, 7) %>%
  pack_rows("FSR", 8, 14) %>%
  pack_rows("CGR", 14, 21) %>%
  pack_rows("CSR", 22, 28) %>%
  pack_rows("Mean Control PMD", 29, 35) %>%
  pack_rows("SD Control PMD", 36, 42)


# + Visualize Risks -------------------------------------------------------

# Visualize risks
risk_df <- risk_df %>%
  pivot_longer(c(FSR, FGR)) %>%
  mutate(control.p = factor(control.p,
                            levels = c("0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45")),
         control.delta_SAM = factor(control.delta_SAM,
                                    levels = c("0.1", "0.15"),
                                    labels = c(expression(paste(delta, " = 0.1")),
                                               expression(paste(delta, " = 0.15")))),
         control.ess_h = factor(control.ess_h,
                                levels = c("45", "90", "180"),
                                labels = c(expression(n[ch*","*e] == 45),
                                           expression(n[ch*","*e] == 90),
                                           expression(n[ch*","*e] == 180))))

shade_df <- data.frame(
  control.delta_SAM = factor(c(0.1, 0.15), levels = c("0.1", "0.15"),
                             labels = c(expression(paste(delta, " = 0.1")),
                                        expression(paste(delta, " = 0.15")))),
  xmin = c(2, 1),
  xmax = c(6, 7),
  ymin = -Inf,
  ymax = Inf)

p_risk <- ggplot(risk_df %>% filter(borrowing == "Borrowing: Yes"),
                 aes(x = control.p, y = value, color = name, shape = gate)) +
  geom_point(size = 2) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  facet_grid(control.delta_SAM ~ control.ess_h,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(y = "Risk Value",
       color = "", shape = "") +
  xlab(~ paste(theta[c], " - ", theta[ch])) +
  geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
  )

p_risk

ggsave( "simulations/sim_OR_conflict_vs_risk_5_19_beta01.jpg", p_risk, width = 10, height = 5)



# + Visualize OC ----------------------------------------------------------

oc_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  select("i", "gate", "borrowing",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p.diff",
         "decision_pr", "proportion_pr", "proportion_ci") %>%
  mutate(decision_pr = factor(decision_pr, levels = c("no-go", "consider", "go"),
                              labels = c("No-Go", "Consider", "Go")),
         true_value.compare_true = ifelse(true_value.compare_true == 0.2, "0.2 (TV)", "0.1 (LRV)"))

oc_new <- data.frame()
oc_tmp <- na.omit(oc_df)
for (i in seq(1, nrow(oc_tmp), 3)) {
  tmp <- oc_tmp[i:(i+2),]
  tmp[tmp$decision_pr == "Go", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"]
  tmp[tmp$decision_pr == "Consider", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] + tmp[tmp$decision_pr == "Consider", "proportion_pr"]
  tmp[tmp$decision_pr == "No-Go", "label_ypos"] <- 1
  oc_new <- rbind(oc_new, tmp)
}

oc_new <- oc_new[, c("true_value.compare_true", "control.delta_SAM", "control.ess_h",
                     "gate", "borrowing",
                     "control.p.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]

oc_new$control.p.diff <- factor(oc_new$control.p.diff)
myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)

oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
oc_new$control.p.diff <- factor(oc_new$control.p.diff,
                                levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15"))
oc_new$control.delta_SAM <- factor(oc_new$control.delta_SAM,
                                   levels = c("0.1", "0.15"),
                                   labels = c(expression(paste(delta, " = 0.1")),
                                              expression(paste(delta, " = 0.15"))))
oc_new$control.ess_h <- factor(oc_new$control.ess_h,
                               levels = c("45", "90", "180"),
                               labels = c(expression(n[ch*","*e] == 45),
                                          expression(n[ch*","*e] == 90),
                                          expression(n[ch*","*e] == 180)))


p_list <- list()
for (i in levels(oc_new$true_value.compare_true)) {
  oc_tmp <- oc_new[oc_new$true_value.compare_true == i, ]

  p_list[[i]] <- ggplot(oc_tmp,
                        aes(x = control.p.diff, y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
              vjust = 1.6,
              fontface = "bold",
              size = 2) +
    scale_fill_manual(name = "Decision", values = myColors)+
    facet_grid(control.ess_h+gate~control.delta_SAM+borrowing,
               labeller = labeller(control.ess_h = label_parsed,
                                   control.delta_SAM = label_parsed,
                                   borrowing = label_parsed)) +
    ggtitle(bquote(Delta == .(i))) +
    labs(x =expression(paste(Delta[c], " = ", theta[c], " - ", theta[ch])), y = "Probability") +
    theme(legend.position = "bottom",
          title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 12)) +
    theme_bw()
}

library(patchwork)
p_oc <- p_list[[1]] + p_list[[2]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        title = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12))
p_oc
ggsave("simulations/sim_OR_zone_size_5_19_beta01.jpg", p_oc, width = 20, height = 15)





# Analysis III ----------------------------------------------------------------

bayes_results_all <- readRDS("simulations/sim_OR_bayes_results_df_5_19_betea001.rds")

# + Check posterior weights -----------------------------------------------

post_params_all <- bayes_results_all$post_params_all %>%
  mutate(across(
    .cols = -c(ends_with(".name")),
    .fns = as.numeric
  )) %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  filter(borrowing == "Borrowing: Yes") %>%
  select("nsim", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p.diff",
         "control.w_prior", "control.w_post") %>%
  mutate(control.p.diff = factor(control.p.diff,
                                 levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15")),
         control.delta_SAM = factor(control.delta_SAM,
                                    levels = c("0.1", "0.15"),
                                    labels = c(expression(paste(delta, " = 0.1")),
                                               expression(paste(delta, " = 0.15")))),
         control.ess_h = factor(control.ess_h,
                                levels = c("45", "90", "180"),
                                labels = c(expression(n[ch*","*e] == 45),
                                           expression(n[ch*","*e] == 90),
                                           expression(n[ch*","*e] == 180))))

shade_df <- data.frame(
  control.delta_SAM = factor(c(0.1, 0.15), levels = c("0.1", "0.15"),
                             labels = c(expression(paste(delta, " = 0.1")),
                                        expression(paste(delta, " = 0.15")))),
  xmin = c(2, 1),
  xmax = c(6, 7),
  ymin = -Inf,
  ymax = Inf)

p1 <- ggplot(post_params_all, aes(x = control.p.diff, y = control.w_prior)) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  facet_grid(control.delta_SAM ~ control.ess_h + gate,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(
    x = expression(paste(theta[c], " - ", theta[ch])),
    y = "Prior Weight",
    title = "Prior Weight of Informative Component of SAM Prior"
  ) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10))

p2 <- ggplot(post_params_all, aes(x = control.p.diff, y = control.w_post)) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  facet_grid(control.delta_SAM ~ control.ess_h + gate,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(x = expression(paste(theta[c], " - ", theta[ch])),
       y = "Posterior Weight",
       title = "Posterior Weight of Informative Component of SAM Prior") +
  ylim(c(0,1)) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10))

library(patchwork)
p1 / p2
ggsave(filename = "simulations/sim_OR_SAM_weight_5_19_beta001.jpg", p1 / p2,
       width = 13, height = 8)


# + Generate summary table for paper --------------------------------------

post_inference_all <-  bayes_results_all$post_inference_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  select("nsim", "borrowing", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p",
         "est2_lalonde") %>%
  pivot_wider(names_from = borrowing, values_from = est2_lalonde) %>%
  mutate(pmd = `Borrowing: Yes` - `Borrowing: No`) %>%
  select(-c(`Borrowing: Yes`, `Borrowing: No`, true_value.compare_true))

# Summary statistics of PMD on control
metrics_tab <- post_inference_all %>%
  group_by(gate, control.delta_SAM, control.ess_h, control.p) %>%
  summarize(mean_pmd = mean(pmd),
            sd_pmd = sd(pmd),
            .groups = "drop")

# risks
risk_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  filter(true_value.compare_true %in% c(0.1, 0.2),
         decision_pr %in% c("go", "no-go")) %>%
  mutate(
    CorrectLabel = case_when(
      true_value.compare_true == 0.2 & decision_pr == "go" ~ "CGR",
      true_value.compare_true == 0.1 & decision_pr == "no-go" ~ "CSR",
      TRUE ~ NA_character_
    ),
    RiskLabel = case_when(
      true_value.compare_true == 0.1 & decision_pr == "go" ~ "FGR",
      true_value.compare_true == 0.2 & decision_pr == "no-go" ~ "FSR",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(cols = c(CorrectLabel, RiskLabel), names_to = "type", values_to = "label") %>%
  filter(!is.na(label)) %>%
  group_by(gate, borrowing, control.delta_SAM, control.ess_h, control.p, label) %>%
  summarize(prop = mean(proportion_pr), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = prop)


metrics_tab <- merge(metrics_tab, risk_df,
                     by = c("gate", "control.delta_SAM", "control.ess_h","control.p"),
                     all = TRUE)

metrics_tab_1 <- metrics_tab %>%
  mutate(mean_pmd = ifelse(borrowing == "Borrowing: No", NA, mean_pmd),
         sd_pmd = ifelse(borrowing == "Borrowing: No", NA, sd_pmd)) %>%
  pivot_longer(c(FGR, FSR, CGR, CSR, mean_pmd, sd_pmd)) %>%
  pivot_wider(names_from = c(borrowing, control.delta_SAM, control.ess_h),
              values_from = value, names_sort = T) %>%
  mutate(name = factor(name, levels = c("FGR", "FSR", "CGR", "CSR", "mean_pmd", "sd_pmd"))) %>%
  select(c(1,3,2,4,10:15)) %>%
  arrange(gate, name, control.p) %>%
  select(c(1:3, 5:10, 4))

# kbl(metrics_tab_1[, -1],
#     escape = F,
#     # format = "latex",
#     row.names = FALSE,
#     digits = 4,
#     col.names = c("", "control.p",
#                   rep(c("45", "90", "180"), 2),
#                   "Borrowing")) %>%
#   kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
#   kable_classic(full_width = FALSE) %>%
#   add_header_above(c(" " = 1, " " = 1,
#                      "Borrowing with delta = 0.1" = 3,
#                      "Borrowing with delta = 0.15" = 3,
#                      "No" = 1)) %>%
#   pack_rows("Gated Control: No", 1, 42) %>%
#   pack_rows("Gated Control: Yes", 43, 84) %>%
#   collapse_rows(columns = 1, valign = "middle")

metrics_tab_1 <- metrics_tab_1 %>%
  pivot_wider(names_from = gate, values_from = 4:10, names_vary = "slowest")

kbl(metrics_tab_1[, -1],
    escape = FALSE,
    row.names = FALSE,
    digits = 4,
    col.names = c(
      "control.p",
      rep(c("45", "90", "180"), 2), "",
      rep(c("45", "90", "180"), 2), ""
    )) %>%
  kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(
    " " = 1,
    "Borrowing with delta = 0.1" = 3,
    "Borrowing with delta = 0.15" = 3,
    "No Borrowing" = 1,
    "Borrowing with delta = 0.1" = 3,
    "Borrowing with delta = 0.15" = 3,
    "No Borrowing" = 1
  )) %>%
  add_header_above(c(
    " " = 1,
    "Gated Control: No" = 7,
    "Gated Control: Yes" = 7
  )) %>%
  pack_rows("FGR", 1, 7) %>%
  pack_rows("FSR", 8, 14) %>%
  pack_rows("CGR", 14, 21) %>%
  pack_rows("CSR", 22, 28) %>%
  pack_rows("Mean Control PMD", 29, 35) %>%
  pack_rows("SD Control PMD", 36, 42)


# + Visualize Risks -------------------------------------------------------

# Visualize risks
risk_df <- risk_df %>%
  pivot_longer(c(FSR, FGR)) %>%
  mutate(control.p = factor(control.p,
                            levels = c("0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45")),
         control.delta_SAM = factor(control.delta_SAM,
                                    levels = c("0.1", "0.15"),
                                    labels = c(expression(paste(delta, " = 0.1")),
                                               expression(paste(delta, " = 0.15")))),
         control.ess_h = factor(control.ess_h,
                                levels = c("45", "90", "180"),
                                labels = c(expression(n[ch*","*e] == 45),
                                           expression(n[ch*","*e] == 90),
                                           expression(n[ch*","*e] == 180))))

shade_df <- data.frame(
  control.delta_SAM = factor(c(0.1, 0.15), levels = c("0.1", "0.15"),
                             labels = c(expression(paste(delta, " = 0.1")),
                                        expression(paste(delta, " = 0.15")))),
  xmin = c(2, 1),
  xmax = c(6, 7),
  ymin = -Inf,
  ymax = Inf)

p_risk <- ggplot(risk_df %>% filter(borrowing == "Borrowing: Yes"),
                 aes(x = control.p, y = value, color = name, shape = gate)) +
  geom_point(size = 2) +
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        group = control.delta_SAM),
    inherit.aes = FALSE,
    fill = "blue", alpha = 0.1
  ) +
  facet_grid(control.delta_SAM ~ control.ess_h,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta_SAM = label_parsed)) +
  labs(y = "Risk Value",
       color = "", shape = "") +
  xlab(~ paste(theta[c], " - ", theta[ch])) +
  geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
  )

p_risk

ggsave( "simulations/sim_OR_conflict_vs_risk_5_19_beta001.jpg", p_risk, width = 10, height = 5)



# + Visualize OC ----------------------------------------------------------

oc_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  select("i", "gate", "borrowing",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p.diff",
         "decision_pr", "proportion_pr", "proportion_ci") %>%
  mutate(decision_pr = factor(decision_pr, levels = c("no-go", "consider", "go"),
                              labels = c("No-Go", "Consider", "Go")),
         true_value.compare_true = ifelse(true_value.compare_true == 0.2, "0.2 (TV)", "0.1 (LRV)"))

oc_new <- data.frame()
oc_tmp <- na.omit(oc_df)
for (i in seq(1, nrow(oc_tmp), 3)) {
  tmp <- oc_tmp[i:(i+2),]
  tmp[tmp$decision_pr == "Go", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"]
  tmp[tmp$decision_pr == "Consider", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] + tmp[tmp$decision_pr == "Consider", "proportion_pr"]
  tmp[tmp$decision_pr == "No-Go", "label_ypos"] <- 1
  oc_new <- rbind(oc_new, tmp)
}

oc_new <- oc_new[, c("true_value.compare_true", "control.delta_SAM", "control.ess_h",
                     "gate", "borrowing",
                     "control.p.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]

oc_new$control.p.diff <- factor(oc_new$control.p.diff)
myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)

oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
oc_new$control.p.diff <- factor(oc_new$control.p.diff,
                                levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15"))
oc_new$control.delta_SAM <- factor(oc_new$control.delta_SAM,
                                   levels = c("0.1", "0.15"),
                                   labels = c(expression(paste(delta, " = 0.1")),
                                              expression(paste(delta, " = 0.15"))))
oc_new$control.ess_h <- factor(oc_new$control.ess_h,
                               levels = c("45", "90", "180"),
                               labels = c(expression(n[ch*","*e] == 45),
                                          expression(n[ch*","*e] == 90),
                                          expression(n[ch*","*e] == 180)))


p_list <- list()
for (i in levels(oc_new$true_value.compare_true)) {
  oc_tmp <- oc_new[oc_new$true_value.compare_true == i, ]

  p_list[[i]] <- ggplot(oc_tmp,
                        aes(x = control.p.diff, y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
              vjust = 1.6,
              fontface = "bold",
              size = 2) +
    scale_fill_manual(name = "Decision", values = myColors)+
    facet_grid(control.ess_h+gate~control.delta_SAM+borrowing,
               labeller = labeller(control.ess_h = label_parsed,
                                   control.delta_SAM = label_parsed,
                                   borrowing = label_parsed)) +
    ggtitle(bquote(Delta == .(i))) +
    labs(x =expression(paste(Delta[c], " = ", theta[c], " - ", theta[ch])), y = "Probability") +
    theme(legend.position = "bottom",
          title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 12)) +
    theme_bw()
}

library(patchwork)
p_oc <- p_list[[1]] + p_list[[2]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        title = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12))
p_oc
ggsave("simulations/sim_OR_zone_size_5_19_beta001.jpg", p_oc, width = 20, height = 15)



#####

oc_all <- bayes_results_all$oc_all
oc_all <- oc_all[, c("i", "true_value.compare_true", "control.n", "borrowing",
                     "treatment.p", "control.p", "control_h.p",
                     "decision_pr", "proportion_pr", "proportion_ci")]

# oc_all$decision_pr <- factor(oc_all$decision_pr, levels = c("no-go", "consider", "go", "NA"))
oc_all$decision_pr <- factor(oc_all$decision_pr, levels = c("no-go", "consider", "go"))
oc_all$true_value.compare_true = ifelse(oc_all$true_value.compare_true == 0.1, "0.1 (LRV)",
                                        ifelse(oc_all$true_value.compare_true == 0.2, "0.2 (TV)",
                                               oc_all$true_value.compare_true))

oc_all$control.p <- as.numeric(oc_all$control.p)
oc_all$control_h.p <- as.numeric(oc_all$control_h.p)
oc_all$control.p.diff <- oc_all$control.p - oc_all$control_h.p
oc_all$borrowing <- factor(oc_all$borrowing,
                           levels = c("No", "Yes"),
                           labels=c("Borrowing: No", "Borrowing: Yes"))

metrics <- oc_all[oc_all$true_value.compare_true %in% c("0.2 (TV)", "0.1 (LRV)"),]
metrics <- oc_all[(oc_all$true_value.compare_true == "0.1 (LRV)" & oc_all$decision_pr == "go") | (oc_all$true_value.compare_true == "0.2 (TV)" & oc_all$decision_pr == "no-go"),]
false_go_risk <- metrics[metrics$true_value.compare_true == "0.1 (LRV)", c("control.n", "borrowing", "proportion_pr", "control.p.diff")]
colnames(false_go_risk)[3] <- "Type I Error (FGR)"
false_stop_risk <- metrics[metrics$true_value.compare_true == "0.2 (TV)", c("control.n", "borrowing", "proportion_pr", "control.p.diff")]
colnames(false_stop_risk)[3] <- "Type II Error (FSR)"
metrics_df <- merge(false_go_risk, false_stop_risk, by = c("control.n", "borrowing", "control.p.diff"))

metrics_long <- reshape2::melt(metrics_df, id.vars = c("control.n", "borrowing", "control.p.diff"),
                               measure.vars = c("Type I Error (FGR)", "Type II Error (FSR)"),
                               variable.name = "risk_type", value.name = "risk_value")


oc_new <- data.frame()
oc_tmp <- na.omit(oc_all)
for (i in seq(1, nrow(oc_tmp), 3)) {
  tmp <- oc_tmp[i:(i+2),]
  tmp[tmp$decision_pr == "go", "label_ypos"] <- tmp[tmp$decision_pr == "go", "proportion_pr"]
  tmp[tmp$decision_pr == "consider", "label_ypos"] <- tmp[tmp$decision_pr == "go", "proportion_pr"] + tmp[tmp$decision_pr == "consider", "proportion_pr"]
  tmp[tmp$decision_pr == "no-go", "label_ypos"] <- 1
  oc_new <- rbind(oc_new, tmp)
}

oc_new <- oc_new[, c("true_value.compare_true", "control.n", "borrowing",
                     "control.p", "control.p.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]
oc_new$control.p.diff <- factor(oc_new$control.p.diff)
myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)
oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)

p_list <- list()
for (i in levels(oc_new$true_value.compare_true)) {
  oc_tmp <- oc_new[oc_new$true_value.compare_true == i, ]
  oc_tmp$control.p.diff <- factor(oc_tmp$control.p.diff,
                                  levels = c("-0.2", "-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15", "0.2"),
                                  labels=c(expression(paste(theta[c], " - ", theta[ch], " = -0.2")),
                                           expression(paste(theta[c], " - ", theta[ch], " = -0.15")),
                                           expression(paste(theta[c], " - ", theta[ch], " = -0.1")),
                                           expression(paste(theta[c], " - ", theta[ch], " = -0.05")),
                                           expression(paste(theta[c], " - ", theta[ch], " = 0")),
                                           expression(paste(theta[c], " - ", theta[ch], " = 0.05")),
                                           expression(paste(theta[c], " - ", theta[ch], " = 0.1")),
                                           expression(paste(theta[c], " - ", theta[ch], " = 0.15")),
                                           expression(paste(theta[c], " - ", theta[ch], " = 0.2"))))
  p_list[[i]] <- ggplot(oc_tmp,
                        aes(x = as.factor(control.n), y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.1)),
              vjust = 1.6,
              fontface = "bold",
              size = 2) +
    scale_fill_manual(name = "Decision", values = myColors)+
    facet_grid(control.p.diff~borrowing,
               labeller = labeller(control.p.diff = label_parsed)) +
    ggtitle(bquote(Delta == .(i))) +
    labs(x = "Sample size per arm", y = "Probability") +
    theme(legend.position = "bottom",
          title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 12)) +
    theme_bw()
}
library(patchwork)
p_oc <- p_list[[1]] + p_list[[2]] + p_list[[3]] +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom",
        title = element_text(size = 15),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,0),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13))
ggsave("simulations/sim_OR_zone_size_5_19.jpg",p_oc, width = 14, height = 14)


oc2 <- oc_all %>%
  select(control.p.diff, control.n, true_value.compare_true, borrowing, decision_pr, proportion_pr) %>%
  mutate(proportion_pr = scales::percent(proportion_pr, 0.01)) %>%
  arrange(decision_pr) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = c(decision_pr, borrowing), values_from = proportion_pr) %>%
  arrange(control.n) %>%
  arrange(control.p.diff)


kbl(oc2[,-1], escape = T, row.names = F, digits = 4,
    format    = "latex",
    booktabs  = T,
    col.names = c("n", "Delta", rep(c("Pr(No-Go)", "Pr(Consider)","Pr(Go)"), 2))) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 2, "Borrowing: No" = 3, "Borrowing: Yes" = 3)) %>%
  pack_rows(index = c("Case 1" = 15,
                      "Case 2" = 15,
                      "Case 3" = 15,
                      "Case 4" = 15,
                      "Case 5" = 15,
                      "Case 6" = 15,
                      "Case 7" = 15,
                      "Case 8" = 15,
                      "Case 9" = 15))







####
metrics_post_dist_all <- bayes_results_all$metrics_post_dist_all
metrics_post_dist_all <- metrics_post_dist_all  %>%
  mutate(across(
    .cols = -c(ends_with(".name")),
    .fns = as.numeric
  ))
metrics_post_dist_all <- metrics_post_dist_all[, c("true_value.compare_true",
                                                   "control.n", "control.delta_SAM",
                                                   "control.ess_h","control.w",
                                                   "treatment.p", "control.p",
                                                   "control_h.p", "bias_avg_rate_diff",
                                                   "sd_avg", "sd_empirical", "cp")]

metrics_post_dist_all$borrowing <- ifelse(is.na(metrics_post_dist_all$control.w),
                                          "Yes", "No")

metrics_post_dist_all$control.p.diff <- metrics_post_dist_all$control.p - metrics_post_dist_all$control_h.p
metrics_post_dist_all$true_value.compare_true = ifelse(metrics_post_dist_all$true_value.compare_true == 0.1, "0.1 (LRV)",
                                                       ifelse(metrics_post_dist_all$true_value.compare_true == 0.2, "0.2 (TV)",
                                                              metrics_post_dist_all$true_value.compare_true))

metrics_df <- metrics_post_dist_all %>%
  select(control.p.diff, control.ess_h, control.delta_SAM,
         true_value.compare_true, borrowing,
         bias_avg_rate_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true)

metrics_df$control.p.diff <- factor(metrics_df$control.p.diff,
                                       levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15"),
                                       labels=c(expression(paste(theta[c], " - ", theta[ch], " = -0.15")),
                                                expression(paste(theta[c], " - ", theta[ch], " = -0.1")),
                                                expression(paste(theta[c], " - ", theta[ch], " = -0.05")),
                                                expression(paste(theta[c], " - ", theta[ch], " = 0")),
                                                expression(paste(theta[c], " - ", theta[ch], " = 0.05")),
                                                expression(paste(theta[c], " - ", theta[ch], " = 0.1")),
                                                expression(paste(theta[c], " - ", theta[ch], " = 0.15"))))

metrics_df$true_value.compare_true <- factor(metrics_df$true_value.compare_true,
                                       levels = c("0.1 (LRV)", "0.2 (TV)"),
                                       labels=c(expression(paste(Delta, " = 0.1")),
                                                expression(paste(Delta, " = 0.2"))))




p1 <- ggplot(metrics_df,
             aes(x = control.ess_h, y = bias_avg_rate_diff, color = borrowing)) +
  geom_point(position = position_dodge(width = 0.1), size = 2) +
  facet_grid(true_value.compare_true~control.p.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.p.diff = label_parsed)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "Average Bias", color= "") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

p3 <- ggplot(metrics_df,
             aes(x = control.n, color = borrowing, y = cp)) +
  geom_point(position = position_dodge(width = 0.1), size = 2) +
  facet_grid(true_value.compare_true~control.p.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.p.diff = label_parsed)) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "Coverage Probability (%)", color = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

metrics_long <- metrics_df %>%
  pivot_longer(c(sd_avg, sd_empirical), names_to = "sd_type", values_to = "sd_value")

p2 <- ggplot(metrics_long,
             aes(x = control.n, y = sd_value, shape = sd_type, color = borrowing)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  ggh4x::facet_grid2(true_value.compare_true~control.p.diff,
                     labeller = labeller(true_value.compare_true = label_parsed,
                                         control.p.diff = label_parsed)) +
  labs(x = "Sample size per arm", y = "Standard Deviation",
       shape = "",
       color = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

library(patchwork)
p <- p1 / p2 / p3 +
  plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')') &
  theme(strip.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,0),
        plot.tag = element_text(size = 14)
        )
p
ggsave(filename = "simulations/sim_OR_dist_check_metrics_5_19_1.jpg", p, width = 15, height = 12)

p <- p1 / p3 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')') &
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 15),
        legend.position = "bottom",
        plot.tag = element_text(size = 17)
  )
ggsave(filename = "simulations/sim_OR_dist_check_metrics_5_19_2.jpg", p, width = 15.5, height = 10)


metrics_df <- metrics_post_dist_all %>%
  select(control.p.diff, control.n,true_value.compare_true, borrowing,
         bias_avg_rate_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = borrowing, values_from = c(bias_avg_rate_diff, sd_avg, sd_empirical, cp)) %>%
  arrange(control.n) %>%
  arrange(control.p.diff)
metrics_df <- metrics_df[, c("control.p.diff", "control.n", "true_value.compare_true",
                             "bias_avg_rate_diff_Borrowing: No", "sd_avg_Borrowing: No",
                             "sd_empirical_Borrowing: No", "cp_Borrowing: No",
                             "bias_avg_rate_diff_Borrowing: Yes", "sd_avg_Borrowing: Yes",
                             "sd_empirical_Borrowing: Yes", "cp_Borrowing: Yes")]

kbl(metrics_df[,-1], escape = F, row.names = F, digits = 4,
    format    = "latex",
    booktabs  = T,
    col.names = c("Sample size", "Delta",
                  "BIAS", "SE-avg","SD-emp", "CP",
                  "BIAS", "SE-avg","SD-emp", "CP")) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 2, "Borrowing: No" = 4, "Borrowing: Yes" = 4)) %>%
  pack_rows(index = c("Case 1" = 15,
                      "Case 2" = 15,
                      "Case 3" = 15,
                      "Case 4" = 15,
                      "Case 5" = 15,
                      "Case 6" = 15,
                      "Case 7" = 15,
                      "Case 8" = 15,
                      "Case 9" = 15))


