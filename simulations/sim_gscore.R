library(kableExtra)
library(LalondeBayesBorrow)
library(dplyr)

# Simulation I ------------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 200, # Sample size for the historical control group.
  ctrl_h_prob = seq(0, 0.4, 0.1), # Probability of zero g-scores in the historical control group.
  ctrl_h_mu = 0.003, # Mean for the non-zero g-scores in the historical control group.
  ctrl_h_sigma = 1, # Standard deviation for the non-zero g-scores in the historical control group.
  stringsAsFactors = FALSE
)

param_grid <- expand.grid(
  trt_n = c(20, 40, 60, 80, 100, 150, 200, 250, 300), # Sample size for the treatment group.
  ctrl_n = c(20, 40, 60, 80, 100, 150, 200, 250, 300), # Sample size for the control group.
  ctrl_prob = seq(0, 0.4, 0.1), # Probability of zero g-scores in the control group.
  trt_mu = 0.0018, # Mean for the non-zero g-scores in the treatment group.
  ctrl_mu = 0.003, # Mean for the non-zero g-scores in the control group.
  trt_sigma = 1, # Standard deviation for the non-zero g-scores in the treatment group.
  ctrl_sigma = 1, # Standard deviation for the non-zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(trt_prob = ctrl_prob + 0.05) %>%
  filter(trt_n == ctrl_n)


param_grid_h <- expand.grid(
  ctrl_h_n = 200, # Sample size for the historical control group.
  ctrl_h_prob = seq(0, 0.4, 0.1), # Probability of zero g-scores in the historical control group.
  ctrl_h_mu = 0.01, # Mean for the non-zero g-scores in the historical control group.
  ctrl_h_sigma = 1, # Standard deviation for the non-zero g-scores in the historical control group.
  stringsAsFactors = FALSE
)

param_grid <- expand.grid(
  trt_n = c(20, 40, 60, 80, 100, 150, 200, 250, 300), # Sample size for the treatment group.
  ctrl_n = c(20, 40, 60, 80, 100, 150, 200, 250, 300), # Sample size for the control group.
  ctrl_prob = seq(0, 0.4, 0.1), # Probability of zero g-scores in the control group.
  trt_mu = 0.005, # Mean for the non-zero g-scores in the treatment group.
  ctrl_mu = 0.01, # Mean for the non-zero g-scores in the control group.
  trt_sigma = 1, # Standard deviation for the non-zero g-scores in the treatment group.
  ctrl_sigma = 1, # Standard deviation for the non-zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(trt_prob = ctrl_prob + 0.05) %>%
  filter(trt_n == ctrl_n)


data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "g-score")
data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "g-score")

prior_params_list <- list(
  Yes = list(control.delta = log(0.85), control.gate = log(0.85), treatment.w = 0, control.w = NULL), # borrowing on the control arm
  No = list(control.delta = log(0.85), treatment.w = 0, control.w = 0) # no borrowing
)

## Run simulations
nsim = 100000
all_historical_est <- lapply(seq_along(data_gen_params_list_h), function(i) {
  data_gen_params <- data_gen_params_list_h[[i]]
  cat("Running", nsim, "simulations for data_gen_params set", i, "\n\n")
  for (arm in names(data_gen_params)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params[[arm]]$n, "\n")
    cat("\t  prob:", data_gen_params[[arm]]$prob, "\n")
    cat("\t  mu:", data_gen_params[[arm]]$mu, "\n")
    cat("\t  sigma:", data_gen_params[[arm]]$sigma, "\n\n")
  }
  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- sapply(data_gen_params, function(x) x$n)
  prob_list <- sapply(data_gen_params, function(x) x$prob)
  mu_list <- sapply(data_gen_params, function(x) x$mu)
  sigma_list <- sapply(data_gen_params, function(x) x$sigma)

  # Generate simulation dataset
  data <- data_gen_gscore(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list,
                          mu_list = mu_list, sigma_list = sigma_list)
  median_est <- data$median_est
  get_true_gscore_median(p=prob_list, mu = mu_list, 1)
  log(get_true_gscore_median(p=prob_list, mu = mu_list, 1))

  return(data.frame(arm = "control_h",
                    n = 200,
                    median_adjusted1 = mean(log(median_est$median_adjusted)),
                    median_adjusted2 = log(mean(median_est$median_adjusted)),
                    median_sd1 = sd(log(median_est$median_adjusted)),
                    median_sd2 = sd(median_est$median_adjusted)/mean(median_est$median_adjusted)))
})


nsim = 10000
metrics_approx_dist <- data.frame()
bayes_results <- list()
for(j in seq_along(all_historical_est)){
  historical_est <- all_historical_est[[j]]
  data_gen_params_h <- data_gen_params_list_h[[j]]

  for (i in seq_along(data_gen_params_list)) {
    data_gen_params <- data_gen_params_list[[i]]
    if(data_gen_params$control$prob == data_gen_params_h$control_h$prob){
      cat("Running", nsim, "simulations for data_gen_params set", i, "\n\n")
      for (arm in names(data_gen_params)) {
        cat("\t Arm:", arm, "\n")
        cat("\t  n:", data_gen_params[[arm]]$n, "\n")
        cat("\t  prob:", data_gen_params[[arm]]$prob, "\n")
        cat("\t  mu:", data_gen_params[[arm]]$mu, "\n")
        cat("\t  sigma:", data_gen_params[[arm]]$sigma, "\n\n")
      }
      n_arms <- length(data_gen_params)
      arm_names <- sapply(data_gen_params, function(x) x$name)
      n_list <- sapply(data_gen_params, function(x) x$n)
      prob_list <- sapply(data_gen_params, function(x) x$prob)
      mu_list <- sapply(data_gen_params, function(x) x$mu)
      sigma_list <- sapply(data_gen_params, function(x) x$sigma)

      # Generate simulation dataset
      data <- data_gen_gscore(n_arms = n_arms,
                              arm_names = arm_names,
                              nsim = nsim, n_list = n_list,
                              prob_list = prob_list,
                              mu_list = mu_list, sigma_list = sigma_list)

      settings <- c(list(true_value = data$true_value), data_gen_params, data_gen_params_h)
      settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
        mutate(across(
          .cols = -c(ends_with(".name")),
          .fns = as.numeric
        ))

      # Assess normal approximation of g-score
      cat("Start assessing the normal approximation of log g-score median estimate.\n")
      metrics_approx_dist <- rbind(metrics_approx_dist,
                                   cbind(settings, eval_gscore_approx_dist(data)))

      median_est <- data$median_est
      median_est$log_median_sd <- median_est$median_sd / median_est$median_adjusted
      median_est$log_median_adjusted <- log(median_est$median_adjusted)

      median_est_h <- data.frame(nsim = 1:nsim,
                                 arm = "control_h",
                                 n = 200,
                                 log_median_adjusted = historical_est$median_adjusted2,
                                 log_median_sd = historical_est$median_sd2)

      median_est <- bind_rows(median_est, median_est_h) %>%
        mutate(log_median_sd = ifelse(log_median_sd == 0, 0.000001, log_median_sd)) %>%
        select(nsim, arm, n, log_median_adjusted, log_median_sd) %>%
        rename(mu_hat = log_median_adjusted,
               s = log_median_sd) %>%
        pivot_wider(names_from = arm,
                    values_from = c(n, mu_hat, s),
                    names_glue = "{arm}.{.value}")

      cat("Start Bayesian analysis.\n")
      for (k in seq_along(prior_params_list)) {
        prior_params <- prior_params_list[[k]]
        borrow_name <- names(prior_params_list)[k]

        settings <- c(list(true_value = data$true_value, borrowing = borrow_name),
                      data_gen_params, data_gen_params_h,
                      prior_params[!names(prior_params) %in% c("treatment.w", "control.w")])
        settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
          mutate(across(
            .cols = -c(ends_with(".name"), borrowing),
            .fns = as.numeric
          ))

        post <- bayesian_lalonde_decision(endpoint = "continuous",
                                          data_summary = median_est,
                                          prior_params = prior_params,
                                          arm_names = list(treatment = "treatment",
                                                           control = "control",
                                                           control_h = "control_h"),
                                          posterior_infer = FALSE, Lalonde_decision = FALSE,
                                          EXP_TRANSFORM = TRUE)

        post$metrics_post_dist <- calc_post_dist_metrics(endpoint = "g-score",
                                                         true_value = data$true_value$compare_true,
                                                         post$post_est_ci)
        post$settings <- settings

        bayes_results[[length(bayes_results)+1]] <- post
      }
    }
  }
}


#
# metrics_approx_dist_fn <- list.files(pattern = "sim1_metrics")
# bayes_results_fn <- list.files(pattern = "sim1_bayes")
#
# metrics_approx_dist_all <- data.frame()
# for (i in metrics_approx_dist_fn) {
#   load(file = i)
#   metrics_approx_dist_all <- bind_rows(metrics_approx_dist_all, metrics_approx_dist)
#   rm(metrics_approx_dist)
# }
#
# bayes_results_all <- list()
# for (i in bayes_results_fn) {
#   load(file = i)
#   bayes_results_all <- c(bayes_results_all, bayes_results)
#   rm(bayes_results)
# }

saveRDS(bayes_results, file = "sim1_gscore_bayes_results_10_1.rds")
saveRDS(metrics_approx_dist, file = "sim1_gscore_metrics_approx_dist_10_1.rds")

# bayes_results <- readRDS(file = "sim1_gscore_bayes_results.rds")
# metrics_approx_dist <- readRDS(file = "sim1_gscore_metrics_approx_dist.rds")

post_params_all <- data.frame()
post_est_ci_all <- data.frame()
metrics_post_dist_all <- data.frame()
for (i in seq_along(bayes_results)) {
  res <- bayes_results[[i]]
  setting <- res$settings
  post_params <- cbind(setting, res$post_params)
  post_est_ci <- cbind(setting, res$post_est_ci)
  metrics_post_dist <- cbind(setting, res$metrics_post_dist)

  post_params_all <- bind_rows(post_params_all, cbind(i, post_params))
  post_est_ci_all <- bind_rows(post_est_ci_all, cbind(i, post_est_ci))
  metrics_post_dist_all <- bind_rows(metrics_post_dist_all, cbind(i, metrics_post_dist))
}

results_all <- list(
  metrics_approx_dist_all = metrics_approx_dist,
  metrics_post_dist_all = metrics_post_dist_all,
  post_params_all = post_params_all,
  post_est_ci_all = post_est_ci_all
)
saveRDS(results_all, file = "sim1_gscore_results_df_10_1.rds")

# Simulation II -----------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 200, # Sample size for the historical control group.
  ctrl_h_prob = 0.1, # Probability of zero g-scores in the historical control group.
  ctrl_h_mu = 0.01, # Mean for the non-zero g-scores in the historical control group.
  ctrl_h_sigma = 1, # Standard deviation for the non-zero g-scores in the historical control group.
  stringsAsFactors = FALSE
)


param_grid1 <- expand.grid(
  trt_n = seq(15, 40, 5), # Sample size for the treatment group.
  trt_prob = 0.1, # Probability of zero g-scores in the treatment group.
  ctrl_prob = 0.1, # Probability of zero g-scores in the control group.
  ctrl_mu = c(0.01 * c(0.6, 0.8), 0.01, 0.01 / c(0.6, 0.8)),
  trt_sigma = 1, # Standard deviation for the non-zero g-scores in the treatment group.
  ctrl_sigma = 1, # Standard deviation for the non-zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_mu = 0.8*ctrl_mu)

param_grid2 <- expand.grid(
  trt_n = seq(15, 40, 5), # Sample size for the treatment group.
  trt_prob = 0.1, # Probability of zero g-scores in the treatment group.
  ctrl_prob = 0.1, # Probability of zero g-scores in the control group.
  ctrl_mu = c(0.01 * c(0.6, 0.8), 0.01, 0.01 / c(0.6, 0.8)),
  trt_sigma = 1, # Standard deviation for the non-zerot g-scores in the treatment group.
  ctrl_sigma = 1, # Standard deviation for the non-zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_mu = 0.5*ctrl_mu)

param_grid <- bind_rows(param_grid1, param_grid2)

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "g-score")
data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "g-score")

prior_params_list <- list(
  Yes = list(treatment.delta = log(0.85), control.delta = log(0.85), treatment.w = 0, control.w = NULL), # borrowing on the control arm
  No = list(treatment.delta = log(0.85), control.delta = log(0.85), treatment.w = 0, control.w = 0) # no borrowing
)


## Run simulations
nsim = 100000
all_historical_est <- lapply(seq_along(data_gen_params_list_h), function(i) {
  data_gen_params <- data_gen_params_list_h[[i]]
  cat("Running", nsim, "simulations for data_gen_params set", i, "\n\n")
  for (arm in names(data_gen_params)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params[[arm]]$n, "\n")
    cat("\t  prob:", data_gen_params[[arm]]$prob, "\n")
    cat("\t  mu:", data_gen_params[[arm]]$mu, "\n")
    cat("\t  sigma:", data_gen_params[[arm]]$sigma, "\n\n")
  }
  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- sapply(data_gen_params, function(x) x$n)
  prob_list <- sapply(data_gen_params, function(x) x$prob)
  mu_list <- sapply(data_gen_params, function(x) x$mu)
  sigma_list <- sapply(data_gen_params, function(x) x$sigma)

  # Generate simulation dataset
  data <- data_gen_gscore(n_arms = n_arms, arm_names = arm_names,
                          nsim = nsim, n_list = n_list, prob_list = prob_list,
                          mu_list = mu_list, sigma_list = sigma_list)
  median_est <- data$median_est
  return(data.frame(arm = "control_h",
                    n = 200,
                    median_adjusted1 = mean(log(median_est$median_adjusted)),
                    median_adjusted2 = log(mean(median_est$median_adjusted)),
                    median_sd1 = sd(log(median_est$median_adjusted)),
                    median_sd2 = sd(median_est$median_adjusted)/mean(median_est$median_adjusted)))
})

saveRDS(all_historical_est, "all_historical_est.rds")
all_historical_est <- readRDS("all_historical_est.rds")

nsim = 10000
metrics_approx_dist <- data.frame()
bayes_results <- list()
for(j in seq_along(all_historical_est)){
  historical_est <- all_historical_est[[j]]
  data_gen_params_h <- data_gen_params_list_h[[j]]

  for (i in seq_along(data_gen_params_list)) {
    data_gen_params <- data_gen_params_list[[i]]
    cat("Running", nsim, "simulations for data_gen_params set", i, "\n\n")
    for (arm in names(data_gen_params)) {
      cat("\t Arm:", arm, "\n")
      cat("\t  n:", data_gen_params[[arm]]$n, "\n")
      cat("\t  prob:", data_gen_params[[arm]]$prob, "\n")
      cat("\t  mu:", data_gen_params[[arm]]$mu, "\n")
      cat("\t  sigma:", data_gen_params[[arm]]$sigma, "\n\n")
    }
    n_arms <- length(data_gen_params)
    arm_names <- sapply(data_gen_params, function(x) x$name)
    n_list <- sapply(data_gen_params, function(x) x$n)
    prob_list <- sapply(data_gen_params, function(x) x$prob)
    mu_list <- sapply(data_gen_params, function(x) x$mu)
    sigma_list <- sapply(data_gen_params, function(x) x$sigma)


    # Generate simulation dataset
    data <- data_gen_gscore(n_arms = n_arms, arm_names = arm_names,
                            nsim = nsim, n_list = n_list, prob_list = prob_list,
                            mu_list = mu_list, sigma_list = sigma_list)

    settings <- c(list(true_value = data$true_value), data_gen_params, data_gen_params_h)
    settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
      mutate(across(
        .cols = -c(ends_with(".name")),
        .fns = as.numeric
      ))

    # Assess normal approximation of g-score
    cat("Start assessing the normal approximation of log g-score median estimate.\n")
    metrics_approx_dist <- rbind(metrics_approx_dist,
                                 cbind(settings, eval_gscore_approx_dist(data)))

    median_est <- data$median_est
    median_est$log_median_sd <- median_est$median_sd / median_est$median_adjusted
    median_est$log_median_adjusted <- log(median_est$median_adjusted)

    median_est_h <- data.frame(nsim = 1:nsim,
                               arm = "control_h",
                               n = 200,
                               log_median_adjusted = historical_est$median_adjusted2,
                               log_median_sd = historical_est$median_sd2)

    median_est <- bind_rows(median_est, median_est_h) %>%
      mutate(log_median_sd = ifelse(log_median_sd == 0, 0.000001, log_median_sd)) %>%
      select(nsim, arm, n, log_median_adjusted, log_median_sd) %>%
      rename(mu_hat = log_median_adjusted,
             s = log_median_sd) %>%
      pivot_wider(names_from = arm,
                  values_from = c(n, mu_hat, s),
                  names_glue = "{arm}.{.value}")

    cat("Start Bayesian analysis.\n")
    for (k in seq_along(prior_params_list)) {
      prior_params <- prior_params_list[[k]]
      borrow_name <- names(prior_params_list)[k]

      settings <- c(list(true_value = data$true_value, borrowing = borrow_name),
                    data_gen_params, data_gen_params_h,
                    prior_params[!names(prior_params) %in% c("treatment.w", "control.w")])
      settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
        mutate(across(
          .cols = -c(ends_with(".name"), borrowing),
          .fns = as.numeric
        ))

      post <- bayesian_lalonde_decision(endpoint = "continuous",
                                        data_summary = median_est,
                                        prior_params = prior_params,
                                        arm_names = list(treatment = "treatment",
                                                         control = "control",
                                                         control_h = "control_h"),
                                        lrv = 0.8, tv = 0.5,
                                        fgr = 0.2, fsr = 0.1,
                                        posterior_infer = T, Lalonde_decision = T,
                                        EXP_TRANSFORM = TRUE)

      post$metrics_post_dist <- calc_post_dist_metrics(endpoint = "g-score",
                                                       true_value = data$true_value$compare_true,
                                                       post$post_est_ci)
      post$settings <- settings

      bayes_results[[length(bayes_results)+1]] <- post
    }
  }
}



saveRDS(bayes_results, file = "sim2_gscore_bayes_results_10_1.rds")
saveRDS(metrics_approx_dist, file = "sim2_gscore_metrics_approx_dist_10_1.rds")

bayes_results <- readRDS("sim2_gscore_bayes_results_10_1.rds")
metrics_approx_dist <- readRDS("sim2_gscore_metrics_approx_dist_10_1.rds")

post_params_all <- data.frame()
post_est_ci_all <- data.frame()
post_inference_all <- data.frame()
metrics_post_dist_all <- data.frame()
oc_all <- data.frame()
for (i in seq_along(bayes_results)) {
  res <- bayes_results[[i]]
  setting <- res$settings
  post_params <- cbind(setting, res$post_params)
  post_est_ci <- cbind(setting, res$post_est_ci)
  post_inference <- cbind(setting, res$post_inference)
  metrics_post_dist <- cbind(setting, res$metrics_post_dist)
  oc <- cbind(setting, obtain_oc(post_inference))

  post_params_all <- bind_rows(post_params_all, cbind(i, post_params))
  post_est_ci_all <- bind_rows(post_est_ci_all, cbind(i, post_est_ci))
  post_inference_all <- bind_rows(post_inference_all, cbind(i, post_inference))
  metrics_post_dist_all <- bind_rows(metrics_post_dist_all, cbind(i, metrics_post_dist))
  oc_all <- bind_rows(oc_all, cbind(i, oc))
}


bayes_results_all <- list(
  metrics_post_dist_all = metrics_post_dist_all,
  post_params_all = post_params_all,
  post_est_ci_all = post_est_ci_all,
  post_inference_all = post_inference_all,
  oc_all = oc_all
)

saveRDS(bayes_results_all, "sim2_gscore_bayes_results_df_10_1.rds")


# Analysis I --------------------------------------------------------------

results_all <- readRDS(file = "sim1_gscore_results_df.rds")

metrics_approx_dist <- results_all$metrics_approx_dist_all
metrics_post_dist <- results_all$metrics_post_dist_all

metrics_approx_dist$borrowing <- "Frequentist"
metrics_post_dist$borrowing <- ifelse(metrics_post_dist$borrowing == "No", "Borrowing: No", "Borrowing: Yes")

combined_dist <- bind_rows(metrics_approx_dist, metrics_post_dist) %>%
  mutate(across(
    .cols = -c(ends_with(".name"), "borrowing"),
    .fns = as.numeric
  ))

combined_dist$true_value.control_h.true_value <- apply(combined_dist, 1, function(row) {
  get_true_gscore_median(p = as.numeric(row["control_h.prob"]), mu = as.numeric(row["control_h.mu"]), 1)
})
combined_dist$control.median.ratio <- combined_dist$true_value.control.true_value / combined_dist$true_value.control_h.true_value
combined_dist$control.mu = round(exp(combined_dist$control.mu), 3)
combined_dist$treatment.mu = round(exp(combined_dist$treatment.mu), 4)
combined_dist$control_h.mu = round(exp(combined_dist$control_h.mu), 3)
combined_dist$control.n = factor(combined_dist$control.n)
combined_dist$Approach <- combined_dist$borrowing
combined_dist$Approach <- factor(combined_dist$Approach, levels = c("Frequentist", "Borrowing: No", "Borrowing: Yes"))

metrics_df <- combined_dist %>%
  select(treatment.prob, control.prob, control_h.prob,
         true_value.compare_true, control.n, Approach,
         bias_avg_median_ratio1, sd_avg, sd_empirical, cp) %>%
  arrange(Approach) %>%
  arrange(treatment.prob, control.n)

metrics_df_long <- metrics_df %>%
  pivot_wider(names_from = Approach, values_from = c(bias_avg_median_ratio1, sd_avg, sd_empirical, cp))

# kbl(metrics_df[, 4:ncol(metrics_df)], escape = F, row.names = F, digits = 3,
#     format    = "latex",
#     longtable = T,
#     booktabs  = T,
#     col.names = c("True Median Ratio\n(Treatment/Control)<sup>1</sup>", "Sample Size<sup>2</sup>", "Borrowing", "BIAS<sup>3</sup>", "SE<sup>4</sup>", "SD<sup>5</sup>", "CP<sup>6</sup>")) %>%
#   kable_classic(full_width = FALSE) %>%
#   collapse_rows(columns = 1:5, valign = "top") %>%
#   pack_rows(index = c("Case 1" = 27,
#                       "Case 2" = 27,
#                       "Case 3" = 27,
#                       "Case 4" = 27,
#                       "Case 5" = 27)) %>%
#   footnote(number = c("True g-score median of treatment arm over true g-score median of control arm, obtained based on the zero-inflated lognormal distribution (determined by proportion of zeros in both arms)",
#                       "Sample size for each arm.",
#                       "Average bias of the estimated median ratio over 10K simulations.",
#                       "Average standard error of the estimated log median ratio over 10K simulations.",
#                       "Empirical standard deviation of the estimated log median ratio.",
#                       "95% coverage probability under normal distribution."))

metrics_df_long <- metrics_df_long[, c("true_value.compare_true", "control.n",
                             "bias_avg_median_ratio1_Frequentist", "sd_avg_Frequentist",
                             "sd_empirical_Frequentist", "cp_Frequentist",
                             "bias_avg_median_ratio1_Borrowing: No", "sd_avg_Borrowing: No",
                             "sd_empirical_Borrowing: No", "cp_Borrowing: No",
                             "bias_avg_median_ratio1_Borrowing: Yes", "sd_avg_Borrowing: Yes",
                             "sd_empirical_Borrowing: Yes", "cp_Borrowing: Yes")]

kbl(metrics_df_long[,-1], escape = F, row.names = F, digits = 3,
    format    = "latex",
    longtable = T,
    booktabs  = T,
    col.names = c("Sample size",
                  "BIAS", "SE-avg","SD-emp", "CP",
                  "BIAS", "SE-avg","SD-emp", "CP",
                  "BIAS", "SE-avg","SD-emp", "CP")) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 1, "Frequentist" = 4, "Borrowing: No" = 4, "Borrowing: Yes" = 4)) %>%
  pack_rows(index = c("Case 1" = 9,
                      "Case 2" = 9,
                      "Case 3" = 9,
                      "Case 4" = 9,
                      "Case 5" = 9))

metrics_df_long <- metrics_df_long[, c("true_value.compare_true", "control.n",
                                       "bias_avg_median_ratio1_Borrowing: No", "sd_avg_Borrowing: No",
                                       "sd_empirical_Borrowing: No", "cp_Borrowing: No",
                                       "bias_avg_median_ratio1_Borrowing: Yes", "sd_avg_Borrowing: Yes",
                                       "sd_empirical_Borrowing: Yes", "cp_Borrowing: Yes")]

kbl(metrics_df_long[,-1], escape = F, row.names = F, digits = 3,
    format    = "latex",
    booktabs  = T,
    col.names = c("Sample size",
                  "BIAS", "SD-avg","SD-emp", "CP",
                  "BIAS", "SD-avg","SD-emp", "CP")) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 1, "Borrowing: No" = 4, "Borrowing: Yes" = 4)) %>%
  pack_rows(index = c("Case 1" = 9,
                      "Case 2" = 9,
                      "Case 3" = 9,
                      "Case 4" = 9,
                      "Case 5" = 9))


metrics_df$facet <- factor(metrics_df$control.prob,
                           levels = c("0", "0.1", "0.2", "0.3", "0.4"),
                           labels=c(expression(paste(p[t], " = 0.05, ", p[c], " = 0.0, ", p[ch], " = 0.0")),
                                    expression(paste(p[t], " = 0.15, ", p[c], " = 0.1, ", p[ch], " = 0.1")),
                                    expression(paste(p[t], " = 0.25, ", p[c], " = 0.2, ", p[ch], " = 0.2")),
                                    expression(paste(p[t], " = 0.35, ", p[c], " = 0.3, ", p[ch], " = 0.3")),
                                    expression(paste(p[t], " = 0.45, ", p[c], " = 0.4, ", p[ch], " = 0.4"))))

p1 <- ggplot(metrics_df, aes(x = control.n, y = bias_avg_median_ratio1, color = Approach)) +
  geom_point(position = position_dodge(width = 0.3), size = 2) +
  facet_grid(~facet, labeller = labeller(facet = label_parsed)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "BIAS", color= "") +
  theme_bw() +
  theme(legend.position = "bottom")

metrics_long <- metrics_df %>%
  pivot_longer(c(sd_avg, sd_empirical), names_to = "sd_type", values_to = "sd_value")

p2 <- ggplot(metrics_long, aes(x = control.n, y = sd_value, shape = sd_type, color = Approach)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  ggh4x::facet_grid2(Approach~facet,
                     scales = "free_y", independent = "y",
                     labeller = labeller(facet = label_parsed)) +
  labs(x = "Sample size per arm", y = "SD", color= "", shape="") +
  theme_bw() +
  theme(legend.position = "bottom")

p3 <- ggplot(metrics_df, aes(x = control.n, color = Approach, y = cp)) +
  geom_point(position = position_dodge(width = 0.3), size = 2) +
  facet_grid(~facet, labeller = labeller(facet = label_parsed)) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "CP", color = "") +
  theme_bw() +
  theme(legend.position = "bottom")

library(patchwork)
p <- p1 / p3 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')') &
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 15),
        legend.position = "bottom",
        # legend.margin=margin(0,0,0,0),
        # legend.box.margin=margin(-5,-5,-5,0),
        plot.tag = element_text(size = 17)
  )
ggsave(filename = "sim1_gscore_dist_check_metrics.jpg", p, width = 18, height = 6.5)


p <- p1 / p2 / p3 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')') &
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 15),
        legend.position = "bottom",
        # legend.margin=margin(0,0,0,0),
        # legend.box.margin=margin(-5,-5,-5,0),
        plot.tag = element_text(size = 17)
  )
ggsave(filename = "sim1_gscore_dist_check_metrics_1.jpg", p, width = 18, height = 12)



# Analysis II -------------------------------------------------------------

results_all <- readRDS("sim2_gscore_bayes_results_df_10_1.rds")

metrics_post_dist <- results_all$metrics_post_dist_all
metrics_post_dist$borrowing <- ifelse(metrics_post_dist$borrowing == "No", "Borrowing: No", "Borrowing: Yes")

combined_dist <- metrics_post_dist %>%
  mutate(across(
    .cols = -c(ends_with(".name"), "borrowing"),
    .fns = as.numeric
  ))

combined_dist$true_value.control_h.true_value <- apply(combined_dist, 1, function(row) {
  get_true_gscore_median(p = as.numeric(row["control_h.prob"]), mu = as.numeric(row["control_h.mu"]), 1)
})

combined_dist$control.median.ratio <- combined_dist$true_value.control.true_value / combined_dist$true_value.control_h.true_value
combined_dist$control.mu = round(exp(combined_dist$control.mu), 3)
combined_dist$treatment.mu = round(exp(combined_dist$treatment.mu), 4)
combined_dist$control_h.mu = round(exp(combined_dist$control_h.mu), 3)
combined_dist$control.n = factor(combined_dist$control.n)
combined_dist$Approach <- combined_dist$borrowing
combined_dist$Approach <- factor(combined_dist$Approach, levels = c("Borrowing: No", "Borrowing: Yes"))

metrics_df <- combined_dist %>%
  select(treatment.mu, control.mu, control_h.mu,
         true_value.compare_true, control.n, Approach,
         bias_avg_median_ratio1, sd_avg, sd_empirical, cp) %>%
  arrange(Approach) %>%
  arrange(treatment.prob, control.n)

metrics_df_long <- metrics_df %>%
  pivot_wider(names_from = Approach, values_from = c(bias_avg_median_ratio1, sd_avg, sd_empirical, cp))

# kbl(metrics_df[, 4:ncol(metrics_df)], escape = F, row.names = F, digits = 3,
#     format    = "latex",
#     longtable = T,
#     booktabs  = T,
#     col.names = c("True Median Ratio\n(Treatment/Control)<sup>1</sup>", "Sample Size<sup>2</sup>", "Borrowing", "BIAS<sup>3</sup>", "SE<sup>4</sup>", "SD<sup>5</sup>", "CP<sup>6</sup>")) %>%
#   kable_classic(full_width = FALSE) %>%
#   collapse_rows(columns = 1:5, valign = "top") %>%
#   pack_rows(index = c("Case 1" = 27,
#                       "Case 2" = 27,
#                       "Case 3" = 27,
#                       "Case 4" = 27,
#                       "Case 5" = 27)) %>%
#   footnote(number = c("True g-score median of treatment arm over true g-score median of control arm, obtained based on the zero-inflated lognormal distribution (determined by proportion of zeros in both arms)",
#                       "Sample size for each arm.",
#                       "Average bias of the estimated median ratio over 10K simulations.",
#                       "Average standard error of the estimated log median ratio over 10K simulations.",
#                       "Empirical standard deviation of the estimated log median ratio.",
#                       "95% coverage probability under normal distribution."))

metrics_df_long <- metrics_df_long[, c("true_value.compare_true", "control.n",
                                       "bias_avg_median_ratio1_Frequentist", "sd_avg_Frequentist",
                                       "sd_empirical_Frequentist", "cp_Frequentist",
                                       "bias_avg_median_ratio1_Borrowing: No", "sd_avg_Borrowing: No",
                                       "sd_empirical_Borrowing: No", "cp_Borrowing: No",
                                       "bias_avg_median_ratio1_Borrowing: Yes", "sd_avg_Borrowing: Yes",
                                       "sd_empirical_Borrowing: Yes", "cp_Borrowing: Yes")]

kbl(metrics_df_long[,-1], escape = F, row.names = F, digits = 3,
    format    = "latex",
    longtable = T,
    booktabs  = T,
    col.names = c("Sample size",
                  "BIAS", "SE-avg","SD-emp", "CP",
                  "BIAS", "SE-avg","SD-emp", "CP",
                  "BIAS", "SE-avg","SD-emp", "CP")) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 1, "Frequentist" = 4, "Borrowing: No" = 4, "Borrowing: Yes" = 4)) %>%
  pack_rows(index = c("Case 1" = 9,
                      "Case 2" = 9,
                      "Case 3" = 9,
                      "Case 4" = 9,
                      "Case 5" = 9))


metrics_df$facet <- factor(metrics_df$control.prob,
                           levels = c("0", "0.1", "0.2", "0.3", "0.4"),
                           labels=c(expression(paste(theta[t], " = 0.05, ", theta[c], " = 0.0, ", theta[ch], " = 0.0")),
                                    expression(paste(theta[t], " = 0.15, ", theta[c], " = 0.1, ", theta[ch], " = 0.1")),
                                    expression(paste(theta[t], " = 0.25, ", theta[c], " = 0.2, ", theta[ch], " = 0.2")),
                                    expression(paste(theta[t], " = 0.35, ", theta[c], " = 0.3, ", theta[ch], " = 0.3")),
                                    expression(paste(theta[t], " = 0.45, ", theta[c], " = 0.4, ", theta[ch], " = 0.4"))))


p1 <- ggplot(metrics_df, aes(x = control.n, y = bias_avg_median_ratio1, color = Approach)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  facet_grid(~facet, labeller = labeller(facet = label_parsed)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "BIAS", color= "") +
  theme_bw() +
  theme(legend.position = "bottom")

metrics_long <- metrics_df %>%
  pivot_longer(c(sd_avg, sd_empirical), names_to = "sd_type", values_to = "sd_value")

p2 <- ggplot(metrics_long, aes(x = control.n, y = sd_value, shape = sd_type, color = Approach)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  ggh4x::facet_grid2(~facet,
                     # scales = "free_y", independent = "y",
                     labeller = labeller(facet = label_parsed)) +
  labs(x = "Sample size per arm", y = "SD", color= "", shape="") +
  theme_bw() +
  theme(legend.position = "bottom")

p3 <- ggplot(metrics_df, aes(x = control.n, color = Approach, y = cp)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  facet_grid(~facet, labeller = labeller(facet = label_parsed)) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "CP", color = "") +
  theme_bw() +
  theme(legend.position = "bottom")

library(patchwork)
p <- p1 / p2 / p3 +
  plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')') &
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.position = "right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,0),
        plot.tag = element_text(size = 16)
  )
ggsave(filename = "sim1_gscore_dist_check_metrics.jpg", p, width = 20, height = 8)


################################
oc_all <- results_all$oc_all
oc_all <- oc_all[, c("i",
                     "true_value.control.true_value",
                     "true_value.compare_true", "control.n", "borrowing",
                     "treatment.prob", "treatment.mu", "control.prob", "control.mu",
                     "control_h.prob", "control_h.mu", "decision_pr", "proportion_pr")] %>%
  mutate(across(
    .cols = -c(ends_with(".name"), "borrowing", "decision_pr"),
    .fns = as.numeric
  ))

oc_all$true_value.control_h.true_value <- apply(oc_all, 1, function(row) {
  get_true_gscore_median(p = as.numeric(row["control_h.prob"]), mu = as.numeric(row["control_h.mu"]), 1)
})
oc_all$decision_pr <- factor(oc_all$decision_pr, levels = c("no-go", "consider", "go"))
oc_all$true_value.compare_true <- round(oc_all$true_value.compare_true,2)
oc_all$true_value.compare_true = ifelse(oc_all$true_value.compare_true == 0.5, "0.5 (TV)",
                                        ifelse(oc_all$true_value.compare_true == 0.8, "0.8 (LRV)",
                                               oc_all$true_value.compare_true))
oc_all$control.median.ratio <- round(oc_all$true_value.control.true_value / oc_all$true_value.control_h.true_value, 2)
oc_all$control.median.ratio <- factor(oc_all$control.median.ratio)

oc_all$control.mu = round(exp(oc_all$control.mu), 3)
oc_all$treatment.mu = round(exp(oc_all$treatment.mu), 4)
oc_all$control_h.mu = round(exp(oc_all$control_h.mu), 3)
oc_all$control.n = factor(oc_all$control.n)

oc_all$borrowing <- factor(oc_all$borrowing,
                           levels = c("No", "Yes"),
                           labels=c("Borrowing: No", "Borrowing: Yes"))

metrics <- oc_all[oc_all$true_value.compare_true %in% c("0.5 (TV)", "0.8 (LRV)"),]
metrics <- oc_all[(oc_all$true_value.compare_true == "0.8 (LRV)" & oc_all$decision_pr == "go") | (oc_all$true_value.compare_true == "0.5 (TV)" & oc_all$decision_pr == "no-go"),]
false_go_risk <- metrics[metrics$true_value.compare_true == "0.8 (LRV)", c("control.n", "borrowing", "proportion_pr", "control.median.ratio")]
colnames(false_go_risk)[3] <- "Type I Error (FGR)"
false_stop_risk <- metrics[metrics$true_value.compare_true == "0.5 (TV)", c("control.n", "borrowing", "proportion_pr", "control.median.ratio")]
colnames(false_stop_risk)[3] <- "Type II Error (FSR)"
metrics_df <- merge(false_go_risk, false_stop_risk, by = c("control.n", "borrowing", "control.median.ratio"))

metrics_long <- reshape2::melt(metrics_df, id.vars = c("control.n", "borrowing", "control.median.ratio"),
                               measure.vars = c("Type I Error (FGR)", "Type II Error (FSR)"),
                               variable.name = "risk_type", value.name = "risk_value")
metrics_long$control.n <- factor(metrics_long$control.n)

p_risk <- ggplot(metrics_long,
                 aes(x = control.median.ratio, y = risk_value, color = risk_type)) +
  geom_point(size = 2) +
  facet_grid(borrowing ~ control.n,
             labeller = labeller(
               control.n = label_parsed)) +
  labs(y = "Risk Value",
       color = "") +
  xlab(~ paste(theta[c], " / ", theta[ch])) +
  geom_hline(aes(yintercept = 0.2, color = "Type I Error (FGR)"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "Type II Error (FSR)"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,0),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 13)
  )
ggsave("sim2_gscore_conflict_vs_risk_10_1.jpg", width = 10, height = 4)


oc_new <- data.frame()
for (i in seq(1, nrow(oc_all), 3)) {
  tmp <- oc_all[i:(i+2),]
  tmp[tmp$decision_pr == "go", "label_ypos"] <- tmp[tmp$decision_pr == "go", "proportion_pr"]
  tmp[tmp$decision_pr == "consider", "label_ypos"] <- tmp[tmp$decision_pr == "go", "proportion_pr"] + tmp[tmp$decision_pr == "consider", "proportion_pr"]
  tmp[tmp$decision_pr == "no-go", "label_ypos"] <- 1
  oc_new <- rbind(oc_new, tmp)
}

oc_new <- oc_new[, c("true_value.compare_true", "control.n", "borrowing",
                     "control.median.ratio",
                     "decision_pr", "proportion_pr", "label_ypos")]
oc_new$control.median.ratio <- factor(oc_new$control.median.ratio)

myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)
oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)

# p_list <- list()
# for (i in levels(oc_new$control.median.ratio)) {
#   oc_tmp <- oc_new[oc_new$control.median.ratio == i, ]
#   oc_tmp$true_value.compare_true <- factor(oc_tmp$true_value.compare_true ,
#                                            levels = c("0.8 (LRV)", "0.5 (TV)"),
#                                            labels = c(expression(paste(Delta, " = 0.8 (LRV)")),
#                                                       expression(paste(Delta, " = 0.5 (TV)"))))
#   p_list[[i]] <- ggplot(oc_tmp,
#                         aes(x = as.factor(control.n), y = proportion_pr, fill = decision_pr)) +
#     geom_bar(stat = "identity", width = 1) +
#     geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.1)), vjust = 1.6,
#               # fontface = "bold",
#               size = 3
#               ) +
#     scale_fill_manual(name = "Decision", values = myColors)+
#     facet_grid(true_value.compare_true~borrowing,
#                labeller = labeller(true_value.compare_true = label_parsed)) +
#     ggtitle(bquote(theta[c] / theta[ch] == .(i))) +
#     labs(x = "Sample size per arm", y = "Probability") +
#     theme(legend.position = "bottom",
#           title = element_text(size = 12),
#           legend.text = element_text(size = 10),
#           axis.text = element_text(size = 9),
#           strip.text = element_text(size = 12)) +
#     theme_bw()
# }
# library(patchwork)
# p_oc <- p_list[[1]] + p_list[[2]] + p_list[[3]] + p_list[[4]] + p_list[[5]] +
#   plot_layout(ncol = 5, guides = "collect") &
#   theme(legend.position = "bottom",
#         title = element_text(size = 14),
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.margin=margin(0,0,0,0),
#         legend.box.margin=margin(-5,-5,-5,0),
#         axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 13))


p_list <- list()
for (i in levels(oc_new$true_value.compare_true)) {
  oc_tmp <- oc_new[oc_new$true_value.compare_true == i, ]
  oc_tmp$control.median.ratio <- factor(oc_tmp$control.median.ratio,
                                     levels = c("0.6", "0.8", "1", "1.25", "1.67"),
                                     labels=c(expression(paste(theta[c], " / ", theta[ch], " = 0.6")),
                                              expression(paste(theta[c], " / ", theta[ch], " = 0.8")),
                                              expression(paste(theta[c], " / ", theta[ch], " = 1")),
                                              expression(paste(theta[c], " / ", theta[ch], " = 1.25")),
                                              expression(paste(theta[c], " / ", theta[ch], " = 1.67"))))
  p_list[[i]] <- ggplot(oc_tmp,
                        aes(x = as.factor(control.n), y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.1)), vjust = 1.6,
              # fontface = "bold",
              size = 3.2) +
    scale_fill_manual(name = "Decision", values = myColors)+
    facet_grid(control.median.ratio~borrowing,
               labeller = labeller(control.median.ratio = label_parsed)) +
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
p_oc <- p_list[[1]] + p_list[[2]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        title = element_text(size = 15),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,0),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13))
ggsave("sim2_gscore_zone_size_10_1.jpg",p_oc, width = 14, height = 8)


# p_list <- list()
# for (i in levels(oc_new$borrowing)) {
#   oc_tmp <- oc_new[oc_new$borrowing == i, ]
#   oc_tmp$control.prob.diff <- factor(oc_tmp$control.prob.diff,
#                                       levels = c("-0.2", "-0.1", "0", "0.1", "0.2"),
#                                       labels=c(expression(paste(theta[c], " - ", theta[ch], " = -0.2")),
#                                                expression(paste(theta[c], " - ", theta[ch], " = -0.1")),
#                                                expression(paste(theta[c], " - ", theta[ch], " = 0")),
#                                                expression(paste(theta[c], " - ", theta[ch], " = 0.1")),
#                                                expression(paste(theta[c], " - ", theta[ch], " = 0.2"))))
#   oc_tmp$true_value.compare_true <- factor(oc_tmp$true_value.compare_true ,
#                                            levels = c("0.1 (LRV)", "0.15", "0.2 (TV)"),
#                                            labels = c(expression(paste(Delta, " = 0.1 (LRV)")),
#                                                       expression(paste(Delta, " = 0.15")),
#                                                       expression(paste(Delta, " = 0.2 (TV)"))))
#
#   p_list[[i]] <- ggplot(oc_tmp,
#                         aes(x = as.factor(control.n), y = proportion_pr, color = decision_pr)) +
#     geom_point() +
#     scale_color_manual(name = "Decision", values = myColors) +
#     facet_grid(true_value.compare_true ~ control.prob.diff,
#                labeller = labeller(control.prob.diff = label_parsed,
#                                    true_value.compare_true = label_parsed)) +
#     ggtitle(i) +
#     labs(x = "Number of patients per arm", y = "Probability") +
#     theme(legend.position = "bottom",
#           title = element_text(size = 12),
#           legend.text = element_text(size = 12),
#           axis.text = element_text(size = 9),
#           strip.text = element_text(size = 13)) +
#     theme_bw()
# }
#
# p_list[[1]] + p_list[[2]] +
#   plot_layout(ncol = 2, guides = "collect") &
#   theme(legend.position = "bottom",
#         title = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         axis.text = element_text(size = 9),
#         strip.text.y = element_text(size = 10),
#         strip.text.x = element_text(size = 11) )

oc2 <- oc_new %>%
  select(control.median.ratio, control.n, true_value.compare_true, borrowing, decision_pr, proportion_pr) %>%
  mutate(proportion_pr = scales::percent(proportion_pr, 0.1)) %>%
  arrange(decision_pr) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = c(decision_pr, borrowing), values_from = proportion_pr) %>%
  arrange(control.n) %>%
  arrange(control.median.ratio)

oc2$`no-go_Borrowing: No` = cell_spec(oc2$`no-go_Borrowing: No`,
                                      background = ifelse(oc2$true_value.compare_true == "0.2 (TV)", "steelblue2", "white"))

oc2$`no-go_Borrowing: Yes` = cell_spec(oc2$`no-go_Borrowing: Yes`,
                                       background = ifelse(oc2$true_value.compare_true == "0.2 (TV)", "steelblue2", "white"))

oc2$`go_Borrowing: No` = cell_spec(oc2$`go_Borrowing: No`,
                                   background = ifelse(oc2$true_value.compare_true == "0.1 (LRV)", "pink", "white"))

oc2$`go_Borrowing: Yes` = cell_spec(oc2$`go_Borrowing: Yes`,
                                    background = ifelse(oc2$true_value.compare_true == "0.1 (LRV)", "pink", "white"))

kbl(oc2[,-1], escape = T, row.names = F, digits = 3,
    format    = "latex",
    booktabs  = T,
    col.names = c("n", "Delta", rep(c("Pr(No-Go)", "Pr(Consider)","Pr(Go)"), 2))) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 2, "Borrowing: No" = 3, "Borrowing: Yes" = 3)) %>%
  pack_rows(index = c("Case 1:" = 12,
                      "Case 2:" = 12,
                      "Case 3:" = 12,
                      "Case 4:" = 12,
                      "Case 5:" = 12))










