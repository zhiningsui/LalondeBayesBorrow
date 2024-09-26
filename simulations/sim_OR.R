library(kableExtra)


# Simulation II -----------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 200, # Sample size for the historical control group.
  ctrl_h_prob = 0.3, # Probability of zero g-scores in the historical control group.
  stringsAsFactors = FALSE
)

param_grid <- expand.grid(
  trt_n = c(20, 40, 60), # Sample size for the treatment group.
  ctrl_prob = seq(0.1, 0.5, 0.1), # Probability of zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_prob = 0.2+ctrl_prob)


data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                               create_data_gen_params, endpoint = "OR")

data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "OR")

prior_params_list <- list(
  Yes = list(treatment.delta = 0.1, control.delta = 0.1, treatment.a = 1, treatment.b = 1, control.a = 1, control.b = 1, treatment.w = 0, control.w = NULL, gate = 0.2), # borrowing on the control arm
  No = list(treatment.delta = 0.1, control.delta = 0.1, treatment.a = 1, treatment.b = 1, control.a = 1, control.b = 1, treatment.w = 0, control.w = 0) # no borrowing
)


nsim = 10000
bayes_results <- list()
historical_est <- all_historical_est[[1]]
data_gen_params_h <- data_gen_params_list_h[[1]]

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
  data <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list)

  # settings <- c(list(true_value = data$true_value), data_gen_params, data_gen_params_h)
  # settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
  #   mutate(across(
  #     .cols = -c(ends_with(".name")),
  #     .fns = as.numeric
  #   ))

  freq <- data$freq
  freq_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                       control_h.count = data_gen_params_h$control_h$n*data_gen_params_h$control_h$prob)
  freq <- cbind(freq, freq_h)
  # freq <- cbind(settings, freq)

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

    post <- bayesian_lalonde_decision(endpoint = "OR",
                                      data_summary = freq,
                                      prior_params = prior_params,
                                      arm_names = list(treatment = "treatment",
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

    bayes_results[[length(bayes_results)+1]] <- post
  }
}


saveRDS(bayes_results, file = "sim_bayes_results_OR.rds")

bayes_results <- readRDS("sim2_bayes_results_2.rds")
metrics_approx_dist <- readRDS("sim2_metrics_approx_dist_2.rds")

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

saveRDS(bayes_results_all, "sim2_bayes_results_all_2.rds")




# Analysis I --------------------------------------------------------------


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
saveRDS(results_all, file = "sim1_results_all_2.rds")

results_all <- readRDS(file = "sim1_results_all_2.rds")
metrics_approx_dist <- results_all$metrics_approx_dist_all
metrics_post_dist <- results_all$metrics_post_dist_all



# metrics_approx_dist$borrowing <- "Frequentist"
# metrics_post_dist$borrowing <- ifelse(metrics_post_dist$borrowing == "No", "Bayesian: No Borrowing", "Bayesian: Borrowing")
# combined_dist <- bind_rows(metrics_approx_dist, metrics_post_dist)
#
# combined_dist$control.mu = round(exp(as.numeric(combined_dist$control.mu)), 3)
# combined_dist$treatment.mu = round(exp(as.numeric(combined_dist$treatment.mu)), 4)
# combined_dist$control_h.mu = round(exp(as.numeric(combined_dist$control_h.mu)), 3)
# combined_dist$control.n = factor(combined_dist$control.n)
# combined_dist$Approach <- combined_dist$borrowing
# combined_dist$Approach <- factor(combined_dist$Approach, levels = c("Frequentist", "Bayesian: No Borrowing", "Bayesian: Borrowing"))
#
#
# # metrics_df <- combined_dist %>%
# #   select(treatment.mu, control.mu, control_h.mu, treatment.prob,
# #          control.prob, control_h.prob, true_value.compare_true, control.n, Approach,
# #          bias_avg_median_ratio1, sd_avg, sd_empirical, cp) %>%
# #   arrange(treatment.mu, control.mu, treatment.prob, control.n)
# #
# # metrics_df <- metrics_df[metrics_df$control.prob == metrics_df$control_h.prob,]
# #
# # kbl(metrics_df, escape = F, row.names = F,
# #     format    = "html",
# #     longtable = T,
# #     booktabs  = T,
# #     col.names = c("Treatment", "Control", "Historical", "Treatment", "Control","Historical", "True Median Ratio\n(Treatment/Control)<sup>1</sup>", "Sample Size<sup>2</sup>", "Approach", "BIAS<sup>3</sup>", "SE<sup>4</sup>", "SD<sup>5</sup>", "CP<sup>6</sup>")) %>%
# #   kable_classic(full_width = FALSE) %>%
# #   add_header_above(c("Mean of positive g-score median" = 3, "Prob of Zero" = 3, " " = 7)) %>%
# #   collapse_rows(columns = 1:8, valign = "top") %>%
# #   footnote(number = c("True g-score median of treatment arm over true g-score median of control arm, obtained based on the zero-inflated lognormal distribution (determined by proportion of zeros in both arms)",
# #                       "Sample size for each arm.",
# #                       "Average bias of the estimated median ratio over 10K simulations.",
# #                       "Average standard error of the estimated log median ratio over 10K simulations.",
# #                       "Empirical standard deviation of the estimated log median ratio.",
# #                       "95% coverage probability under normal distribution."))
#
# metrics_df <- combined_dist %>%
#   select(treatment.prob, control.prob, control_h.prob,
#          true_value.compare_true, control.n, Approach,
#          bias_avg_median_ratio1, sd_avg, sd_empirical, cp) %>%
#   arrange(treatment.prob, control.n)
#
# kbl(metrics_df, escape = F, row.names = F,
#     format    = "html",
#     longtable = T,
#     booktabs  = T,
#     col.names = c("Treatment", "Control","Historical", "True Median Ratio\n(Treatment/Control)<sup>1</sup>", "Sample Size<sup>2</sup>", "Approach", "BIAS<sup>3</sup>", "SE<sup>4</sup>", "SD<sup>5</sup>", "CP<sup>6</sup>")) %>%
#   kable_classic(full_width = FALSE) %>%
#   add_header_above(c("Prob of Zero" = 3, " " = 7)) %>%
#   collapse_rows(columns = 1:5, valign = "top") %>%
#   footnote(number = c("True g-score median of treatment arm over true g-score median of control arm, obtained based on the zero-inflated lognormal distribution (determined by proportion of zeros in both arms)",
#                       "Sample size for each arm.",
#                       "Average bias of the estimated median ratio over 10K simulations.",
#                       "Average standard error of the estimated log median ratio over 10K simulations.",
#                       "Empirical standard deviation of the estimated log median ratio.",
#                       "95% coverage probability under normal distribution.")) %>%
#   save_kable("sim1_metrics_table.png")
#
# metrics_df$probs <- paste0("Treatment Prob of Zero = ", metrics_df$treatment.prob, ",\nControl Prob of Zero = ", metrics_df$control.prob)
# p1 <- ggplot(metrics_df, aes(x = control.n, y = bias_avg_median_ratio1, color = Approach)) +
#   geom_point() +
#   facet_grid(~probs) +
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
#   labs(x = "Sample Size", y = "Bias", title = "Average Bias") +
#   theme_bw() +
#   theme(legend.position = "bottom")
#
#
# metrics_long <- metrics_df %>%
#   pivot_longer(c(sd_avg, sd_empirical), names_to = "sd_type", values_to = "sd_value")
#
# p2 <- ggplot(metrics_long, aes(x = control.n, y = sd_value, color = sd_type)) +
#   geom_point() +
#   facet_grid(Approach~probs) +
#   labs(x = "Sample Size", y = "SD", title = "Average SD vs. Empirical SD") +
#   theme_bw() +
#   theme(legend.position = "bottom")
#
# p3 <- ggplot(metrics_df, aes(x = control.n, color = Approach, y = cp)) +
#   geom_point() +
#   facet_grid(~probs) +
#   geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
#   labs(x = "Sample Size", y = "Coverage probability", title = "Coverage Probability") +
#   theme_bw() +
#   theme(legend.position = "bottom")
#
#
# p <- p1 / p2 / p3 +
#   plot_layout(heights = c(1,3,1))
#
# ggsave(filename = "sim1_dist_check_metrics.jpg", p, width = 12, height = 12)
#
#
#
#

metrics_post_dist$control.mu = round(exp(as.numeric(metrics_post_dist$control.mu)), 3)
metrics_post_dist$treatment.mu = round(exp(as.numeric(metrics_post_dist$treatment.mu)), 4)
metrics_post_dist$control_h.mu = round(exp(as.numeric(metrics_post_dist$control_h.mu)), 3)
metrics_post_dist$control.n = factor(metrics_post_dist$control.n)
metrics_post_dist$Approach <- metrics_post_dist$borrowing
metrics_post_dist$Approach <- ifelse(metrics_post_dist$borrowing == "No", "Yes", "No")


metrics_df <- metrics_post_dist %>%
  select(treatment.prob, control.prob, control_h.prob,
         true_value.compare_true, control.n, Approach,
         bias_avg_median_ratio1, sd_avg, sd_empirical, cp) %>%
  arrange(treatment.prob, control.n)

kbl(metrics_df, escape = F, row.names = F,
    format    = "html",
    longtable = T,
    booktabs  = T,
    col.names = c("Treatment", "Control","Historical", "True Median Ratio\n(Treatment/Control)<sup>1</sup>", "Sample Size<sup>2</sup>", "Borrowing", "BIAS<sup>3</sup>", "SE<sup>4</sup>", "SD<sup>5</sup>", "CP<sup>6</sup>")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c("Prob of Zero" = 3, " " = 7)) %>%
  collapse_rows(columns = 1:5, valign = "top") %>%
  footnote(number = c("True g-score median of treatment arm over true g-score median of control arm, obtained based on the zero-inflated lognormal distribution (determined by proportion of zeros in both arms)",
                      "Sample size for each arm.",
                      "Average bias of the estimated median ratio over 10K simulations.",
                      "Average standard error of the estimated log median ratio over 10K simulations.",
                      "Empirical standard deviation of the estimated log median ratio.",
                      "95% coverage probability under normal distribution.")) %>%
  save_kable("sim1_metrics_table.png")

metrics_df$probs <- paste0("Treatment Prob of Zero = ", metrics_df$treatment.prob, ",\nControl Prob of Zero = ", metrics_df$control.prob)
p1 <- ggplot(metrics_df, aes(x = control.n, y = bias_avg_median_ratio1, color = Approach)) +
  geom_point() +
  facet_grid(~probs) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Sample Size", y = "Bias", title = "Average Bias", color= "Borrowing") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(filename = "sim1_check_post_dist_bias.jpg", p1, width = 10, height = 3.6)

metrics_post_dist$Approach <- ifelse(metrics_post_dist$borrowing == "No", "No", "Yes")
metrics_df <- metrics_post_dist %>%
  select(treatment.prob, control.prob, control_h.prob,
         true_value.compare_true, control.n, Approach,
         bias_avg_median_ratio1, sd_avg, sd_empirical, cp) %>%
  arrange(treatment.prob, control.n)
metrics_df$probs <- paste0("Treatment Prob of Zero = ", metrics_df$treatment.prob, ",\nControl Prob of Zero = ", metrics_df$control.prob)
metrics_long <- metrics_df %>%
  pivot_longer(c(sd_avg, sd_empirical), names_to = "sd_type", values_to = "sd_value")

p2 <- ggplot(metrics_long, aes(x = control.n, y = sd_value, color = sd_type)) +
  geom_point() +
  ggh4x::facet_grid2(Approach~probs,
                     scales = "free_y", independent = "y",
             labeller = labeller(
               Approach = function(x) paste("Borrowing: ", x))) +
  labs(x = "Sample Size", y = "SD", title = "Average SD vs. Empirical SD") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(filename = "sim1_check_post_dist_sd.jpg", p2, width = 12, height = 4)

p3 <- ggplot(metrics_df, aes(x = control.n, color = Approach, y = cp)) +
  geom_point() +
  facet_grid(~probs) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  labs(x = "Sample Size", y = "Coverage probability", color = "Borrowing",
       title = "Coverage Probability") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(filename = "sim1_check_post_dist_cp.jpg", p3, width = 10, height = 4)


p <- p1 / p2 / p3 +
  plot_layout(heights = c(1,3,1))

ggsave(filename = "sim1_dist_check_metrics.jpg", p, width = 12, height = 12)







metrics_post_dist <- metrics_post_dist %>%
  mutate(across(
    .cols = -c(ends_with(".name"), "borrowing"),
    .fns = as.numeric
  ))

metrics_post_dist$control.mu = round(exp(as.numeric(metrics_post_dist$control.mu)), 3)
metrics_post_dist$treatment.mu = round(exp(as.numeric(metrics_post_dist$treatment.mu)), 4)
metrics_post_dist$control_h.mu = round(exp(as.numeric(metrics_post_dist$control_h.mu)), 3)
metrics_post_dist$control.n = factor(metrics_post_dist$control.n)

metrics_df <- metrics_post_dist %>%
  select(treatment.mu, control.mu, control_h.mu, treatment.prob,
         control.prob, control_h.prob, true_value.compare_true, control.n, borrowing,
         bias_avg_median_ratio1, sd_avg, sd_empirical, cp) %>%
  arrange(treatment.mu, control.mu, treatment.prob, control.n)

kbl(metrics_df, escape = F, row.names = F,
    format    = "html",
    longtable = T,
    booktabs  = T,
    col.names = c("Treatment", "Control", "Historical", "Treatment", "Control","Historical", "True Median Ratio\n(Treatment/Control)<sup>1</sup>", "Sample Size<sup>2</sup>", "Borrowing", "BIAS<sup>3</sup>", "SE<sup>4</sup>", "SD<sup>5</sup>", "CP<sup>6</sup>")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c("Mean of positive g-score median" = 3, "Prob of Zero" = 3, " " = 7)) %>%
  collapse_rows(columns = 1:8, valign = "top") %>%
  footnote(number = c("True g-score median of treatment arm over true g-score median of control arm, obtained based on the zero-inflated lognormal distribution (determined by proportion of zeros in both arms)",
                      "Sample size for each arm.",
                      "Average bias of the estimated median ratio over 10K simulations.",
                      "Average standard error of the estimated log median ratio over 10K simulations.",
                      "Empirical standard deviation of the estimated log median ratio.",
                      "95% coverage probability under normal distribution."))


p1 <- ggplot(metrics_df, aes(x = control.n, y = bias_avg_median_ratio1, color = borrowing)) +
  geom_point() +
  ggh4x::facet_grid2(~treatment.prob + control.prob, labeller = label_both,
                     scales = "free_y", independent = "y") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Sample Size", y = "Bias", title = "Average Bias") +
  theme_bw()

metrics_long <- metrics_df %>%
  pivot_longer(c(sd_avg, sd_empirical), names_to = "sd_type", values_to = "sd_value")

p2 <- ggplot(metrics_long, aes(x = control.n, y = sd_value, color = borrowing, shape = sd_type)) +
  geom_point() +
  facet_grid(~treatment.prob + control.prob, labeller = label_both) +
  labs(x = "Sample Size", y = "SD", title = "Average SD vs. Empirical SD") +
  theme_bw() +
  theme(legend.position = "bottom")

p3 <- ggplot(metrics_df, aes(x = control.n, color = borrowing, y = cp)) +
  geom_point() +
  facet_grid(~treatment.prob + control.prob, labeller = label_both) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  labs(x = "Sample Size", y = "Coverage probability", title = "Coverage Probability") +
  theme_bw()

p <- p1 / p2 / p3

ggsave(filename = "sim1_dist_check_post_metrics.jpg", p, width = 8, height = 8)




# Analysis II -------------------------------------------------------------

bayes_results_all <- readRDS("sim2_bayes_results_all_2.rds")
oc_all <- bayes_results_all$oc_all

oc_all <- oc_all[, c("i", "true_value.compare_true", "control.n", "borrowing",
                     "treatment.prob", "treatment.mu", "control.prob", "control.mu",
                     "control_h.prob", "control_h.mu", "decision_pr", "proportion_pr")]

oc_all$decision_pr <- factor(oc_all$decision_pr, levels = c("no-go", "consider", "go"))
oc_all$true_value.compare_true <- round(as.numeric(oc_all$true_value.compare_true),2)
oc_all$true_value.compare_true = ifelse(oc_all$true_value.compare_true == 0.5, "0.5 (TV)",
                                        ifelse(oc_all$true_value.compare_true == 0.8, "0.8 (LRV)",
                                               oc_all$true_value.compare_true))
oc_all$control.mu <- as.numeric(oc_all$control.mu)
oc_all$control_h.mu <- as.numeric(oc_all$control_h.mu)
oc_all$control.mu.logdiff <- as.character(round(exp(oc_all$control.mu - oc_all$control_h.mu), 2))
oc_all$control.mu.diff <- paste0("log(", as.character(oc_all$control.mu.logdiff), ")")

oc_all$borrowing <- factor(oc_all$borrowing,
                           levels = c("No", "Yes"),
                           labels=c("Borrowing: No", "Borrowing: Yes"))
oc_all$control.mu.logdiff <- as.factor(oc_all$control.mu.logdiff)

metrics <- oc_all[oc_all$true_value.compare_true %in% c("0.5 (TV)", "0.8 (LRV)"),]
metrics <- oc_all[(oc_all$true_value.compare_true == "0.8 (LRV)" & oc_all$decision_pr == "go") | (oc_all$true_value.compare_true == "0.5 (TV)" & oc_all$decision_pr == "no-go"),]

false_go_risk <- metrics[metrics$true_value.compare_true == "0.8 (LRV)", c("control.n", "borrowing", "proportion_pr", "control.mu.logdiff")]
colnames(false_go_risk)[3] <- "false_go_risk"
false_stop_risk <- metrics[metrics$true_value.compare_true == "0.5 (TV)", c("control.n", "borrowing", "proportion_pr", "control.mu.logdiff")]
colnames(false_stop_risk)[3] <- "false_stop_risk"
metrics_df <- merge(false_go_risk, false_stop_risk, by = c("control.n", "borrowing", "control.mu.logdiff"))

metrics_long <- reshape2::melt(metrics_df, id.vars = c("control.n", "borrowing", "control.mu.logdiff"),
                     measure.vars = c("false_go_risk", "false_stop_risk"),
                     variable.name = "risk_type", value.name = "risk_value")

ggplot(metrics_long,
       aes(x = control.mu.logdiff, y = risk_value, color = risk_type)) +
  geom_point() +
  facet_grid(borrowing ~ control.n,
             labeller = labeller(
               control.n = function(x) paste("n = ", x))) +
  labs(title = "FGR and FSR vs. Prior-data Conflict",
       y = "Risk Value",
       color = "Risk Type") +
  xlab(~ paste(theta[c], " / ", theta[ch])) +
  geom_hline(aes(yintercept = 0.2, color = "false_go_risk"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "false_stop_risk"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom", strip.text.x = element_text(size = 10))

ggsave("sim2_pd_conflict_vs_risk.jpg", width = 9, height = 3.5)

ggplot(metrics_long[metrics_long$borrowing == "Borrowing: Yes",],
       aes(x = control.mu.logdiff, y = risk_value, color = risk_type)) +
  geom_point() +
  facet_grid(borrowing ~ control.n,
             labeller = labeller(
               control.n = function(x) paste("n = ", x))) +
  labs(title = "FGR and FSR vs. Prior-data Conflict",
       y = "Risk Value",
       color = "Risk Type") +
  xlab(~ paste(theta[c], " / ", theta[ch])) +
  geom_hline(aes(yintercept = 0.2, color = "false_go_risk"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "false_stop_risk"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom", strip.text.x = element_text(size = 10))

ggsave("sim2_pd_conflict_vs_risk_borrowing.jpg", width = 9, height = 2.5)


oc_new <- data.frame()
for (i in seq(1, nrow(oc_all), 3)) {
  tmp <- oc_all[i:(i+2),]
  tmp[tmp$decision_pr == "go", "label_ypos"] <- tmp[tmp$decision_pr == "go", "proportion_pr"]
  tmp[tmp$decision_pr == "consider", "label_ypos"] <- tmp[tmp$decision_pr == "go", "proportion_pr"] + tmp[tmp$decision_pr == "consider", "proportion_pr"]
  tmp[tmp$decision_pr == "no-go", "label_ypos"] <- 1
  oc_new <- rbind(oc_new, tmp)
}


# oc_new <- oc_new[oc_new$control.mu.logdiff == "1", c("true_value.compare_true", "control.n", "borrowing", "control.prob", "control.mu.logdiff",
#                                                      "decision_pr", "proportion_pr", "label_ypos")]

oc_new <- oc_new[, c("true_value.compare_true", "control.n", "borrowing",
                     "control.prob", "control.mu.logdiff",
                     "decision_pr", "proportion_pr", "label_ypos")]

myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)

lines <- data.frame(true_value.compare_true = c("0.5 (TV)", "0.8 (LRV)"),
                    Z = c(0.9, 0.2))
p_list <- list()
for (i in levels(oc_new$control.mu.logdiff)) {
  p_list[[i]] <- ggplot(oc_new[oc_new$control.mu.logdiff == i, ],
         aes(x = as.factor(control.n), y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    # geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.1)), vjust = 1.6,
    #           size = 1.2, fontface = "bold") +
    scale_fill_manual(name = "Decision", values = myColors)+
    facet_grid(true_value.compare_true~borrowing,
               labeller = labeller(
                 true_value.compare_true = function(x) paste("True median ratio =", x))) +
    geom_hline(data = lines, aes(yintercept = Z), linetype = "dashed") +
    ggtitle(bquote(theta[c] / theta[ch] == .(i))) +
    labs(x = "Number of patients per arm", y = "Probability") +
    theme(legend.position = "bottom",
          title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 12)) +
    theme_bw()
}
library(patchwork)
p_list[[1]] + p_list[[2]] + p_list[[3]] + p_list[[4]] + p_list[[5]] +
  plot_annotation(title = "Probability of Each Decision Under Various Sample Sizes") +
  plot_layout(ncol = 5, guides = "collect") &
  theme(legend.position = "bottom",
        title = element_text(size = 10),
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 9),
        strip.text.y = element_text(size = 9),
        strip.text.x = element_text(size = 9) )
ggsave("sim2_zone_size.jpg", width = 15, height = 5.5)


p_list <- list()
for (i in levels(oc_new$borrowing)) {
  oc_tmp <- oc_new[oc_new$borrowing == i, ]
  oc_tmp$control.mu.logdiff <- factor(oc_tmp$control.mu.logdiff,
                                      levels = c("0.6", "0.8", "1", "1.25", "1.67"),
                                      labels=c(expression(paste(theta[c], "/", theta[ch], " = 0.6")),
                                               expression(paste(theta[c], "/", theta[ch], " = 0.8")),
                                               expression(paste(theta[c], "/", theta[ch], " = 1")),
                                               expression(paste(theta[c], "/", theta[ch], " = 1.25")),
                                               expression(paste(theta[c], "/", theta[ch], " = 1.67"))))
  oc_tmp$true_value.compare_true <- paste0("True median ratio = ", oc_tmp$true_value.compare_true)
  p_list[[i]] <- ggplot(oc_tmp,
                        aes(x = as.factor(control.n), y = proportion_pr, color = decision_pr)) +
    geom_point() +
    scale_color_manual(name = "Decision", values = myColors) +
    facet_grid(true_value.compare_true ~ control.mu.logdiff,
               labeller = labeller(control.mu.logdiff = label_parsed)) +
    ggtitle(i) +
    labs(x = "Number of patients per arm", y = "Probability") +
    theme(legend.position = "bottom",
          title = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 9),
          strip.text = element_text(size = 13)) +
    theme_bw()
}

p_list[[1]] + p_list[[2]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 9),
        strip.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 11) )

ggsave("sim2_zone_size_2.jpg", width = 12.5, height = 5.8)

oc2 <- oc_new %>%
  select(-c(label_ypos, control.prob)) %>%
  mutate(proportion_pr = scales::percent(proportion_pr, 0.1)) %>%
  arrange(decision_pr) %>%
  arrange(control.mu.logdiff) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = c(decision_pr, borrowing), values_from = proportion_pr) %>%
  arrange(control.n)

oc2$`no-go_Borrowing: No` = cell_spec(oc2$`no-go_Borrowing: No`,
                           background = ifelse(oc2$true_value.compare_true == "0.5 (TV)", "steelblue2", "white"))

oc2$`no-go_Borrowing: Yes` = cell_spec(oc2$`no-go_Borrowing: Yes`,
                            background = ifelse(oc2$true_value.compare_true == "0.5 (TV)", "steelblue2", "white"))

oc2$`go_Borrowing: No` = cell_spec(oc2$`go_Borrowing: No`,
                        background = ifelse(oc2$true_value.compare_true == "0.8 (LRV)", "pink", "white"))

oc2$`go_Borrowing: Yes` = cell_spec(oc2$`go_Borrowing: Yes`,
                         background = ifelse(oc2$true_value.compare_true == "0.8 (LRV)", "pink", "white"))

for (i in levels(oc2$control.mu.logdiff)) {
  kbl(oc2[oc2$control.mu.logdiff == i, c(2,1,4:ncol(oc2))], escape = F,
      col.names = c("Sample Size", "True value",
                    rep(c("No-Go", "Consider","Go"), 2))) %>%
    kable_classic(full_width = FALSE) %>%
    add_header_above(c(" " = 2, "Borrowing: No" = 3, "Borrowing: Yes" = 3)) %>%
    column_spec(1, bold = T) %>%
    collapse_rows(columns = 1:2, valign = "top")  %>%
    save_kable(paste0("oc_", i, ".png"),
               zoom = 1.5)

}









