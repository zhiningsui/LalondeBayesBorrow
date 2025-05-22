library(kableExtra)
library(LalondeBayesBorrow)
library(dplyr)
library(tidyverse)

# Simulation ------------------------------------------------------------

# Define historical control data
param_grid_h <- expand.grid(
  ctrl_h_n = 180,
  ctrl_h_p = 0.04,
  ctrl_h_mu = -0.2,
  ctrl_h_sigma = 0.3,
  stringsAsFactors = FALSE
)

# Define current trial configurations
param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_mu = seq(-0.35, -0.05, 0.05),
  trt_sigma = 0.3,
  ctrl_sigma = 0.3,
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_mu = ctrl_mu - 0.1)

param_grid2 <- expand.grid(
  trt_n = c(45),
  ctrl_mu = seq(-0.35, -0.05, 0.05),
  trt_sigma = 0.4,
  ctrl_sigma = 0.4,
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_mu = ctrl_mu - 0.2)

param_grid <- bind_rows(param_grid1, param_grid2)
param_grid$ctrl_p <- -0.2 * param_grid$ctrl_mu
param_grid$trt_p <- -0.2 * param_grid$trt_mu

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "DpR")

data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "DpR")

data_gen_params_h <- data_gen_params_list_h[[1]]

## Run simulations
nsim = 10000
bayes_results <- list()
for (i in seq_along(data_gen_params_list)) {
  start_time_i <- Sys.time()

  data_gen_params <- data_gen_params_list[[i]]

  cat("Running", nsim, "simulations for data_gen_params set", i, "\n\n")

  for (arm in names(data_gen_params)) {
    cat("\t Arm:", arm, "\n")
    cat("\t  n:", data_gen_params[[arm]]$n, "\n")
    cat("\t  p:", data_gen_params[[arm]]$p, "\n")
    cat("\t  mu:", data_gen_params[[arm]]$mu, "\n")
    cat("\t  sigma:", data_gen_params[[arm]]$sigma, "\n\n")
  }

  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- sapply(data_gen_params, function(x) x$n)
  p_list <- sapply(data_gen_params, function(x) x$p)
  mu_list <- sapply(data_gen_params, function(x) x$mu)
  sigma_list <- sapply(data_gen_params, function(x) x$sigma)

  # Generate simulation dataset
  data <- data_gen_DpR(n_arms = n_arms, arm_names = arm_names,
                       nsim = nsim, n_list = n_list, p_list = p_list,
                       mu_list = mu_list, sigma_list = sigma_list)

  hist(data$data$value[data$data$arm == "control"])
  hist(data$data$value[data$data$arm == "treatment"])

  settings <- c(list(true_value = data$true_value), data_gen_params, data_gen_params_h)
  settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
    mutate(
      across(
        .cols = -c(ends_with(".name")),
        .fns = as.numeric
        )
      )

  summary <- data$summary
  summary_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                          control_h.mu_hat = data_gen_params_h$control_h$mu,
                          control_h.s = data_gen_params_h$control_h$sigma / sqrt(data_gen_params_h$control_h$n))
  summary <- cbind(summary, summary_h)

  # Define borrowing scenarios
  prior_params_list <- list(
    # borrowing on the control arm
    Yes = list(treatment.theta0 = 0, treatment.s0 = 100, treatment.w = 0,
               control.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
               control.s0 = as.numeric(sqrt(summary_h["control_h.n"] * (summary_h["control_h.s"])^2)),
               control.w = NULL, control.delta = 0.15, control.ess_h = n_list[["control"]]),
    # no borrowing
    No = list(treatment.theta0 = 0, treatment.s0 = 100, treatment.w = 0,
              control.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
              control.s0 = as.numeric(sqrt(summary_h["control_h.n"] * (summary_h["control_h.s"])^2)),
              control.w = 0, control.delta = 0.15, control.ess_h = n_list[["control"]])
  )

  cat("Start Bayesian analysis.\n")
  for (k in seq_along(prior_params_list)) {
    start_time_k <- Sys.time()

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
                                      data_summary = summary,
                                      prior_params = prior_params,
                                      arm_names = list(treatment = "treatment",
                                                       control = "control",
                                                       control_h = "control_h"),
                                      lrv = -0.1, tv = -0.2,
                                      fgr = 0.2, fsr = 0.1,
                                      posterior_infer = T,
                                      Lalonde_decision = T)


    post$metrics_post_dist <- calc_post_dist_metrics(endpoint = "continuous",
                                                     true_value = data$true_value$compare_true,
                                                     post$post_est_ci)
    post$settings <- settings

    bayes_results[[length(bayes_results)+1]] <- post

    end_time_k <- Sys.time()
    cat("Time for Bayesian analysis with borrowing", borrow_name, "=", round(difftime(end_time_k, start_time_k, units = "secs"), 2), "seconds\n\n")
  }

  end_time_i <- Sys.time()
  cat("Total time for data_gen_params set", i, "=", round(difftime(end_time_i, start_time_i, units = "secs"), 2), "seconds\n\n")
}

saveRDS(bayes_results, file = "simulations/sim_DpR_constrained_bayes_results.rds")
bayes_results <- readRDS(file = "simulations/sim_DpR_constrained_bayes_results.rds")

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

saveRDS(bayes_results_all, "simulations/sim_DpR_constrained_bayes_results_df.rds")


# Analysis -------------------------------------------------------------

bayes_results_all <- readRDS("simulations/sim_DpR_constrained_bayes_results_df.rds")

## Check posterior parameters
post_params_all <- bayes_results_all$post_params_all
post_params_all <- post_params_all  %>%
  mutate(across(
    .cols = -c(borrowing, ends_with(".name")),
    .fns = as.numeric
  ))

post_params_all$control.mean.diff <- post_params_all$control.mu - post_params_all$control_h.mu
post_params_all$control.mean.diff <- factor(post_params_all$control.mean.diff,
                                            levels = c("-0.2", "-0.15", "-0.1","-0.05", "0", "0.05", "0.1", "0.15", "0.2"))

post_params_all$borrowing <- factor(post_params_all$borrowing,
                                    levels = c("No", "Yes"),
                                    labels=c("Borrowing: No", "Borrowing: Yes"))


post_params_all$control.n <- factor(post_params_all$control.n,
                                    levels = c("20", "30", "40", "50", "60"),
                                    labels=c(expression(paste(n, " = 20")),
                                             expression(paste(n, " = 30")),
                                             expression(paste(n, " = 40")),
                                             expression(paste(n, " = 50")),
                                             expression(paste(n, " = 60"))))

p1 <- post_params_all %>%
  filter(borrowing == "Borrowing: Yes") %>% # Use filter() for subsetting
  ggplot(aes(x = control.mean.diff, y = control.w_prior)) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  annotate("rect", xmin = 2, xmax = 8,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "blue") +
  facet_grid(~control.n,
             labeller = labeller(control.n = label_parsed)) +
  labs(x = expression(paste(theta[c], " - ", theta[ch])),
       y = "Prior Weight",
       title = "Prior Weight of Informative Component of SAM Prior") +
  theme_bw() +
  theme(strip.text.x = element_text(size=10))

p2 <- ggplot(post_params_all %>% subset(borrowing =="Borrowing: Yes" ),
             aes(x = control.mean.diff, y = control.w_post)) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  annotate("rect", xmin = 2, xmax = 8,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "blue") +
  facet_grid(~control.n,
             labeller = labeller(control.n = label_parsed)) +
  labs(x = expression(paste(theta[c], " - ", theta[ch])),
       y = "Posterior Weight",
       title = "Posterior Weight of Informative Component of SAM Prior (CSD = 0.15)") +
  theme_bw() +
  theme(strip.text.x = element_text(size=10))

p1 / p2

ggsave(filename = "simulations/sim_DpR_constrained_SAM_weight.jpg", p1 / p2,
       width = 13, height = 5)

###
metrics_post_dist_all <- bayes_results_all$metrics_post_dist_all
metrics_post_dist_all <- metrics_post_dist_all %>%
  mutate(across(
    .cols = -c(borrowing, ends_with(".name")),
    .fns = as.numeric
  ))

metrics_post_dist_all <- metrics_post_dist_all[, c("i", "true_value.compare_true",
                                                   "control.n", "borrowing",
                                                   "treatment.mu", "control.mu",
                                                   "control_h.mu",
                                                   "delta", "bias_avg_mean_diff",
                                                   "sd_avg", "sd_empirical", "cp")]

metrics_post_dist_all$control.mean.diff <- metrics_post_dist_all$control.mu - metrics_post_dist_all$control_h.mu

metrics_post_dist_all$true_value.compare_true <- round(metrics_post_dist_all$true_value.compare_true,2)
metrics_post_dist_all$true_value.compare_true = ifelse(metrics_post_dist_all$true_value.compare_true == -0.20, "-0.2 (TV)",
                                        ifelse(metrics_post_dist_all$true_value.compare_true == -0.10, "-0.1 (LRV)",
                                               "-0.15"))

metrics_post_dist_all$borrowing <- factor(metrics_post_dist_all$borrowing,
                                          levels = c("No", "Yes"),
                                          labels=c("Borrowing: No", "Borrowing: Yes"))

metrics_df <- metrics_post_dist_all %>%
  select(control.mean.diff,
         control.n,
         true_value.compare_true,
         borrowing,
         bias_avg_mean_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true, control.n)

metrics_df$control.mean.diff <- factor(metrics_df$control.mean.diff,
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

metrics_df$true_value.compare_true <- factor(metrics_df$true_value.compare_true,
                                             levels = c("-0.1 (LRV)", "-0.15", "-0.2 (TV)"),
                                                        labels=c(expression(paste(Delta, " = -0.1")),
                                                                 expression(paste(Delta, " = -0.15")),
                                                                 expression(paste(Delta, " = -0.2"))))

metrics_df$control.n <- factor(metrics_df$control.n)

p1 <- ggplot(metrics_df,
             aes(x = control.n, y = bias_avg_mean_diff, color = borrowing)) +
  geom_point(position = position_dodge(width = 0.1), size = 2) +
  facet_grid(true_value.compare_true~control.mean.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.mean.diff = label_parsed)) +
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
  facet_grid(true_value.compare_true~control.mean.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.mean.diff = label_parsed)) +
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
  ggh4x::facet_grid2(true_value.compare_true~control.mean.diff,
                     labeller = labeller(true_value.compare_true = label_parsed,
                                         control.mean.diff = label_parsed)) +
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
ggsave(filename = "simulations/sim_DpR_constrained_dist_check_metrics_1.jpg", p, width = 15, height = 12)

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
ggsave(filename = "simulations/sim_DpR_constrained_dist_check_metrics_2.jpg", p, width = 15.5, height = 10)


metrics_df <- metrics_post_dist_all %>%
  select(control.mean.diff, control.n,true_value.compare_true, borrowing,
         bias_avg_mean_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = borrowing, values_from = c(bias_avg_mean_diff, sd_avg, sd_empirical, cp)) %>%
  arrange(control.n) %>%
  arrange(control.mean.diff)

metrics_df <- metrics_df[, c("control.mean.diff", "control.n", "true_value.compare_true",
                             "bias_avg_mean_diff_Borrowing: No", "sd_avg_Borrowing: No",
                             "sd_empirical_Borrowing: No", "cp_Borrowing: No",
                             "bias_avg_mean_diff_Borrowing: Yes", "sd_avg_Borrowing: Yes",
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




################################
oc_all <- bayes_results_all$oc_all
oc_all <- oc_all[, c("i", "true_value.compare_true",
                     "control.n", "borrowing",
                     "treatment.mu", "control.mu", "control_h.mu",
                     "decision_pr", "proportion_pr", "proportion_ci")] %>%
  mutate(across(
    .cols = -c(ends_with(".name"), "borrowing", "decision_pr"),
    .fns = as.numeric
  ))

oc_all$decision_pr <- factor(oc_all$decision_pr, levels = c("no-go", "consider", "go"))

oc_all$true_value.compare_true = ifelse(oc_all$true_value.compare_true == -0.2, "-0.2 (TV)",
                                        ifelse(oc_all$true_value.compare_true == -0.1, "-0.1 (LRV)",
                                               "-0.15"))
oc_all$true_value.compare_true <- factor(oc_all$true_value.compare_true,
                                         levels = c("-0.1 (LRV)", "-0.15", "-0.2 (TV)"))


oc_all$control.mean.diff <- oc_all$control.mu - oc_all$control_h.mu
oc_all$borrowing <- factor(oc_all$borrowing,
                           levels = c("No", "Yes"),
                           labels=c("Borrowing: No", "Borrowing: Yes"))


metrics <- oc_all[oc_all$true_value.compare_true %in% c("-0.2 (TV)", "-0.1 (LRV)"), ]
metrics <- oc_all[(oc_all$true_value.compare_true == "-0.1 (LRV)" & oc_all$decision_pr == "go") | (oc_all$true_value.compare_true == "-0.2 (TV)" & oc_all$decision_pr == "no-go"),]
false_go_risk <- metrics[metrics$true_value.compare_true == "-0.1 (LRV)", c("control.n", "borrowing", "proportion_pr", "control.mean.diff")]
colnames(false_go_risk)[3] <- "Type I Error (FGR)"
false_stop_risk <- metrics[metrics$true_value.compare_true == "-0.2 (TV)", c("control.n", "borrowing", "proportion_pr", "control.mean.diff")]
colnames(false_stop_risk)[3] <- "Type II Error (FSR)"
metrics_df <- merge(false_go_risk, false_stop_risk, by = c("control.n", "borrowing", "control.mean.diff"))

metrics_long <- reshape2::melt(metrics_df, id.vars = c("control.n", "borrowing", "control.mean.diff"),
                               measure.vars = c("Type I Error (FGR)", "Type II Error (FSR)"),
                               variable.name = "risk_type", value.name = "risk_value")
metrics_long$control.n <- factor(metrics_long$control.n)

p_risk <- ggplot(metrics_long,
                 aes(x = control.mean.diff, y = risk_value, color = risk_type)) +
  geom_point(size = 2) +
  annotate("rect", xmin = -0.15, xmax = 0.15,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "blue") +
  facet_grid(borrowing ~ control.n,
             labeller = labeller(
               control.n = function(x) paste("n = ", x))) +
  labs(y = "Risk Value",
       color = "") +
  xlab(~ paste(theta[c], " - ", theta[ch])) +
  geom_hline(aes(yintercept = 0.2, color = "Type I Error (FGR)"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "Type II Error (FSR)"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
  )
ggsave("simulations/sim_DpR_constrained_conflict_vs_risk.jpg", p_risk, width = 13, height = 5)


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
                     "control.mean.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]
oc_new$control.mean.diff <- factor(oc_new$control.mean.diff)
myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)
oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true,
                                         levels = c("-0.1 (LRV)", "-0.15", "-0.2 (TV)"))


p_list <- list()
for (i in levels(oc_new$true_value.compare_true)) {
  oc_tmp <- oc_new[oc_new$true_value.compare_true == i, ]
  oc_tmp$control.mean.diff <- factor(oc_tmp$control.mean.diff,
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
    facet_grid(control.mean.diff~borrowing,
               labeller = labeller(control.mean.diff = label_parsed)) +
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

ggsave("simulations/sim_DpR_constrained_zone_size.jpg", p_oc, width = 14, height = 14)

oc2 <- oc_new %>%
  select(control.mean.diff, control.n, true_value.compare_true, borrowing, decision_pr, proportion_pr) %>%
  mutate(proportion_pr = scales::percent(proportion_pr, 0.01)) %>%
  arrange(decision_pr) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = c(decision_pr, borrowing), values_from = proportion_pr) %>%
  arrange(control.n) %>%
  arrange(control.mean.diff)

# oc2$`no-go_Borrowing: No` = cell_spec(oc2$`no-go_Borrowing: No`,
#                                       background = ifelse(oc2$true_value.compare_true == "-0.2 (TV)", "steelblue2", "white"))
#
# oc2$`no-go_Borrowing: Yes` = cell_spec(oc2$`no-go_Borrowing: Yes`,
#                                        background = ifelse(oc2$true_value.compare_true == "-0.2 (TV)", "steelblue2", "white"))
#
# oc2$`go_Borrowing: No` = cell_spec(oc2$`go_Borrowing: No`,
#                                    background = ifelse(oc2$true_value.compare_true == "-0.1 (LRV)", "pink", "white"))
#
# oc2$`go_Borrowing: Yes` = cell_spec(oc2$`go_Borrowing: Yes`,
#                                     background = ifelse(oc2$true_value.compare_true == "-0.1 (LRV)", "pink", "white"))

kbl(oc2[,-1], escape = T, row.names = F, digits = 3,
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












