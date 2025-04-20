library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)


# Simulations -------------------------------------------------------------

param_grid_h <- expand.grid(
  ctrl_h_n = 200, # Sample size for the historical control group.
  ctrl_h_prob = 0.3, # Probability of zero g-scores in the historical control group.
  stringsAsFactors = FALSE
)

param_grid1 <- expand.grid(
  trt_n = c(20, 40, 60), # Sample size for the treatment group.
  ctrl_prob = seq(0.1, 0.5, 0.1), # Probability of zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_prob = 0.1 + ctrl_prob)

param_grid2 <- expand.grid(
  trt_n = c(20, 40, 60), # Sample size for the treatment group.
  ctrl_prob = seq(0.1, 0.5, 0.1), # Probability of zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_prob = 0.15 + ctrl_prob)

param_grid3 <- expand.grid(
  trt_n = c(20, 40, 60), # Sample size for the treatment group.
  ctrl_prob = seq(0.1, 0.5, 0.1), # Probability of zero g-scores in the control group.
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n,
         trt_prob = 0.2 + ctrl_prob)

param_grid <- bind_rows(param_grid1, param_grid2, param_grid3)

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                               create_data_gen_params, endpoint = "OR")

data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "OR")

prior_params_list <- list(
  Yes = list(treatment.a = 0.01, treatment.b = 0.01, control.delta = 0.1, control.a = 0.01, control.b = 0.01, treatment.w = 0, control.w = NULL), # borrowing on the control arm
  No = list(treatment.a = 0.01, treatment.b = 0.01, control.delta = 0.1, control.a = 0.01, control.b = 0.01, treatment.w = 0, control.w = 0) # no borrowing
)

nsim = 10000
bayes_results <- list()
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
  freq <- data$freq
  freq_h <- data.frame(control_h.n = data_gen_params_h$control_h$n,
                       control_h.count = data_gen_params_h$control_h$n*data_gen_params_h$control_h$prob)
  freq <- cbind(freq, freq_h)

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

    post <- bayesian_lalonde_decision(endpoint = "binary",
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

# for (i in seq_along(bayes_results)) {
#   res <- bayes_results[[i]]$post_est_ci %>%
#     na.omit()
#
#   metrics_post_dist <- calc_post_dist_metrics(endpoint = "OR",
#                                               true_value = bayes_results[[i]]$settings$true_value.compare_true,
#                                               res)
#   bayes_results[[i]]$metrics_post_dist_no_na <- metrics_post_dist
# }


# saveRDS(bayes_results, file = "sim_OR_bayes_results_10_10.rds")
# bayes_results <- readRDS("sim_OR_bayes_results_10_10.rds")

post_params_all <- data.frame()
post_est_ci_all <- data.frame()
post_inference_all <- data.frame()
metrics_post_dist_all <- data.frame()
metrics_post_dist_no_na_all <- data.frame()
oc_all <- data.frame()
for (i in seq_along(bayes_results)) {
  res <- bayes_results[[i]]
  setting <- res$settings
  post_params <- cbind(setting, res$post_params)
  post_est_ci <- cbind(setting, res$post_est_ci)
  post_inference <- cbind(setting, res$post_inference)
  metrics_post_dist <- cbind(setting, res$metrics_post_dist)

  metrics_post_dist_no_na <- calc_post_dist_metrics(endpoint = "OR",
                                                    true_value = setting$true_value.compare_true,
                                                    post_est_ci = res$post_est_ci, remove.na = T)
  metrics_post_dist_no_na <- cbind(setting, metrics_post_dist_no_na)

  post_inference$decision_pr <- ifelse(is.na(post_inference$decision_pr), post_inference$decision_ci, post_inference$decision_pr)
  post_inference$decision_ci <- ifelse(is.na(post_inference$decision_ci), post_inference$decision_pr, post_inference$decision_ci)

  oc <- cbind(setting, obtain_oc(post_inference))

  post_params_all <- bind_rows(post_params_all, cbind(i, post_params))
  post_est_ci_all <- bind_rows(post_est_ci_all, cbind(i, post_est_ci))
  post_inference_all <- bind_rows(post_inference_all, cbind(i, post_inference))
  metrics_post_dist_all <- bind_rows(metrics_post_dist_all, cbind(i, metrics_post_dist))
  metrics_post_dist_no_na_all <- bind_rows(metrics_post_dist_no_na_all, cbind(i,metrics_post_dist_no_na))
  oc_all <- bind_rows(oc_all, cbind(i, oc))
}

bayes_results_all <- list(
  metrics_post_dist_all = metrics_post_dist_all,
  metrics_post_dist_no_na_all = metrics_post_dist_no_na_all,
  post_params_all = post_params_all,
  post_est_ci_all = post_est_ci_all,
  post_inference_all = post_inference_all,
  oc_all = oc_all
)

saveRDS(bayes_results_all, "sim_OR_bayes_results_df_10_10.rds")


# Analysis ----------------------------------------------------------------

bayes_results_all <- readRDS("sim_OR_bayes_results_df_10_10.rds")

## Check posterior parameters
post_params_all <- bayes_results_all$post_params_all
post_params_all <- post_params_all  %>%
  mutate(across(
    .cols = -c(borrowing, ends_with(".name")),
    .fns = as.numeric
  ))

post_params_all$control.prob.diff <- post_params_all$control.prob - post_params_all$control_h.prob
post_params_all$true_value.compare_true = ifelse(post_params_all$true_value.compare_true == 0.1, "0.1 (LRV)",
                                                       ifelse(post_params_all$true_value.compare_true == 0.2, "0.2 (TV)",
                                                              post_params_all$true_value.compare_true))
post_params_all$borrowing <- factor(post_params_all$borrowing,
                                    levels = c("No", "Yes"),
                                    labels=c("Borrowing: No", "Borrowing: Yes"))

ggplot(post_params_all,
       aes(x = control.prob.diff, y = control.w_post)) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, outlier.alpha = 0.5) +
  facet_grid(true_value.compare_true~control.n, labeller = label_both) +
  labs(x = "Ratio of the mean of g-score > 0 comparing current to historical control",
       y = "Posterior weight of borrowing",
       color = "Prob of g-score = 0 in current control",
       title = "Mean of g-score > 0 in current control = 0.01, delta for SAM = log(0.85)") +
  theme_bw() +
  theme(legend.position = "bottom", strip.text.x = element_text(size=10))


post_params_all$zeroweight <- ifelse(post_params_all$control.w_post == 0, "zero", "nonzero")
table(post_params_all$control.prob.diff, post_params_all$zeroweight)

#
# post_param <- combined_post_params[combined_post_params$control.prob == "0.1", ]
# post_params <- post_param[1, c("control.w1_post", "control.mu1_post", "control.sigma1_post", "control.mu2_post", "control.sigma2_post")]
# names(post_params)[1] <- "control.w_post"
#
# df <- readRDS("sim_rslt_gscore_decisions_7_8.rds")
# obtain_ESS <- function(post_params, sigma){
#
#   post <- convert_RBesT_mix(post = post_params, endpoint = "g-score")
#   ess_post <- ess(post, sigma = sigma)
# }
#
#
#
metrics_post_dist_all <- bayes_results_all$metrics_post_dist_no_na_all
metrics_post_dist_no_na_all <- bayes_results_all$metrics_post_dist_no_na_all
metrics_post_dist_all <- metrics_post_dist_all  %>%
  mutate(across(
    .cols = -c(borrowing, ends_with(".name")),
    .fns = as.numeric
  ))
metrics_post_dist_all <- metrics_post_dist_all[, c("i", "true_value.compare_true",
                                                   "control.n", "borrowing",
                                                   "treatment.prob", "control.prob",
                                                   "control_h.prob",
                                                   "delta", "bias_avg_rate_diff",
                                                   "sd_avg", "sd_empirical", "cp")]
metrics_post_dist_all$control.prob.diff <- metrics_post_dist_all$control.prob - metrics_post_dist_all$control_h.prob
metrics_post_dist_all$true_value.compare_true = ifelse(metrics_post_dist_all$true_value.compare_true == 0.1, "0.1 (LRV)",
                                                       ifelse(metrics_post_dist_all$true_value.compare_true == 0.2, "0.2 (TV)",
                                                              metrics_post_dist_all$true_value.compare_true))
metrics_post_dist_all$borrowing <- factor(metrics_post_dist_all$borrowing,
                                          levels = c("No", "Yes"),
                                          labels=c("Borrowing: No", "Borrowing: Yes"))
metrics_df <- metrics_post_dist_all %>%
  select(control.prob.diff, control.n,true_value.compare_true, borrowing,
         bias_avg_rate_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true, control.n)

metrics_df$control.prob.diff <- factor(metrics_df$control.prob.diff,
                                       levels = c("-0.2", "-0.1", "0", "0.1", "0.2"),
                                       labels=c(expression(paste(theta[c], " - ", theta[ch], " = -0.2")),
                                                expression(paste(theta[c], " - ", theta[ch], " = -0.1")),
                                                expression(paste(theta[c], " - ", theta[ch], " = 0")),
                                                expression(paste(theta[c], " - ", theta[ch], " = 0.1")),
                                                expression(paste(theta[c], " - ", theta[ch], " = 0.2"))))
metrics_df$true_value.compare_true <- factor(metrics_df$true_value.compare_true,
                                       levels = c("0.1 (LRV)", "0.15", "0.2 (TV)"),
                                       labels=c(expression(paste(Delta, " = 0.1")),
                                                expression(paste(Delta, " = 0.15")),
                                                expression(paste(Delta, " = 0.2"))))
metrics_df$control.n <- factor(metrics_df$control.n)

p1 <- ggplot(metrics_df,
             aes(x = control.n, y = bias_avg_rate_diff, color = borrowing)) +
  geom_point(position = position_dodge(width = 0.1), size = 2) +
  facet_grid(true_value.compare_true~control.prob.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.prob.diff = label_parsed)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "BIAS", color= "") +
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
  facet_grid(true_value.compare_true~control.prob.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.prob.diff = label_parsed)) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
  labs(x = "Sample size per arm", y = "CP", color = "") +
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
  ggh4x::facet_grid2(true_value.compare_true~control.prob.diff,
                     labeller = labeller(true_value.compare_true = label_parsed,
                                         control.prob.diff = label_parsed)) +
  labs(x = "Sample size per arm", y = "SD",
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
ggsave(filename = "sim_OR_dist_check_metrics_10_10.jpg", p, width = 9.5, height = 11.5)

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
ggsave(filename = "sim_OR_dist_check_metrics.jpg", p, width = 8.5, height = 8.5)



metrics_df <- metrics_post_dist_all %>%
  select(control.prob.diff, control.n,true_value.compare_true, borrowing,
         bias_avg_rate_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = borrowing, values_from = c(bias_avg_rate_diff, sd_avg, sd_empirical, cp)) %>%
  arrange(control.n) %>%
  arrange(control.prob.diff)
metrics_df <- metrics_df[, c("control.prob.diff", "control.n", "true_value.compare_true",
                             "bias_avg_rate_diff_Borrowing: No", "sd_avg_Borrowing: No",
                             "sd_empirical_Borrowing: No", "cp_Borrowing: No",
                             "bias_avg_rate_diff_Borrowing: Yes", "sd_avg_Borrowing: Yes",
                             "sd_empirical_Borrowing: Yes", "cp_Borrowing: Yes")]

kbl(metrics_df[,-1], escape = F, row.names = F, digits = 3,
    format    = "latex",
    booktabs  = T,
    col.names = c("Sample size", "Delta",
                  "BIAS", "SE-avg","SD-emp", "CP",
                  "BIAS", "SE-avg","SD-emp", "CP")) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 2, "Borrowing: No" = 4, "Borrowing: Yes" = 4)) %>%
  pack_rows(index = c("Case 1" = 9,
                      "Case 2" = 9,
                      "Case 3" = 9,
                      "Case 4" = 9,
                      "Case 5" = 9))


oc_all <- bayes_results_all$oc_all
oc_all <- oc_all[, c("i", "true_value.compare_true", "control.n", "borrowing",
                     "treatment.prob", "control.prob", "control_h.prob",
                     "decision_pr", "proportion_pr", "proportion_ci")]
# oc_new <- data.frame()
# for (i in unique(oc_all$i)) {
#   oc_tmp <- oc_all[oc_all$i == i,]
#   # Check if there are any NA values in the 'decision_pr' column
#   if (sum(is.na(oc_tmp$decision_pr)) > 0) {
#     oc_na <- oc_tmp[is.na(oc_tmp$decision_pr), ]
#
#     if (sum(is.na(oc_na$proportion_pr)) < sum(is.na(oc_na$proportion_ci))) {
#       oc_tmp$proportion <- oc_tmp$proportion_pr
#     } else if (sum(is.na(oc_na$proportion_ci)) < sum(is.na(oc_na$proportion_pr))) {
#       oc_tmp$proportion <- oc_tmp$proportion_ci
#     } else {
#       if (oc_na$proportion_ci < oc_na$proportion_pr) {
#         oc_tmp$proportion <- oc_tmp$proportion_ci
#       } else {
#         oc_tmp$proportion <- oc_tmp$proportion_pr
#       }
#     }
#   }
#
#
#   }else{
#     oc_new <- cbind(oc_new, oc_tmp)
#   }
#
# }
# oc_all$proportion <- ifelse(!is.na(oc_all$decision_pr), oc_all$proportion_pr,
#                             ifelse(is.na(oc_all$proportion_ci), oc_all$proportion_pr, oc_all$proportion_ci))


# oc_all$decision_pr <- factor(oc_all$decision_pr, levels = c("no-go", "consider", "go", "NA"))
oc_all$decision_pr <- factor(oc_all$decision_pr, levels = c("no-go", "consider", "go"))
oc_all$true_value.compare_true = ifelse(oc_all$true_value.compare_true == 0.1, "0.1 (LRV)",
                                        ifelse(oc_all$true_value.compare_true == 0.2, "0.2 (TV)",
                                               oc_all$true_value.compare_true))

oc_all$control.prob <- as.numeric(oc_all$control.prob)
oc_all$control_h.prob <- as.numeric(oc_all$control_h.prob)
oc_all$control.prob.diff <- oc_all$control.prob - oc_all$control_h.prob
oc_all$borrowing <- factor(oc_all$borrowing,
                           levels = c("No", "Yes"),
                           labels=c("Borrowing: No", "Borrowing: Yes"))

metrics <- oc_all[oc_all$true_value.compare_true %in% c("0.2 (TV)", "0.1 (LRV)"),]
metrics <- oc_all[(oc_all$true_value.compare_true == "0.1 (LRV)" & oc_all$decision_pr == "go") | (oc_all$true_value.compare_true == "0.2 (TV)" & oc_all$decision_pr == "no-go"),]
false_go_risk <- metrics[metrics$true_value.compare_true == "0.1 (LRV)", c("control.n", "borrowing", "proportion_pr", "control.prob.diff")]
colnames(false_go_risk)[3] <- "Type I Error (FGR)"
false_stop_risk <- metrics[metrics$true_value.compare_true == "0.2 (TV)", c("control.n", "borrowing", "proportion_pr", "control.prob.diff")]
colnames(false_stop_risk)[3] <- "Type II Error (FSR)"
metrics_df <- merge(false_go_risk, false_stop_risk, by = c("control.n", "borrowing", "control.prob.diff"))

metrics_long <- reshape2::melt(metrics_df, id.vars = c("control.n", "borrowing", "control.prob.diff"),
                     measure.vars = c("Type I Error (FGR)", "Type II Error (FSR)"),
                     variable.name = "risk_type", value.name = "risk_value")

p_risk <- ggplot(metrics_long,
       aes(x = control.prob.diff, y = risk_value, color = risk_type)) +
  geom_point(size = 2) +
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
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,0),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
        )
ggsave("sim_OR_conflict_vs_risk.jpg", width = 7, height = 4)


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
                     "control.prob", "control.prob.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]
oc_new$control.prob.diff <- factor(oc_new$control.prob.diff)
myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)
oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
# lines <- data.frame(true_value.compare_true = factor(c("0.1 (LRV)", "0.2 (TV)"),
#                                                      levels = c("0.1 (LRV)","0.2 (TV)"),
#                                                      labels = c(expression(paste(Delta, " = 0.1 (LRV)")),
#                                                                 expression(paste(Delta, " = 0.2 (TV)")))),
#                     Z = c(0.9, 0.2))

# p_list <- list()
# for (i in levels(oc_new$control.prob.diff)) {
#   oc_tmp <- oc_new[oc_new$control.prob.diff == i, ]
#   oc_tmp$true_value.compare_true <- factor(oc_tmp$true_value.compare_true ,
#                                            levels = c("0.1 (LRV)", "0.15", "0.2 (TV)"),
#                                            labels = c(expression(paste(Delta, " = 0.1 (LRV)")),
#                                                       expression(paste(Delta, " = 0.15")),
#                                                       expression(paste(Delta, " = 0.2 (TV)"))))
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
#     ggtitle(bquote(theta[c]-theta[ch] == .(i))) +
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
#

p_list <- list()
for (i in levels(oc_new$true_value.compare_true)) {
  oc_tmp <- oc_new[oc_new$true_value.compare_true == i, ]
  oc_tmp$control.prob.diff <- factor(oc_tmp$control.prob.diff,
                                     levels = c("-0.2", "-0.1", "0", "0.1", "0.2"),
                                     labels=c(expression(paste(theta[c], " - ", theta[ch], " = -0.2")),
                                              expression(paste(theta[c], " - ", theta[ch], " = -0.1")),
                                              expression(paste(theta[c], " - ", theta[ch], " = 0")),
                                              expression(paste(theta[c], " - ", theta[ch], " = 0.1")),
                                              expression(paste(theta[c], " - ", theta[ch], " = 0.2"))))
  p_list[[i]] <- ggplot(oc_tmp,
                        aes(x = as.factor(control.n), y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.1)), vjust = 1.6,
              # fontface = "bold",
              size = 3.2) +
    scale_fill_manual(name = "Decision", values = myColors)+
    facet_grid(control.prob.diff~borrowing,
               labeller = labeller(control.prob.diff = label_parsed)) +
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
ggsave("sim_OR_zone_size.jpg",p_oc, width = 14, height = 8)


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



oc2 <- oc_all %>%
  select(control.prob.diff, control.n, true_value.compare_true, borrowing, decision_pr, proportion_pr) %>%
  mutate(proportion_pr = scales::percent(proportion_pr, 0.1)) %>%
  arrange(decision_pr) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = c(decision_pr, borrowing), values_from = proportion_pr) %>%
  arrange(control.n) %>%
  arrange(control.prob.diff)

# oc2$`no-go_Borrowing: No` = cell_spec(oc2$`no-go_Borrowing: No`,
#                                       background = ifelse(oc2$true_value.compare_true == "0.2 (TV)", "steelblue2", "white"))
#
# oc2$`no-go_Borrowing: Yes` = cell_spec(oc2$`no-go_Borrowing: Yes`,
#                                        background = ifelse(oc2$true_value.compare_true == "0.2 (TV)", "steelblue2", "white"))
#
# oc2$`go_Borrowing: No` = cell_spec(oc2$`go_Borrowing: No`,
#                                    background = ifelse(oc2$true_value.compare_true == "0.1 (LRV)", "pink", "white"))
#
# oc2$`go_Borrowing: Yes` = cell_spec(oc2$`go_Borrowing: Yes`,
#                                     background = ifelse(oc2$true_value.compare_true == "0.1 (LRV)", "pink", "white"))

kbl(oc2[,-1], escape = T, row.names = F, digits = 3,
    format    = "latex",
    booktabs  = T,
    col.names = c("n", "Delta", rep(c("Pr(No-Go)", "Pr(Consider)","Pr(Go)"), 2))) %>%
  kable_styling(bootstrap_options = c("condensed")) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 2, "Borrowing: No" = 3, "Borrowing: Yes" = 3)) %>%
  pack_rows(index = c("Case 1:" = 9,
                      "Case 2:" = 9,
                      "Case 3:" = 9,
                      "Case 4:" = 9,
                      "Case 5:" = 9))


