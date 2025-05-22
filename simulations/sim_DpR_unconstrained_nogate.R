rm(list = ls())
library(kableExtra)
library(LalondeBayesBorrow)
library(dplyr)
library(tidyverse)

# Simulation ------------------------------------------------------------

# Define historical control data
param_grid_h <- expand.grid(
  ctrl_h_n = 180,
  ctrl_h_mu = -0.2,
  ctrl_h_sigma = 0.4,
  stringsAsFactors = FALSE
)

# Define current trial configurations
param_grid1 <- expand.grid(
  trt_n = c(45),
  ctrl_mu = seq(-0.35, -0.05, 0.05),
  trt_sigma = 0.4,
  ctrl_sigma = 0.4,
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

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "continuous")

data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "continuous")

## Run simulations
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

  # Generate simulation dataset
  data <- data_gen_continuous(arm_names = arm_names,
                              nsim = nsim, n_list = n_list,
                              mu_list = mu_list, sigma_list = sigma_list,
                              seed = 123)

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

  prior_grid1 <- expand.grid(
    treatment.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    treatment.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    treatment.w = 0,
    control.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    control.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    control.w = c(NA, 0),
    control.delta_gate = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  ) %>%
    mutate(control.delta_SAM = control.delta_gate)

  prior_grid2 <- expand.grid(
    treatment.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    treatment.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    treatment.w = 0,
    control.theta0 = as.numeric(summary_h["control_h.mu_hat"]),
    control.s0 = as.numeric(data_gen_params_h$control_h$sigma),
    control.w = c(NA, 0),
    control.delta_SAM = c(0.1, 0.15),
    control.ess_h = c(45, 90, 180),
    stringsAsFactors = FALSE
  )

  prior_grid <- bind_rows(prior_grid1, prior_grid2)

  # Convert to list and replace NA with NULL
  prior_params_list <- apply(prior_grid, 1, function(row) {
    out <- as.list(row)
    if (is.na(out$control.w)) out$control.w <- NULL
    if (is.na(out$control.delta_gate)) out$control.delta_gate <- NULL
    out
  })

  cat("proportion of |conflict| < 0.1 =", sum(abs(summary$control.mu_hat-summary$control_h.mu_hat) < 0.1)/nsim, "\n")
  cat("proportion of |conflict| < 0.15 =", sum(abs(summary$control.mu_hat-summary$control_h.mu_hat) < 0.15)/nsim, "\n")

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

    post <- bayesian_lalonde_decision(endpoint = "continuous",
                                      data_summary = summary,
                                      prior_params = prior_params,
                                      arm_names = c(treatment = "treatment",
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
    post$data_summary <- cbind(i, summary)

    bayes_results[[length(bayes_results)+1]] <- post

    end_time_k <- Sys.time()
    cat("Time for Bayesian analysis with prior parameter list", k, "=", round(difftime(end_time_k, start_time_k, units = "secs"), 2), "seconds\n\n")
  }

  end_time_i <- Sys.time()
  cat("Total time for data_gen_params set", i, "=", round(difftime(end_time_i, start_time_i, units = "secs"), 2), "seconds\n\n")
}

saveRDS(bayes_results, file = "simulations/sim_DpR_unconstrained_bayes_results_5_19.rds")

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

saveRDS(bayes_results_all, "simulations/sim_DpR_unconstrained_bayes_results_df_5_19.rds")


# Analysis -------------------------------------------------------------

bayes_results_all <- readRDS("simulations/sim_DpR_unconstrained_bayes_results_df_5_19.rds")

# + Check posterior weights -----------------------------------------------

post_params_all <- bayes_results_all$post_params_all %>%
  mutate(across(
    .cols = -c(ends_with(".name")),
    .fns = as.numeric
  )) %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.mu.diff = control.mu - control_h.mu) %>%
  filter(borrowing == "Borrowing: Yes") %>%
  select("nsim", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.mu.diff",
         "control.w_prior", "control.w_post") %>%
  mutate(control.mu.diff = factor(control.mu.diff,
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

p1 <- ggplot(post_params_all, aes(x = control.mu.diff, y = control.w_prior)) +
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

p2 <- ggplot(post_params_all, aes(x = control.mu.diff, y = control.w_post)) +
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
  theme(strip.text.x = element_text(size=10))

library(patchwork)
p1 / p2
ggsave(filename = "simulations/sim_DpR_unconstrained_SAM_weight_5_19.jpg", p1 / p2,
       width = 13, height = 8)


# + Generate summary table for paper --------------------------------------

post_inference_all <-  bayes_results_all$post_inference_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.mu.diff = control.mu - control_h.mu) %>%
  select("nsim", "borrowing", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.mu",
         "est2_lalonde") %>%
  pivot_wider(names_from = borrowing, values_from = est2_lalonde) %>%
  mutate(pmd = `Borrowing: Yes` - `Borrowing: No`) %>%
  select(-c(`Borrowing: Yes`, `Borrowing: No`, true_value.compare_true))


# Summary statistics of PMD
metrics_tab <- post_inference_all %>%
  group_by(gate, control.delta_SAM, control.ess_h, control.mu) %>%
  summarize(mean_pmd = mean(pmd),
            sd_pmd = sd(pmd),
            .groups = "drop")

# risks
risk_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.mu.diff = control.mu - control_h.mu) %>%
  filter(true_value.compare_true %in% c(-0.1, -0.2),
         decision_pr %in% c("go", "no-go")) %>%
  mutate(
    CorrectLabel = case_when(
      true_value.compare_true == -0.2 & decision_pr == "go" ~ "CGR",
      true_value.compare_true == -0.1 & decision_pr == "no-go" ~ "CSR",
      TRUE ~ NA_character_
    ),
    RiskLabel = case_when(
      true_value.compare_true == -0.1 & decision_pr == "go" ~ "FGR",
      true_value.compare_true == -0.2 & decision_pr == "no-go" ~ "FSR",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(cols = c(CorrectLabel, RiskLabel), names_to = "type", values_to = "label") %>%
  filter(!is.na(label)) %>%
  group_by(gate, borrowing, control.delta_SAM, control.ess_h, control.mu, label) %>%
  summarize(prop = mean(proportion_pr), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = prop)

metrics_tab <- merge(metrics_tab, risk_df,
                     by = c("gate", "control.delta_SAM", "control.ess_h","control.mu"),
                     all = TRUE)

metrics_tab_1 <- metrics_tab %>%
  mutate(mean_pmd = ifelse(borrowing == "Borrowing: No", NA, mean_pmd),
         sd_pmd = ifelse(borrowing == "Borrowing: No", NA, sd_pmd)) %>%
  pivot_longer(c(FGR, FSR, CGR, CSR, mean_pmd, sd_pmd)) %>%
  pivot_wider(names_from = c(borrowing, control.delta_SAM, control.ess_h),
              values_from = value, names_sort = T) %>%
  mutate(name = factor(name, levels = c("FGR", "FSR", "CGR", "CSR", "mean_pmd", "sd_pmd"))) %>%
  select(c(1,3,2,4,10:15)) %>%
  arrange(gate, name, control.mu) %>%
  select(c(1:3, 5:10, 4))

# kbl(metrics_tab_1,
#     escape = F,
#     # format = "latex",
#     row.names = FALSE,
#     digits = 4,
#     col.names = c("Metric", "control.mu",
#                   rep(c("45", "90", "180"), 2),
#                   "Borrowing")) %>%
#   kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
#   kable_classic(full_width = FALSE) %>%
#   add_header_above(c(" " = 1, " " = 1,
#                      "Borrowing with delta = 0.1" = 3,
#                      "Borrowing with delta = 0.15" = 3,
#                      "No" = 1)) %>%
#   collapse_rows(columns = 1, valign = "middle")

metrics_tab_1 <- metrics_tab_1 %>%
  pivot_wider(names_from = gate, values_from = 4:10, names_vary = "slowest")

kbl(metrics_tab_1[, -1],
    escape = FALSE,
    row.names = FALSE,
    digits = 4,
    col.names = c(
      "control.mu",
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


# Visualize risks
risk_df <- risk_df %>%
  pivot_longer(c(FSR, FGR)) %>%
  mutate(control.mu = factor(control.mu,
                                 levels = c("-0.35", "-0.3", "-0.25", "-0.2", "-0.15", "-0.1", "-0.05")),
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
                 aes(x = control.mu, y = value, color = name, shape = gate)) +
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
  xlab(~ paste(theta[c])) +
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

ggsave( "simulations/sim_DpR_unconstrained_conflict_vs_risk_5_19.jpg", p_risk, width = 10, height = 5)




# + Generate plots for OC -------------------------------------------------

oc_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.mu.diff = control.mu - control_h.mu) %>%
  select("i", "gate", "borrowing",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.mu.diff",
         "decision_pr", "proportion_pr", "proportion_ci") %>%
  mutate(decision_pr = factor(decision_pr, levels = c("no-go", "consider", "go"),
                              labels = c("No-Go", "Consider", "Go")),
         true_value.compare_true = ifelse(true_value.compare_true == -0.2, "-0.2 (TV)", "-0.1 (LRV)"))


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
                     "control.mu.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]

oc_new$control.mu.diff <- factor(oc_new$control.mu.diff)
myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_new$decision_pr)

oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
# oc_new$control.mu.diff <- factor(oc_new$control.mu.diff,
#                                  levels = c("-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15"),
#                                  labels=c(expression(paste(theta[c], " - ", theta[ch], " = -0.15")),
#                                           expression(paste(theta[c], " - ", theta[ch], " = -0.1")),
#                                           expression(paste(theta[c], " - ", theta[ch], " = -0.05")),
#                                           expression(paste(theta[c], " - ", theta[ch], " = 0")),
#                                           expression(paste(theta[c], " - ", theta[ch], " = 0.05")),
#                                           expression(paste(theta[c], " - ", theta[ch], " = 0.1")),
#                                           expression(paste(theta[c], " - ", theta[ch], " = 0.15"))))

oc_new$control.mu.diff <- factor(oc_new$control.mu.diff,
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
                        aes(x = control.mu.diff, y = proportion_pr, fill = decision_pr)) +
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
ggsave("simulations/sim_DpR_unconstrained_zone_size_5_19.jpg", p_oc, width = 20, height = 15)

#
# oc2 <- oc_new %>%
#   select(control.mu.diff, control.n, true_value.compare_true, borrowing, decision_pr, proportion_pr) %>%
#   mutate(proportion_pr = scales::percent(proportion_pr, 0.01)) %>%
#   arrange(decision_pr) %>%
#   arrange(true_value.compare_true) %>%
#   arrange(borrowing) %>%
#   pivot_wider(names_from = c(decision_pr, borrowing), values_from = proportion_pr) %>%
#   arrange(control.n) %>%
#   arrange(control.mu.diff)
#
# # oc2$`no-go_Borrowing: No` = cell_spec(oc2$`no-go_Borrowing: No`,
# #                                       background = ifelse(oc2$true_value.compare_true == "-0.2 (TV)", "steelblue2", "white"))
# #
# # oc2$`no-go_Borrowing: Yes` = cell_spec(oc2$`no-go_Borrowing: Yes`,
# #                                        background = ifelse(oc2$true_value.compare_true == "-0.2 (TV)", "steelblue2", "white"))
# #
# # oc2$`go_Borrowing: No` = cell_spec(oc2$`go_Borrowing: No`,
# #                                    background = ifelse(oc2$true_value.compare_true == "-0.1 (LRV)", "pink", "white"))
# #
# # oc2$`go_Borrowing: Yes` = cell_spec(oc2$`go_Borrowing: Yes`,
# #                                     background = ifelse(oc2$true_value.compare_true == "-0.1 (LRV)", "pink", "white"))
#
# kbl(oc2[,-1], escape = T, row.names = F, digits = 3,
#     format    = "latex",
#     booktabs  = T,
#     col.names = c("n", "Delta", rep(c("Pr(No-Go)", "Pr(Consider)","Pr(Go)"), 2))) %>%
#   kable_styling(bootstrap_options = c("condensed")) %>%
#   kable_classic(full_width = FALSE) %>%
#   add_header_above(c(" " = 2, "Borrowing: No" = 3, "Borrowing: Yes" = 3)) %>%
#   pack_rows(index = c("Case 1" = 15,
#                       "Case 2" = 15,
#                       "Case 3" = 15,
#                       "Case 4" = 15,
#                       "Case 5" = 15,
#                       "Case 6" = 15,
#                       "Case 7" = 15,
#                       "Case 8" = 15,
#                       "Case 9" = 15))
#



######

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

metrics_post_dist_all$control.mu.diff <- metrics_post_dist_all$control.mu - metrics_post_dist_all$control_h.mu

metrics_post_dist_all$true_value.compare_true <- round(metrics_post_dist_all$true_value.compare_true,2)
metrics_post_dist_all$true_value.compare_true = ifelse(metrics_post_dist_all$true_value.compare_true == -0.20, "-0.2 (TV)",
                                                       ifelse(metrics_post_dist_all$true_value.compare_true == -0.10, "-0.1 (LRV)",
                                                              "-0.15"))

metrics_post_dist_all$borrowing <- factor(metrics_post_dist_all$borrowing,
                                          levels = c("No", "Yes"),
                                          labels=c("Borrowing: No", "Borrowing: Yes"))

metrics_df <- metrics_post_dist_all %>%
  select(control.mu.diff,
         control.n,
         true_value.compare_true,
         borrowing,
         bias_avg_mean_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true, control.n)

metrics_df$control.mu.diff <- factor(metrics_df$control.mu.diff,
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
  facet_grid(true_value.compare_true~control.mu.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.mu.diff = label_parsed)) +
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
  facet_grid(true_value.compare_true~control.mu.diff,
             labeller = labeller(true_value.compare_true = label_parsed,
                                 control.mu.diff = label_parsed)) +
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
  ggh4x::facet_grid2(true_value.compare_true~control.mu.diff,
                     labeller = labeller(true_value.compare_true = label_parsed,
                                         control.mu.diff = label_parsed)) +
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
ggsave(filename = "simulations/sim_DpR_unconstrained_dist_check_metrics_1.jpg", p, width = 15, height = 12)

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
ggsave(filename = "simulations/sim_DpR_unconstrained_dist_check_metrics_2.jpg", p, width = 15.5, height = 10)


metrics_df <- metrics_post_dist_all %>%
  select(control.mu.diff, control.n,true_value.compare_true, borrowing,
         bias_avg_mean_diff, sd_avg, sd_empirical, cp) %>%
  arrange(true_value.compare_true) %>%
  arrange(borrowing) %>%
  pivot_wider(names_from = borrowing, values_from = c(bias_avg_mean_diff, sd_avg, sd_empirical, cp)) %>%
  arrange(control.n) %>%
  arrange(control.mu.diff)

metrics_df <- metrics_df[, c("control.mu.diff", "control.n", "true_value.compare_true",
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




