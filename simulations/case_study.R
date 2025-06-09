rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(scales)

# 1. Sample Size Determination (No Conflict) ------------------------------------

ctrl_n_values <- seq(20, 50, by = 1)    # Try a range of control sample sizes
control_orr <- 0.27                       # Control = historical = no conflict
treatment_orr <- control_orr + 0.2        # TV = 0.2

scenarios <- lapply(ctrl_n_values, function(ctrl_n) {
  trt_n <- 2 * ctrl_n

  tibble(
    total_n = ctrl_n + trt_n,
    trt_n = trt_n,
    ctrl_n = ctrl_n,
    trt_p = treatment_orr,
    ctrl_p = control_orr
  )
}) %>% bind_rows()


param_grid_h <- expand.grid(
  ctrl_h_n = 637,
  ctrl_h_p = 0.27,
  stringsAsFactors = FALSE
)

data_gen_params_list <- lapply(apply(scenarios, 1, as.list),
                               create_data_gen_params, endpoint = "binary")

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "binary")

data_gen_params_h <- data_gen_params_list_h[[1]]

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
  }

  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- sapply(data_gen_params, function(x) x$n)
  prob_list <- sapply(data_gen_params, function(x) x$p)

  # Generate simulation dataset
  data <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list, seed = 123)

  summary <- data$summary
  summary_h <- data.frame(control_h.n = param_grid_h$ctrl_h_n,
                          control_h.count = round(param_grid_h$ctrl_h_n * param_grid_h$ctrl_h_p))
  summary <- cbind(summary, summary_h)

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
    control.w = NA, # borrowing
    control.delta_gate = c(0.1),
    control.ess_h = c(1, 2) * unique(summary$control.n),
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

  # No borrowing
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
    control.w = 0,
    stringsAsFactors = FALSE
  )

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
                                      Lalonde_decision = T)

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

saveRDS(bayes_results, file = "simulations/casestudy_OR_bayes_results_1.rds")
bayes_results <- readRDS(file = "simulations/casestudy_OR_bayes_results.rds")

results <- process_sim_results(bayes_results)
results


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

saveRDS(bayes_results_all, "simulations/casestudy_OR_bayes_results_df_1.rds")


# + Visualization ---------------------------------------------------------

bayes_results_all <- readRDS("simulations/casestudy_OR_bayes_results_df_1.rds")

oc_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No")) %>%
  subset(is.na(control.delta_SAM) | control.delta_SAM != 0.15) %>%
  select("i" , "borrowing",
         "treatment.n", "control.n", "control.ess_h",
         "control.p", "control_h.p",
         "decision_pr", "proportion_pr", "proportion_ci") %>%
  mutate(decision_pr = factor(decision_pr, levels = c("no-go", "consider", "go"),
                              labels = c("No-Go", "Consider", "Go"))) %>%
  mutate(ess.ratio = ifelse(borrowing == "Borrowing: Yes", control.ess_h / control.n, 0))

oc_df <- oc_df %>%
  subset(ess.ratio != 1.5 &
           control.p == control_h.p &
           treatment.n == 2 * control.n &
           control.n >= 20) %>%
  mutate(n_ratios = factor(ess.ratio,
                           levels = c("0", "1", "2"),
                           labels = c(expression(n[hc*","*e] == 0 * " (No Borrowing)"),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:1')),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:2')))),


  )

p1 <- ggplot(oc_df ,
                   aes(x = control.n, y = proportion_pr, fill = decision_pr)) +
  geom_area(color = "black", linewidth = 0.1, position = "stack") +
  scale_fill_manual(
    values = c("No-Go" = "red", "Consider" = "#F0E442", "Go" = "#009E73"),
    name = "Decision"
  ) +
  facet_grid(~ n_ratios,
             labeller = label_parsed) +
  labs(y = "Decision Probability", x = expression("Concurrent Control Arm Sample Size (" * n[cc] * ")")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave("simulations/casestudy_zone_size_vs_sample_size.jpg", p1, width = 8, height = 3.5)

oc_subset <- oc_df %>%
  subset(decision_pr == "Consider" & proportion_pr <= 0.3) %>%
  arrange(control.n) %>%
  group_by(ess.ratio)


# 2. OC for Varying Treatment Effects (Without Conflict) ---------------------------------------------

treatment_orr <- seq(0.22, 0.52, by = 0.05)      # Try a range of treatment ORR
control_orr <- 0.27                       # Control = historical = no conflict

scenarios <- lapply(treatment_orr, function(trt_p) {
  tibble(
    trt_n = 60,
    ctrl_n = 30,
    trt_p = trt_p,
    ctrl_p = control_orr
  )
}) %>% bind_rows()


param_grid_h <- expand.grid(
  ctrl_h_n = 637,
  ctrl_h_p = 0.27,
  stringsAsFactors = FALSE
)


data_gen_params_list <- lapply(apply(scenarios, 1, as.list),
                               create_data_gen_params, endpoint = "binary")

data_gen_params_list_h <- lapply(apply(param_grid_h, 1, as.list),
                                 create_data_gen_params, endpoint = "binary")

data_gen_params_h <- data_gen_params_list_h[[1]]

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
  }

  n_arms <- length(data_gen_params)
  arm_names <- sapply(data_gen_params, function(x) x$name)
  n_list <- sapply(data_gen_params, function(x) x$n)
  prob_list <- sapply(data_gen_params, function(x) x$p)

  # Generate simulation dataset
  data <- data_gen_binary(n_arms = n_arms, arm_names = arm_names, nsim = nsim,
                          n_list = n_list, prob_list = prob_list, seed = 123)

  summary <- data$summary
  summary_h <- data.frame(control_h.n = param_grid_h$ctrl_h_n,
                          control_h.count = round(param_grid_h$ctrl_h_n * param_grid_h$ctrl_h_p))
  summary <- cbind(summary, summary_h)

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
    control.w = NA, # borrowing
    control.delta_gate = 0.1,
    control.ess_h = c(1, 2) * unique(summary$control.n),
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

  prior_params_list[[3]] <- list(
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
                                      Lalonde_decision = T)

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


saveRDS(bayes_results, file = "simulations/casestudy_OR_bayes_results_vary_Delta.rds")

bayes_results <- readRDS(file = "simulations/casestudy_OR_bayes_results_vary_Delta.rds")
results <- process_sim_results(bayes_results)

oc_df <- results$oc_all %>%
  mutate(ess.ratio = ifelse(borrow == "Yes", control.ess_h / control.n, 0)) %>%
  mutate(n_ratios = factor(ess.ratio,
                           levels = c("0", "1", "2"),
                           labels = c(expression(n[hc*","*e] == 0 * " (No Borrowing)"),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:1')),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:2')))))


p2 <- plot_oc(oc_df, x_var = "treatment.p", facet_formula = ~ n_ratios, plot_type = "smooth") +
  labs(y = "Decision Probability", x = expression("Treatment ORR (" * theta[t] * ")")) +
  scale_x_continuous(breaks= seq(0.22, 0.52, by = 0.05) )  +
  theme(
    legend.position = "bottom",
    title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )


ggsave("simulations/casestudy_zone_size_vs_trt_orr.jpg", p2, width = 8, height = 3.5)


# 3. Real Analysis (Without Conflict) -------------------------------------

summary <- expand.grid(
  control.n = 30,
  treatment.n = 60,
  treatment.count = 29,
  control.count = c(6, 8, 10),
  control_h.n = 637,
  control_h.count = 172
)
summary$nsim <- 1:nrow(summary)

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
  control.w = c(0,NA),
  control.delta_gate = 0.1,
  control.ess_h = c(1, 2) * unique(summary$control.n),
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


bayes_results <- list()
for (k in seq_along(prior_params_list)) {
  prior_params <- prior_params_list[[k]]

  settings <- prior_params

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


  post$settings <- settings
  post$data_summary <- cbind(i, summary)

  bayes_results[[length(bayes_results)+1]] <- post
}


post_params_all <- data.frame()
post_est_ci_all <- data.frame()
post_inference_all <- data.frame()
data_summary_all <- data.frame()
for (i in seq_along(bayes_results)) {
  res <- bayes_results[[i]]
  setting <- res$settings
  post_params <- cbind(setting, res$post_params)
  post_est_ci <- cbind(setting, res$post_est_ci)
  post_inference <- cbind(setting, res$post_inference)

  post_inference$decision_pr <- ifelse(is.na(post_inference$decision_pr), post_inference$decision_ci, post_inference$decision_pr)
  post_inference$decision_ci <- ifelse(is.na(post_inference$decision_ci), post_inference$decision_pr, post_inference$decision_ci)

  data_summary_all <- bind_rows(data_summary_all, res$data_summary)
  post_params_all <- bind_rows(post_params_all, cbind(i, post_params))
  post_est_ci_all <- bind_rows(post_est_ci_all, cbind(i, post_est_ci))
  post_inference_all <- bind_rows(post_inference_all, cbind(i, post_inference))
}


post_params_all <- post_params_all %>%
  mutate(borrow = ifelse(is.na(control.w), "Yes", "No"))

post_params_1 <- post_params_all[post_params_all$control.ess_h == 30, ]
post_params_2 <- post_params_all[post_params_all$control.ess_h == 60, ]

post_inference_1 <- post_inference_all[post_inference_all$control.ess_h == 30, ]
post_inference_2 <- post_inference_all[post_inference_all$control.ess_h == 60, ]

library(ggplot2)
library(dplyr)
library(tidyr)

#   bm <- RBesT::mixbeta(inf = c(post_params_borrow$control.w_post, post_params_borrow$control.a1_post, post_params_borrow$control.b1_post),
#                 inf2 = c(1 - post_params_borrow$control.w_post, post_params_borrow$control.a2_post, post_params_borrow$control.b2_post))
#
# library(RBesT)
#
# plot(bm)

plot_posterior <- function(post_params, post_inference, data_summary) {
  x_vals <- seq(0, 1, length.out = 1000)

  # Extract parameter values
  post_params_borrow <- post_params %>% filter(borrow == "Yes")
  post_params_noborrow <- post_params %>% filter(borrow == "No")

  # Compute densities
  df <- tibble(x = x_vals) %>%
    mutate(
      control_noborrow = dbeta(x, post_params_noborrow$control.a2_post, post_params_noborrow$control.b2_post),
      treatment = dbeta(x, post_params_noborrow$treatment.a2_post, post_params_noborrow$treatment.b2_post),
      control_borrow1 = dbeta(x, post_params_borrow$control.a1_post, post_params_borrow$control.b1_post),
      control_borrow2 = dbeta(x, post_params_borrow$control.a2_post, post_params_borrow$control.b2_post),
      control_mixture = post_params_borrow$control.w_post * control_borrow1 +
        (1 - post_params_borrow$control.w_post) * control_borrow2
    )

  # Reshape to long format for ggplot
  df_long <- df %>%
    select(x, control_noborrow, control_mixture, treatment) %>%
    pivot_longer(-x, names_to = "curve", values_to = "density") %>%
    mutate(curve = recode(curve,
                          control_noborrow = "Concurrent Control",
                          control_mixture = "Hybrid Control (SAM Prior)",
                          treatment = "Experimental Arm"))

  # Define fill colors with transparency
  fill_colors <- c(
    "Concurrent Control" = rgb(0, 0, 1, 0.3),
    "Hybrid Control (SAM Prior)" = rgb(1, 0, 0, 0.3),
    "Experimental Arm" = rgb(0, 1, 0, 0.3)
  )

  line_colors <- c(
    "Concurrent Control" = "blue",
    "Hybrid Control (SAM Prior)" = "red",
    "Experimental Arm" = "green"
  )

  # Get summary values
  post_inference <- post_inference %>%
    subset(is.na(control.w))
  data_summary <- unique(data_summary)
  n_trt <- data_summary$treatment.n
  y_trt <- data_summary$treatment.count
  n_ctrl <- data_summary$control.n
  y_ctrl <- data_summary$control.count
  n_ctrl_h <- data_summary$control_h.n
  y_ctrl_h <- data_summary$control_h.count
  decision <- post_inference$decision_pr

  prob_text <- paste0(
    "P(Minimal) = ", sprintf("%.3f", post_inference$pr_m),"\n",
    "P(Lower) = ", sprintf("%.3f", post_inference$pr_l), "\n",
    "P(Target) = ", sprintf("%.3f", post_inference$pr_t), "\n\n",
    "Decision: ", stringr::str_to_title(decision)
  )

  # response_text <- paste0(
  #   "Experimental: ", y_trt, "/", n_trt, "\n",
  #   "Concurrent Control: ", y_ctrl, "/", n_ctrl, "\n",
  #   "Historical Control: ", y_ctrl_h, "/", n_ctrl_h
  # )

  response_text <- paste0(
    "Concurrent Control: ", y_ctrl, "/", n_ctrl
    )

  # Plot
  ggplot(df_long, aes(x = x, y = density, fill = curve, color = curve)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line(size = 1.2) +
    scale_color_manual(values = line_colors) +
    scale_fill_manual(values = fill_colors) +
    labs(
      title = response_text,
      x = "ORR",
      y = "Density",
      color = "",
      fill = ""
    ) +
    annotate("text", x = 0.98, y = max(df_long$density) * 0.9,
             hjust = 1, vjust = 1, size = 4.5, fontface = "bold",
             label = prob_text, color = "black") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
}

post_params_1 <- post_params_all[post_params_all$control.ess_h == 30, ]
post_params_2 <- post_params_all[post_params_all$control.ess_h == 60, ]

post_inference_1 <- post_inference_all[post_inference_all$control.ess_h == 30, ]
post_inference_2 <- post_inference_all[post_inference_all$control.ess_h == 60, ]

p1_list <- list()
p2_list <- list()
for(i in 1:3){
  p1_list[[i]] <- plot_posterior(post_params = post_params_1[post_params_1$nsim == i,],
                                 post_inference = post_inference_1[post_inference_1$nsim == i,],
                                 data_summary = data_summary_all[data_summary_all$nsim == i,])

  p2_list[[i]] <- plot_posterior(post_params_2[post_params_2$nsim == i,],
                                 post_inference_2[post_inference_2$nsim == i,],
                                 data_summary = data_summary_all[data_summary_all$nsim == i,])

}

p1 <- wrap_plots(p1_list) +
  plot_layout(nrow = 1, guides = 'collect',
              axes = "collect") +
  plot_annotation(title = expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:1')),
                  theme = theme(plot.title = element_text(size = 18)))
p1 <- wrap_elements(p1)

p2 <- wrap_plots(p2_list) +
  plot_layout(nrow = 1, guides = 'collect',
              axes = "collect") +
  plot_annotation(title = expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:2')),
                  theme = theme(plot.title = element_text(size = 18)))

p2 <- wrap_elements(p2)

p_combined <- p1/p2

ggsave("simulations/casestudy_analysis.jpg", p_combined, width = 17, height = 8)




############


