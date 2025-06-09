rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)

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

saveRDS(bayes_results, file = "simulations/sim_OR_bayes_results_new.rds")


results <- process_sim_results(bayes_results)

saveRDS(results, "simulations/sim_OR_bayes_results_df_new.rds")

# Analysis ----------------------------------------------------------------

results <- readRDS("simulations/sim_OR_bayes_results_df_new.rds")

# + Generate summary table --------------------------------------

# Create the summary table for PMD
pmd_summary <- create_pmd_summary(post_inference_all = results$post_inference_all)

# Create Table 2
pmd_tab <- pmd_summary %>%
  pivot_longer(c(mean_pmd, sd_pmd)) %>%
  select(-control.delta_gate) %>%
  pivot_wider(names_from = c(control.delta_SAM, control.ess_h),
              values_from = value, names_sort = T) %>%
  mutate(name = factor(name, levels = c("mean_pmd", "sd_pmd"))) %>%
  select(c(12, 3, 13:18)) %>%
  arrange(name, control.p)

kbl(pmd_tab,
    escape = F,
    # format = "latex",
    row.names = FALSE,
    digits = 4,
    col.names = c("Metric", "control.p",
                  rep(c("45", "90", "180"), 2))) %>%
  kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 1, " " = 1,
                     "Borrowing with delta = 0.1" = 3,
                     "Borrowing with delta = 0.15" = 3)) %>%
  collapse_rows(columns = 1, valign = "middle")


# + Visualize risks -------------------------------------------------------

risk_df <- create_risk_df(oc_all = results$oc_all, lrv = 0.1, tv = 0.2)

risk_df1 <- risk_df %>%
  pivot_longer(c(FSR, FGR)) %>%
  mutate(control.p = factor(control.p,
                            levels = c("0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5")),
         control.delta_SAM = factor(control.delta_SAM,
                                    levels = c("0.1", "0.15"),
                                    labels = c(expression(paste(delta, " = 0.1")),
                                               expression(paste(delta, " = 0.15")))),
         control.ess_h = factor(control.ess_h,
                                levels = c("45", "90", "180"),
                                labels = c(expression(n[hc*","*e] == 45),
                                           expression(n[hc*","*e] == 90),
                                           expression(n[hc*","*e] == 180))),
         borrow = paste0("Borrowing: ", borrow))

shade_df <- data.frame(
  control.delta_SAM = factor(c(0.1, 0.15),
                             levels = c("0.1", "0.15"),
                             labels = c(expression(paste(delta, " = 0.1")),
                                        expression(paste(delta, " = 0.15")))),
  xmin = c(3, 2),
  xmax = c(7, 8),
  ymin = -Inf,
  ymax = Inf)

p_risk_borrow <- ggplot(risk_df1 %>% filter(borrow == "Borrowing: Yes"),
                        aes(x = control.p, y = value, color = name)) +
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
       color = "") +
  xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
  geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
  )

p_risk_noborrow <- risk_df1 %>%
  subset(borrow == "Borrowing: No") %>%
  select(-c(control.delta_SAM, control.ess_h)) %>%
  distinct() %>%
  mutate(control.ess_h = factor(0, labels = expression(n[hc*","*e] == 0 * " (No Borrowing)"))) %>%
  ggplot(aes(x = control.p, y = value, color = name)) +
  geom_point(size = 2) +
  facet_grid( ~ control.ess_h,
              labeller = labeller(control.ess_h = label_parsed)) +
  labs(y = "Risk Value",
       color = "") +
  xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
  geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
  )

library(patchwork)

top <- p_risk_noborrow + guide_area() +
  plot_layout(widths = c(1, 2), guides = 'collect')

p_risk <- (p_risk_borrow + theme(legend.position = "none")) / top  +
  plot_layout(heights = c(2, 1), guides = 'collect', axis_titles = 'collect')

ggsave("simulations/sim_OR_conflict_vs_risk_paper.jpg", p_risk, width = 12, height = 8)

# + Generate plots for OC -------------------------------------------------

oc_data <- results$oc_all %>%
  mutate(borrow = paste0("Borrowing: ", borrow)) %>%
  select("borrow",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p",
         "decision_pr", "proportion_pr") %>%
  mutate(true_value.compare_true = ifelse(true_value.compare_true == 0.2, "0.2 (TV)", "0.1 (LRV)"))


oc_borrow <- oc_data %>% subset(borrow == "Borrowing: Yes")
oc_noborrow <- oc_data %>%
  subset(borrow == "Borrowing: No") %>%
  select(-c(control.delta_SAM, control.ess_h)) %>%
  distinct() %>%
  mutate(control.delta_SAM = NA,
         control.ess_h = 0)

oc_new <- bind_rows(oc_borrow, oc_noborrow)
oc_new$control.ess_h <- factor(as.character(oc_new$control.ess_h),
                               levels = c("0", "45", "90", "180"),
                               labels = c(
                                 expression(n[hc*","*e] == 0 * " (No Borrowing)"),
                                 expression(n[hc*","*e] == 45),
                                 expression(n[hc*","*e] == 90),
                                 expression(n[hc*","*e] == 180)
                               ))

oc_new$control.p <- factor(oc_new$control.p,
                           levels = c("0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5"))
oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
oc_new$control.delta_SAM <- factor(oc_new$control.delta_SAM,
                                   levels = c("0.1", "0.15"),
                                   labels = c(expression(paste(delta, " = 0.1")),
                                              expression(paste(delta, " = 0.15"))))

oc_borrow <- subset(oc_new, borrow == "Borrowing: Yes")
oc_noborrow <- subset(oc_new, borrow == "Borrowing: No")

p_list <- list()
for (i in levels(oc_borrow$true_value.compare_true)) {

  oc_tmp_noborrow <- oc_noborrow %>% filter(true_value.compare_true == i)
  oc_tmp_borrow <- oc_borrow %>% filter(true_value.compare_true == i)

  oc_tmp_borrow1 <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[2])
  oc_tmp_borrow2a <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[3])
  oc_tmp_borrow2b <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[4])

  p_noborrow <- plot_oc(oc_data = oc_tmp_noborrow, x_var = "control.p",
                        facet_formula = ~ control.ess_h ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")"))


  make_plot <- function(data, legend_pos = "none") {
    ggplot(data, aes(x = control.p, y = proportion_pr, fill = decision_pr)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
                vjust = 1.6, fontface = "bold", size = 4.2) +
      scale_fill_manual(name = "Decision", values = myColors) +
      facet_grid(control.delta_SAM ~ control.ess_h,
                 labeller = labeller(control.ess_h = label_parsed,
                                     control.delta_SAM = label_parsed)) +
      labs(y = "Probability") +
      xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
      theme_bw(base_size = 14) +
      theme(legend.position = legend_pos,
            strip.text = element_text(size = 16),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 15))
  }

  p_borrow1 <-  plot_oc(oc_data = oc_tmp_borrow1, x_var = "control.p",
                        facet_formula = control.delta_SAM ~ control.ess_h ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")"))

  p_borrow2a <- plot_oc(oc_data = oc_tmp_borrow2a, x_var = "control.p",
                        facet_formula = control.delta_SAM ~ control.ess_h ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
    theme(legend.position = "none")

  p_borrow2b <- plot_oc(oc_data = oc_tmp_borrow2b, x_var = "control.p",
                        facet_formula = control.delta_SAM ~ control.ess_h ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
    theme(legend.position = "none")

  layout <- "
  AC
  AC
  AC
  BC
  BC
  DE
  DE
  DE
  DE
  "
  p_combined <- p_noborrow + guide_area() + p_borrow1 + p_borrow2a + p_borrow2b +
    plot_layout(design = layout, guides = "collect") +
    plot_annotation(title = bquote(Delta == .(i))) &
    theme(plot.title = element_text(size = 23))

  p_list[[i]] <- wrap_elements(full = p_combined)
}

library(patchwork)
p_oc <- wrap_plots(p_list, ncol = 1)

ggsave("simulations/sim_OR_zone_size_paper.jpg", p_oc, width = 15, height = 22)



