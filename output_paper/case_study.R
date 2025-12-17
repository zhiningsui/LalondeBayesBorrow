rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)
library(ggplot2)
library(dplyr)
library(scales)

# 1. Sample Size Determination (No Conflict) ------------------------------------

ctrl_n_values <- seq(20, 50, by = 1)      # Try a range of control sample sizes
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

saveRDS(bayes_results, file = "output_paper/casestudy_OR_bayes_results_1.rds")

bayes_results <- readRDS(file = "output_paper/casestudy_OR_bayes_results_1.rds")

results <- process_sim_results(bayes_results)

# + Visualization ---------------------------------------------------------


oc_df <- results$oc_all %>%
  mutate(ess.ratio = ifelse(borrow == "Yes", control.ess_h / control.n, 0),
         n_ratios = factor(ess.ratio,
                           levels = c("0", "1", "2"),
                           labels = c(expression(n[hc*","*e] == 0 * " (No Borrowing)"),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:1')),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:2')))),


  )

p1 <- plot_oc(oc_data = oc_df, x_var = "control.n", facet_formula = ~ n_ratios, plot_type = "smooth") +
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

ggsave("output_paper/casestudy_zone_size_vs_sample_size.jpg", p1, width = 8, height = 3.5)

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

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}
make_key <- function(i, meth) sprintf("design%03d__%s", i, safe_name(meth))


for (i in seq_along(data_gen_params_list)) {
  start_time_i <- Sys.time()

  data_gen_params <- data_gen_params_list[[i]]
  cat("\n===========================================================\n")
  cat("Design set", i, "\n")
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
                          n_list = n_list, prob_list = prob_list, seed = 1000)

  summary <- data$summary %>%
    mutate(control_h.n = data_gen_params_h$control_h$n,
           control_h.count = round(data_gen_params_h$control_h$n * data_gen_params_h$control_h$p))

  cat("Start Bayesian evaluation analysis on this design...\n")
  prior_grid <- expand.grid(
    method = "SAM",
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
    control.ess_h = c(1, 3, 5, 7) * unique(summary$control.n),
    stringsAsFactors = FALSE
  ) %>%
    mutate(control.delta_SAM = control.delta_gate)

  # prior_grid <- expand.grid(
  #   method = "SAM",
  #   treatment.a0 = 1,
  #   treatment.b0 = 1,
  #   treatment.a = 1,
  #   treatment.b = 1,
  #   treatment.w = 0,
  #   control.a0 = 1,
  #   control.b0 = 1,
  #   control.a = 1,
  #   control.b = 1,
  #   control.w = NA, # borrowing
  #   control.delta_gate = 0.1,
  #   control.ess_h = c(1, 2) * unique(summary$control.n),
  #   stringsAsFactors = FALSE
  # ) %>%
  #   mutate(control.delta_SAM = control.delta_gate)


  # Convert to list and replace NA with NULL
  prior_params_list <- apply(prior_grid, 1, function(row) {
    out <- as.list(row)
    if (is.na(out$control.w)) out$control.w <- NULL
    if (is.na(out$control.delta_gate)) out$control.delta_gate <- NULL
    out
  })
  names(prior_params_list) <- paste0("gated_SAM.", seq_along(prior_params_list))

  prior_params_list[["NoBorrow"]] <- list(
    method = "SAM",
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
    meth <- names(prior_params_list)[k]
    key <- make_key(i, meth)

    settings <- c(list(true_value = data$true_value),
                  data_gen_params, data_gen_params_h,
                  prior_params)

    settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE)
    settings_num <- unlist(apply(settings, 2, num_or_null))
    settings <- as.data.frame(c(settings[!names(settings) %in% names(settings_num)], settings_num))

    # settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE) %>%
    #   mutate(across(
    #     .cols = -c(ends_with(".name")),
    #     .fns = as.numeric
    #   ))

    post <- bayesian_lalonde_decision(
      endpoint = "binary",
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
                                    calc_post_dist_metrics(
                                      endpoint = "OR",
                                      true_value = unique(post$post_est_ci$true_value.compare_true),
                                      post_est_ci = post$post_est_ci))

    end_time_k <- Sys.time()
    post$elapsed_time_secs <- as.numeric(difftime(end_time_k, start_time_k, units = "secs"))

    bayes_results[[key]] <- post
    cat("Completed Bayesian analysis: ", key, "\n\n")
  }
  end_time_i <- Sys.time()
  cat("Total time for design set", i, "=", round(as.numeric(difftime(end_time_i, start_time_i, units="secs")), 2), "seconds\n")
}


saveRDS(bayes_results, file = "casestudy_OR_bayes_results_2_1017.rds")

bayes_results <- readRDS("casestudy_OR_bayes_results_2_1017.rds")
results <- process_sim_results(bayes_results)

oc_df <- results$oc_all  %>%
  mutate(ess.ratio = ifelse(borrow == "Yes", control.ess_h / control.n, 0),
         n_ratios = factor(ess.ratio,
                           levels = c("0", "1", "3", "5", "7"),
                           labels = c(expression(n[hc*","*e] == 0 * " (No Borrowing)"),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:1')),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:3')),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:5')),
                                      expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ' = 2:1:7'))))) %>%
  filter(ess.ratio != 7)


p2 <- plot_oc(oc_df, x_var = "treatment.p", facet_formula = ~ n_ratios, plot_type = "smooth") +
  labs(y = "Decision Probability", x = expression("Treatment ORR (" * theta[t] * ")")) +
  scale_x_continuous(breaks= seq(0.22, 0.52, by = 0.05) )  +
  theme(
    legend.position = "right",
    # title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )


ggsave("casestudy_zone_size_vs_trt_orr_1017.jpg", p2, width = 9.5, height = 3.5)


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
  method = "SAM",
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
  control.ess_h = c(1, 3, 5) * unique(summary$control.n),
  stringsAsFactors = FALSE
) %>%
  mutate(control.delta_SAM = control.delta_gate)

prior_params_list <- apply(prior_grid, 1, function(row) {
  out <- as.list(row)
  if (is.na(out$control.w)) out$control.w <- NULL
  if (is.na(out$control.delta_gate)) out$control.delta_gate <- NULL
  out
})
names(prior_params_list) <- paste0("gated_SAM.", seq_along(prior_params_list))

prior_params_list[["NoBorrow"]] <- list(
  method = "SAM",
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


bayes_results <- list()
for (k in seq_along(prior_params_list)) {

  prior_params <- prior_params_list[[k]]
  meth <- names(prior_params_list)[k]
  # key <- make_key(i, meth)

  settings <- prior_params

  settings <- as.data.frame(t(unlist(settings)), stringsAsFactors = FALSE)
  settings_num <- unlist(apply(settings, 2, num_or_null))
  settings <- as.data.frame(c(settings[!names(settings) %in% names(settings_num)], settings_num))

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

  bayes_results[[length(bayes_results)+1]] <- post
}


results <- process_sim_results(bayes_results)

post_params_all <- results$post_params_all
post_inference_all <- results$post_inference_all

plot_list <- list()
for (ess in unique(post_params_all$control.ess_h)) {
  if (is.na(ess)) next
  post_params <- post_params_all %>%
    filter(control.ess_h == ess | is.na(control.ess_h))
  post_inference <- post_inference_all %>%
    filter(control.ess_h == ess | is.na(control.ess_h))

  p_list <- list()
  for(i in 1:3){
    p_list[[i]] <- plot_posterior(post_params = post_params[post_params$nsim == i,],
                                   post_inference = post_inference[post_inference$nsim == i & post_inference$borrow == "Yes", ],
                                   endpoint = "binary",
                                   title_format = "Concurrent Control: {control.count}/{control.n}") +
      guides(colour = guide_legend(ncol = 1))

  }
  ratio = paste0(" = 2:1:", ess / 30)
  p1 <- wrap_plots(p_list) +
    plot_layout(nrow = 1, guides = 'collect',
                axes = "collect") +
    plot_annotation(
      # title = expression(paste(n[t], ' : ', n[cc], ' : ', n[hc*','*e], ratio)),
      title = bquote(n[t] : n[cc] : n[hc*","*e] ~ .(ratio)),
      theme = theme(plot.title = element_text(size = 18)))

  p1 <- wrap_elements(p1)

  plot_list[[as.character(ess)]] <- p1
}

p_combined <- wrap_elements(plot_list$`30`/plot_list$`90`/plot_list$`150`)

ggsave("casestudy_analysis_1017.jpg", p_combined, width = 17, height = 16)



fig_real_A <- wrap_elements((plot_spacer() + p2 + plot_spacer()) +
                      plot_layout(widths = c(0.1, 1.3, 0.2)))

design <- "AAAAAAAA###
           BBBBBBBBBB#
           BBBBBBBBBB#
           BBBBBBBBBB#
           BBBBBBBBBB#
           BBBBBBBBBB#
           BBBBBBBBBB#
           BBBBBBBBBB#"

fig_real <- fig_real_A / p_combined +
  plot_annotation(tag_levels = 'A') +
  plot_layout(
    # design = design
    heights = c(1.1,3.7)) &
  theme(plot.tag = element_text(size = 25))

ggsave("casestudy_analysis_new.jpg", fig_real, width = 16, height = 14)



design <- "
AAAAAA##
BBBBBBBB
BBBBBBBB
BBBBBBBB
BBBBBBBB
BBBBBBBB
BBBBBBBB
"

fig_real <-
  wrap_elements(p2 + plot_annotation(tag = "A")) /
  wrap_elements(p_combined + plot_annotation(tag = "B")) +
  plot_layout(design = design) &  # collect & align legends
  theme(
    # tag styling and alignment (top-left of each block)
    plot.tag = element_text(size = 25, face = "bold"),
    plot.tag.position = c(0, 1),  # x=0 (left), y=1 (top)
    # legend alignment (shared because of guides='collect')
    legend.position = "right",    # or "bottom"
    legend.justification = c(0, 1),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )



############


