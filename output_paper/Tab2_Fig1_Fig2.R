rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)
library(patchwork)

# Analysis ----------------------------------------------------------------

bayes_results <- readRDS("sim_OR_bayes_results_1015.rds") # Obtained from sim_OR.R
bayes_results <- bayes_results[sort(names(bayes_results))]

results <- process_sim_results(bayes_results)

# + Generate summary table --------------------------------------

post_inference_all <- results$post_inference_all %>%
  mutate(method = ifelse(is.na(control.delta_gate), "Naive-SAM", method))

post_inference_all_unique <- post_inference_all[, !duplicated(as.list(post_inference_all))]

pmd_results <- create_pmd_summary(post_inference_all_unique)

# Prior parameters
priors_param_tab <- pmd_results$priors %>%
  mutate(method = ifelse(prior_name == "NoBorrow", "NoBorrow", method),
         method = factor(method, levels = c("SAM", "Naive-SAM", "DPP", "NoBorrow")),
         control.delta = ifelse(is.na(control.delta_gate), control.delta_SAM, control.delta_gate)) %>%
  arrange(method, control.delta, control.ess_h) %>%
  select(prior_name, control.delta, control.ess_h) %>%
  remove_rownames()

prior_order <- priors_param_tab$prior_name
priors_param_tab

# Data generation scenarios
scn_params <- pmd_results$scenarios

scn_params_pmd <- scn_params[grepl("control", colnames(scn_params))] %>%
  distinct()


# PMD table
pmd_tab <- pmd_results$pmd_summary %>%
  select(design_scenario.control, prior_name, mean_pmd, sd_pmd) %>%
  pivot_longer(c(mean_pmd, sd_pmd),
               names_to = "metric") %>%
  left_join(scn_params_pmd) %>%
  mutate(prior_name = factor(prior_name, levels = prior_order)) %>%
  pivot_wider(names_from = prior_name,
              values_from = value,
              names_sort = TRUE
  ) %>%
  rename(control.p = true_value.control.rate_true) %>%
  select(where(~ n_distinct(.) > 1)) %>%
  select(-design_scenario.control)

# Risks
risk_df <- create_risk_df(oc_all = results$oc_all, lrv = 0.1, tv = 0.2) %>%
  select(prior_name, control.p, FSR, FGR, CSR, CGR) %>%
  pivot_longer(c(FSR, FGR, CSR, CGR), names_to = "metric") %>%
  na.omit()

# FGR and FSR table
risk_tab <- risk_df %>%
  filter(metric %in% c("FSR", "FGR")) %>%
  mutate(prior_name = factor(prior_name, levels = prior_order)) %>%
  pivot_wider(
    names_from = prior_name,
    values_from = value,
    names_sort = T
  )

# CGR and CSR table
priors_param_tab1 <- priors_param_tab %>%
  mutate(
    method = factor(sub("\\..*$", "", prior_name),
                    levels = c("gated_SAM", "naive_SAM", "DPP", "NoBorrow"),
                    labels = c("Proposed SAM", "Naive SAM", "DPP", "No Borrowing"))
  ) %>%
  mutate(
    method2 = paste0(method, "_", control.ess_h)
  )

risk_tab2 <- risk_df %>%
  filter(metric %in% c("CSR", "CGR")) %>%
  left_join(priors_param_tab1) %>%
  mutate(
    control.delta = factor(control.delta,
                           levels = c("0.1", "0.15"),
                           labels = c(expression(paste(delta, " = 0.1")),
                                      expression(paste(delta, " = 0.15")))),
    control.ess_h = factor(control.ess_h,
                           levels = c("0", "45", "90", "180"),
                           labels = c(expression(n[hc*","*e] == 0),
                                      expression(n[hc*","*e] == 45),
                                      expression(n[hc*","*e] == 90),
                                      expression(n[hc*","*e] == 180)))
  ) %>%
  select(where(~ n_distinct(.) > 1))


ggplot(risk_tab2 %>% filter(method == "Proposed SAM"),
       aes(x = control.p, y = value,
           color = control.ess_h, linetype = control.delta)) +
  geom_line() +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "solid")) +
  facet_grid( ~ metric,
              labeller = labeller(control.ess_h = label_parsed,
                                  control.delta = label_parsed)) +
  labs(y = "Risk Value",
       color = "")




ggplot(risk_tab2,
       aes(x = control.p, y = value, color = control.ess_h, linetype = method)) +
  geom_line() +
  facet_grid(control.delta ~ metric,
             labeller = labeller(control.ess_h = label_parsed,
                                 control.delta = label_parsed)) +
  labs(y = "Risk Value",
       color = "")

xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
  geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
  geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
  coord_cartesian(ylim = c(ymin, ymax)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13)
  )



# Calibrated tau table
tau_tab <- pmd_results$priors %>%
  select(prior_name, tau_TV, tau_LRV) %>%
  distinct() %>%
  mutate(prior_name = factor(prior_name, levels = prior_order)) %>%
  pivot_longer(c(tau_TV, tau_LRV),
               names_to = "metric") %>%
  pivot_wider(names_from = prior_name,
              values_from = value,
              names_sort = TRUE) %>%
  mutate(control.p = NA)

# MSE table
mse_summary <- post_inference_all_unique %>%
  mutate(sq_err = (true_value.compare_true - est_compare_lalonde)^2) %>%
  group_by(prior_name, true_value.compare_true, true_value.control.rate_true) %>%
  summarise(
    mse = mean(sq_err),
    .groups = "drop"
  )

mse_tab <- mse_summary %>%
  mutate(prior_name = factor(prior_name, levels = prior_order)) %>%
  pivot_wider(
    names_from = prior_name,
    values_from = mse,
    names_sort = T
  ) %>%
  mutate(metric = paste0("mse_", true_value.compare_true)) %>%
  rename(control.p = true_value.control.rate_true) %>%
  arrange(control.p)


# Table 2
tab_2 <- bind_rows(tau_tab, pmd_tab, risk_tab, mse_tab) %>%
  select(metric, control.p, all_of(prior_order)) %>%
  mutate(metric = factor(metric, levels = c("tau_TV", "tau_LRV", "FSR", "FGR", "mean_pmd", "sd_pmd", "mse_0.1", "mse_0.2"))) %>%
  arrange(metric, control.p)

kbl(tab_2,
    escape = T,
    format = "latex",
    row.names = FALSE,
    digits = 4,
    format.args = list(nsmall = 4, scientific = FALSE),
    col.names = c("Metric", "control.p",
                  priors_param_tab[priors_param_tab$prior_name %in% prior_order, "control.ess_h"]),
    booktabs = T) %>%
  kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 1,
                     " " = 1,
                     "delta = 0.1" = 3,
                     "delta = 0.15" = 3,
                     "delta = 0.1" = 1,
                     "delta = 0.15" = 1,
                     "delta = 0.1" = 3,
                     "delta = 0.15" = 3,
                     " " = 1)) %>%
  add_header_above(c(" " = 1,
                     " " = 1,
                     "Proposed SAM" = 6,
                     "Naive SAM" = 2,
                     "DPP" = 6,
                     "No Borrow" = 1)) %>%
  collapse_rows(columns = 1, valign = "middle")






#
#
# # + Visualize risks -------------------------------------------------------
#
#
# # Generate Figure 1
# risk_df1 <- risk_df %>%
#   left_join(priors_param_tab %>% rownames_to_column("prior_name")) %>%
#   mutate(
#     method = factor(sub("\\..*$", "", prior_name),
#                     levels = c("gated_SAM", "naive_SAM", "DPP", "NoBorrow"),
#                     labels = c("Proposed SAM", "Naive SAM", "DPP", "No Borrowing")),
#     control.p = factor(control.p,
#                        levels = c("0.1", "0.15", "0.2", "0.25",
#                                   "0.3", "0.35", "0.4", "0.45", "0.5")),
#     control.delta = factor(control.delta,
#                            levels = c("0.1", "0.15"),
#                            labels = c(expression(paste(delta, " = 0.1")),
#                                       expression(paste(delta, " = 0.15")))),
#     control.ess_h = ifelse(method == "Naive SAM", "180",
#                            ifelse(method == "No Borrowing", "0", control.ess_h)),
#     control.ess_h = factor(control.ess_h,
#                            levels = c("0", "45", "90", "180"),
#                            labels = c(expression(n[hc*","*e] == 0),
#                                       expression(n[hc*","*e] == 45),
#                                       expression(n[hc*","*e] == 90),
#                                       expression(n[hc*","*e] == 180)))
#
#     ) %>%
#   select(where(~ n_distinct(.) > 1))
#
#
#
#
#
# shade_df <- data.frame(
#   control.delta = factor(c(0.1, 0.15),
#                              levels = c("0.1", "0.15"),
#                              labels = c(expression(paste(delta, " = 0.1")),
#                                         expression(paste(delta, " = 0.15")))),
#   xmin = c(3, 2),
#   xmax = c(7, 8),
#   ymin = -Inf,
#   ymax = Inf)
#
# ymin <- floor(min(risk_df1$value, na.rm = TRUE) * 100) / 100
# ymax <- ceiling(max(risk_df1$value, na.rm = TRUE) * 100) / 100
#
# p_risk_borrow1 <- ggplot(
#   risk_df1 %>% filter(method %in% c("Proposed SAM", "Naive SAM")),
#   aes(x = control.p, y = value, color = Metric)
# ) +
#   geom_point(size = 2) +
#   geom_rect(
#     data = shade_df,
#     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
#         group = control.delta),
#     inherit.aes = FALSE,
#     fill = "blue", alpha = 0.1
#   ) +
#   facet_grid(control.delta ~ method + control.ess_h,
#              labeller = labeller(control.ess_h = label_parsed,
#                                  control.delta = label_parsed)) +
#   labs(y = "Risk Value",
#        color = "") +
#   xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
#   geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
#   geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
#   coord_cartesian(ylim = c(ymin, ymax)) +
#   theme_bw() +
#   theme(legend.position = "none",
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 13)
#   )
#
#
# p_risk_borrow2 <- ggplot(
#   risk_df1 %>% filter(method %in% c("DPP")),
#   aes(x = control.p, y = value, color = Metric)
# ) +
#   geom_point(size = 2) +
#   geom_rect(
#     data = shade_df,
#     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
#         group = control.delta),
#     inherit.aes = FALSE,
#     fill = "blue", alpha = 0.1
#   ) +
#   facet_grid(control.delta ~ method + control.ess_h,
#              labeller = labeller(control.ess_h = label_parsed,
#                                  control.delta = label_parsed)) +
#   labs(y = "Risk Value",
#        color = "") +
#   xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
#   geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
#   geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
#   coord_cartesian(ylim = c(ymin, ymax)) +
#   theme_bw()+
#   theme(legend.position = "none",
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 13)
#   )
#
# p_risk_noborrow <- ggplot(
#   risk_df1 %>% filter(method == "No Borrowing"),
#   aes(x = control.p, y = value, color = Metric)
# ) +
#   geom_point(size = 2) +
#   facet_grid( ~ method + control.ess_h,
#               labeller = labeller(control.ess_h = label_parsed)) +
#   labs(y = "Risk Value",
#        color = "") +
#   xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
#   geom_hline(aes(yintercept = 0.2, color = "FGR"), linetype = "dashed") +  # Horizontal line for false go risk
#   geom_hline(aes(yintercept = 0.1, color = "FSR"), linetype = "dashed") +
#   coord_cartesian(ylim = c(ymin, ymax)) +
#   theme_bw()+
#   theme(legend.position = "bottom",
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 13)
#   )
#
# p_risk_noborrow1 <- (p_risk_noborrow / guide_area()) +
#   plot_layout(heights = c(1.45, 1),
#               guides = 'collect')
#
# p_top <- p_risk_borrow1 + theme()
# p_bottom <- p_risk_borrow2 + p_risk_noborrow1 +
#   plot_layout(widths = c(3,1), guides = 'collect', axis_titles = 'collect')
#
# p_risk <- p_top / p_bottom +
#   plot_layout(heights = c(1, 1))
#
# ggsave("sim_OR_conflict_vs_risk_paper.jpg", p_risk, width = 15, height = 10)
#

# + Generate plots for OC -------------------------------------------------

oc_data <- results$oc_all %>%
  select(design_scenario, prior_name, borrow, decision_pr, proportion_pr) %>%
  left_join(scn_params) %>%
  left_join(priors_param_tab) %>%
  mutate(prior_name = factor(prior_name, levels = prior_order),
         control.ess_h_code = ifelse(prior_name == "NoBorrow", 0,
                                     ifelse(str_starts(prior_name, "naive_SAM"), 180, control.ess_h)),
         model_prefix = case_when(
           str_starts(prior_name, "gated_SAM") ~ "gSAM: ",
           str_starts(prior_name, "naive_SAM") ~ "Naive SAM: ",
           str_starts(prior_name, "DPP")       ~ "DPP-EB: ",
           TRUE                                    ~ ""
         ),
         # facet_ess_label = case_when(
         #   prior_name == "NoBorrow" ~ '"No Borrowing"',
         #   control.ess_h_code %in% c(45, 90, 180) ~
         #     paste0('paste("', model_prefix, '", n[hc*","*e] == ', control.ess_h_code, ')'),
         #   TRUE ~ '"Unknown"'
         # ),
         facet_ess_label = case_when(
           prior_name == "NoBorrow" ~ 'bold("No Borrowing")',
           control.ess_h_code %in% c(45, 90, 180) ~
             paste0('paste(bold("', model_prefix, '"), n[hc*","*e] == ', control.ess_h_code, ')'),
           TRUE ~ '"Unknown"'
         ),
         borrow = paste0("Borrowing: ", borrow),
         true_value.compare_true = factor(ifelse(true_value.compare_true == 0.2, "0.2 (TV)", "0.1 (LRV)")),
         true_value.control.rate_true = factor(true_value.control.rate_true,
                                               levels = c("0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5")),
         control.delta = factor(control.delta,
                                levels = c("0.1", "0.15"),
                                labels = c(expression(paste(delta, " = 0.1")),
                                           expression(paste(delta, " = 0.15"))))
  )




oc_plots <- list()
for (i in levels(oc_data$true_value.compare_true)) {
  oc_list <- list()
  oc_noborrow <- oc_data %>% filter(prior_name == "NoBorrow", true_value.compare_true == i, )

  p_noborrow <- plot_oc(oc_data = oc_noborrow, x_var = "true_value.control.rate_true",
                        facet_formula = ~ facet_ess_label ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")"))

  oc_list[["NoBorrow"]] <- p_noborrow

  for (nm in c("gated_SAM", "naive_SAM", "DPP")) {
    oc_borrow <- oc_data %>% filter(true_value.compare_true == i, str_starts(prior_name, nm))

    for (ess in sort(unique(oc_borrow$control.ess_h_code))) {
      oc_borrow1 <- oc_borrow %>% filter(control.ess_h_code == ess)

      p_borrow1 <-  plot_oc(oc_data = oc_borrow1, x_var = "true_value.control.rate_true",
                            facet_formula = control.delta ~ facet_ess_label ) +
        xlab(~ paste("Concurrent Control ORR (", theta[c], ")"))

      oc_list[[paste0(nm, "_", ess)]] <- p_borrow1

    }
  }
  oc_plots[[i]] <- oc_list
}

p_list <- list()
for (i in names(oc_plots)) {

  oc_lrv <- oc_plots[[i]]

  layout <- "
  ACDE
  ACDE
  ACDE
  BCDE
  BCDE
  FGHI
  FGHI
  FGHI
  FGHI
  FGHI
  "

  p_combined <- oc_lrv$NoBorrow + guide_area() + oc_lrv$gated_SAM_45 + oc_lrv$gated_SAM_90 + oc_lrv$gated_SAM_180 +
    oc_lrv$naive_SAM_180 + oc_lrv$DPP_45 + oc_lrv$DPP_90 + oc_lrv$DPP_180 +
    plot_layout(design = layout, guides = "collect") +
    plot_annotation(
      title = paste0("\u0394 = ", i),
      theme = theme(plot.title = element_text(size = 30, face = "bold"))
    )

  p_list[[i]] <- wrap_elements(full = p_combined)
}

p_list <- lapply(p_list, function(p) p + theme(plot.margin = margin(10, 0, 30, 0)))
p_oc <- wrap_plots(p_list, ncol = 1)

p_oc <- wrap_plots(p_list, ncol = 1)

ggsave("sim_OR_zone_size_paper_all.jpg", p_oc, width = 28, height = 28)


p_list <- list()
for (i in names(oc_plots)) {

  oc_lrv <- oc_plots[[i]]

  layout <- "
  AB
  CD
  "

  p_combined <- oc_lrv$naive_SAM_180 + oc_lrv$DPP_45 + oc_lrv$DPP_90 + oc_lrv$DPP_180 +
    plot_layout(design = layout, guides = "collect") +
    plot_annotation(
      title = paste0("\u0394 = ", i),
      theme = theme(plot.title = element_text(size = 30, face = "bold"),
                    legend.position = "bottom")
    )

  p_list[[i]] <- wrap_elements(full = p_combined)
}

p_list <- lapply(p_list, function(p) p + theme(plot.margin = margin(10, 0, 30, 0)))

p_oc <- wrap_plots(p_list, ncol = 1)

ggsave("sim_OR_zone_size_paper_other_methods.jpg", p_oc, width = 15, height = 24)



p_list <- list()
for (i in names(oc_plots)) {

  oc_lrv <- oc_plots[[i]]

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
  DE
  DE
  "

  p_combined <- oc_lrv$NoBorrow + guide_area() + oc_lrv$gated_SAM_45 + oc_lrv$gated_SAM_90 + oc_lrv$gated_SAM_180 +
    plot_layout(design = layout, guides = "collect") +
    plot_annotation(
      title = paste0("\u0394 = ", i),
      theme = theme(plot.title = element_text(size = 25, face = "bold"))
    )

  p_list[[i]] <- wrap_elements(full = p_combined)
}

p_oc <- wrap_plots(p_list, ncol = 1)

ggsave("sim_OR_zone_size_paper_gatedSAM.jpg", p_oc, width = 15, height = 24)



# + Generate plots for OC -------------------------------------------------

oc_borrow <- subset(oc_data, borrow == "Borrowing: Yes")
oc_noborrow <- subset(oc_data, borrow == "Borrowing: No")

# Generate Figure 2
# Proposed SAM
p_list <- list()
for (i in levels(oc_borrow$true_value.compare_true)) {
  oc_noborrow <- oc_noborrow %>% filter(true_value.compare_true == i, )
  oc_borrow <- oc_borrow %>% filter(true_value.compare_true == i)

  oc_borrow1 <- oc_borrow %>% filter(control.ess_h == levels(oc_data$control.ess_h)[2])
  oc_borrow2a <- oc_borrow %>% filter(control.ess_h == levels(oc_data$control.ess_h)[3])
  oc_borrow2b <- oc_borrow %>% filter(control.ess_h == levels(oc_data$control.ess_h)[4])

  p_noborrow <- plot_oc(oc_data = oc_noborrow, x_var = "true_value.control.rate_true",
                        facet_formula = ~ control.ess_h ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")"))

  p_borrow1 <-  plot_oc(oc_data = oc_borrow1, x_var = "true_value.control.rate_true",
                        facet_formula = control.delta_SAM ~ control.ess_h ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")"))

  p_borrow2a <- plot_oc(oc_data = oc_borrow2a, x_var = "true_value.control.rate_true",
                        facet_formula = control.delta_SAM ~ control.ess_h ) +
    xlab(~ paste("Concurrent Control ORR (", theta[c], ")")) +
    theme(legend.position = "none")

  p_borrow2b <- plot_oc(oc_data = oc_borrow2b, x_var = "true_value.control.rate_true",
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

ggsave("output_paper/sim_OR_zone_size_paper.jpg", p_oc, width = 15, height = 22)



