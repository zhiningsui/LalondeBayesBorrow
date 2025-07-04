rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)
library(patchwork)

# Analysis ----------------------------------------------------------------

bayes_results <- readRDS("output_paper/sim_OR_bayes_results.rds") # Obtained from sim_OR.R
results <- process_sim_results(bayes_results)

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

# Generate Figure 1
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


top <- p_risk_noborrow + guide_area() +
  plot_layout(widths = c(1, 2), guides = 'collect')

p_risk <- (p_risk_borrow + theme(legend.position = "none")) / top  +
  plot_layout(heights = c(2, 1), guides = 'collect', axis_titles = 'collect')

ggsave("output_paper/sim_OR_conflict_vs_risk_paper.jpg", p_risk, width = 12, height = 8)


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

# Generate Figure 2
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

ggsave("output_paper/sim_OR_zone_size_paper.jpg", p_oc, width = 15, height = 22)



