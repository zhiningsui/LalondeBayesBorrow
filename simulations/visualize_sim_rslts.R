rm(list = ls())
library(LalondeBayesBorrow)
library(tidyverse)
library(kableExtra)


# OR endpoint -------------------------------------------------------------

bayes_results_all <- readRDS("simulations/sim_OR_bayes_results_df_5_23.rds")

# + Generate summary table --------------------------------------

post_inference_all <-  bayes_results_all$post_inference_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  select("nsim", "borrowing", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p",
         "est2_lalonde") %>%
  pivot_wider(names_from = borrowing, values_from = est2_lalonde) %>%
  mutate(pmd = `Borrowing: Yes` - `Borrowing: No`) %>%
  select(-c(`Borrowing: Yes`, `Borrowing: No`, true_value.compare_true))

# Summary statistics of PMD
metrics_tab <- post_inference_all %>%
  subset(gate == "Gated Control: Yes" & control.delta_SAM != 0.05) %>%
  group_by(control.delta_SAM, control.ess_h, control.p) %>%
  summarize(mean_pmd = mean(pmd),
            sd_pmd = sd(pmd),
            .groups = "drop")

metrics_tab_1 <- metrics_tab %>%
  pivot_longer(c(mean_pmd, sd_pmd)) %>%
  pivot_wider(names_from = c(control.delta_SAM, control.ess_h),
              values_from = value, names_sort = T) %>%
  mutate(name = factor(name, levels = c("mean_pmd", "sd_pmd"))) %>%
  select(c(2,1,3:8)) %>%
  arrange(name, control.p)


kbl(metrics_tab_1,
    escape = F,
    format = "latex",
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

risk_df <- bayes_results_all$oc_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.p.diff = control.p - control_h.p) %>%
  filter(true_value.compare_true %in% c(0.1, 0.2),
         decision_pr %in% c("go", "no-go")) %>%
  mutate(
    CorrectLabel = case_when(
      true_value.compare_true == 0.2 & decision_pr == "go" ~ "CGR",
      true_value.compare_true == 0.1 & decision_pr == "no-go" ~ "CSR",
      TRUE ~ NA_character_
    ),
    RiskLabel = case_when(
      true_value.compare_true == 0.1 & decision_pr == "go" ~ "FGR",
      true_value.compare_true == 0.2 & decision_pr == "no-go" ~ "FSR",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(cols = c(CorrectLabel, RiskLabel), names_to = "type", values_to = "label") %>%
  filter(!is.na(label)) %>%
  group_by(gate, borrowing, control.delta_SAM, control.ess_h, control.p, label) %>%
  summarize(prop = mean(proportion_pr), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = prop)


risk_df1 <- risk_df %>%
  pivot_longer(c(FSR, FGR)) %>%
  subset(gate == "Gated Control: Yes" & control.delta_SAM !=0.05) %>%
  mutate(control.p = factor(control.p,
                            levels = c("0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5")),
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
  control.delta_SAM = factor(c(0.1, 0.15),
                             levels = c("0.1", "0.15"),
                             labels = c(expression(paste(delta, " = 0.1")),
                                        expression(paste(delta, " = 0.15")))),
  xmin = c(3, 2),
  xmax = c(7, 8),
  ymin = -Inf,
  ymax = Inf)

p_risk_borrow <- ggplot(risk_df1 %>% filter(borrowing == "Borrowing: Yes"),
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
  subset(borrowing == "Borrowing: No") %>%
  select(-c(control.delta_SAM, control.ess_h)) %>%
  distinct() %>%
  mutate(control.ess_h = factor(0, labels = expression(n[ch*","*e] == 0 * " (No Borrowing)"))) %>%
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

# p_risk <- top / (p_risk_borrow + theme(legend.position = "none"))  +
#   plot_layout(heights = c(1, 2), guides = 'collect', axis_titles = 'collect')

p_risk <- (p_risk_borrow + theme(legend.position = "none")) / top  +
  plot_layout(heights = c(2, 1), guides = 'collect', axis_titles = 'collect')


ggsave( "simulations/sim_OR_conflict_vs_risk_paper.jpg", p_risk, width = 12, height = 8)

# + Generate plots for OC -------------------------------------------------

oc_df <- bayes_results_all$oc_all %>%
  subset(control.delta_gate %in% c(0.1, 0.15)) %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         control.p.diff = control.p - control_h.p) %>%
  select("i", "borrowing",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.p.diff",
         "decision_pr", "proportion_pr", "proportion_ci") %>%
  mutate(decision_pr = factor(decision_pr, levels = c("no-go", "consider", "go"),
                              labels = c("No-Go", "Consider", "Go")),
         true_value.compare_true = ifelse(true_value.compare_true == 0.2, "0.2 (TV)", "0.1 (LRV)"))


oc_new <- data.frame()
oc_tmp <- na.omit(oc_df)
for (i in seq(1, nrow(oc_tmp), 3)) {
  tmp <- oc_tmp[i:(i+2),]
  tmp[tmp$decision_pr == "Go", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] * 0.995
  tmp[tmp$decision_pr == "Consider", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] + tmp[tmp$decision_pr == "Consider", "proportion_pr"] + 0.02
  tmp[tmp$decision_pr == "No-Go", "label_ypos"] <- 1 + 0.012
  oc_new <- rbind(oc_new, tmp)
}

oc_new <- oc_new[, c("true_value.compare_true", "control.delta_SAM", "control.ess_h",
                     "borrowing",
                     "control.p.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]

oc_borrow <- oc_new %>% subset(borrowing == "Borrowing: Yes")
oc_noborrow <- oc_new %>%
  subset(borrowing == "Borrowing: No") %>%
  select(-c(control.delta_SAM, control.ess_h)) %>%
  distinct() %>%
  mutate(control.delta_SAM = NA,
         control.ess_h = 0)

oc_new <- bind_rows(oc_borrow, oc_noborrow)
oc_new$control.ess_h <- factor(as.character(oc_new$control.ess_h),
                               levels = c("0", "45", "90", "180"),
                               labels = c(
                                 expression(n[ch*","*e] == 0 * " (No Borrowing)"),
                                 expression(n[ch*","*e] == 45),
                                 expression(n[ch*","*e] == 90),
                                 expression(n[ch*","*e] == 180)
                               ))

oc_new$control.p.diff <- factor(oc_new$control.p.diff,
                                levels = c("-0.2", "-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15", "0.2"))
oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
oc_new$control.delta_SAM <- factor(oc_new$control.delta_SAM,
                                          levels = c("0.1", "0.15"),
                                          labels = c(expression(paste(delta, " = 0.1")),
                                                     expression(paste(delta, " = 0.15"))))

myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_borrow$decision_pr)

oc_borrow <- subset(oc_new, borrowing == "Borrowing: Yes")
oc_noborrow <- subset(oc_new, borrowing == "Borrowing: No")

p_list <- list()
for (i in levels(oc_borrow$true_value.compare_true)) {

  oc_tmp_noborrow <- oc_noborrow %>% filter(true_value.compare_true == i)
  oc_tmp_borrow <- oc_borrow %>% filter(true_value.compare_true == i)

  oc_tmp_borrow1 <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[2])
  oc_tmp_borrow2a <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[3])
  oc_tmp_borrow2b <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[4])

  p_noborrow <- ggplot(oc_tmp_noborrow,
                       aes(x = control.p.diff, y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
              vjust = 1.6, fontface = "bold", size = 4.2) +
    scale_fill_manual(name = "Decision", values = myColors) +
    facet_grid( ~ control.ess_h,
               labeller = labeller(control.ess_h = label_parsed)) +
    labs(x = expression(Delta[c] == theta[c] - theta[ch]), y = "Probability") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15))

  make_plot <- function(data, legend_pos = "none") {
    ggplot(data, aes(x = control.p.diff, y = proportion_pr, fill = decision_pr)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
                vjust = 1.6, fontface = "bold", size = 4.2) +
      scale_fill_manual(name = "Decision", values = myColors) +
      facet_grid(control.delta_SAM ~ control.ess_h,
                 labeller = labeller(control.ess_h = label_parsed,
                                     control.delta_SAM = label_parsed)) +
      labs(x = expression(Delta[c] == theta[c] - theta[ch]), y = "Probability") +
      theme_bw(base_size = 14) +
      theme(legend.position = legend_pos,
            strip.text = element_text(size = 16),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 15))
  }

  p_borrow1 <- make_plot(oc_tmp_borrow1, legend_pos = "bottom")
  p_borrow2a <- make_plot(oc_tmp_borrow2a)
  p_borrow2b <- make_plot(oc_tmp_borrow2b)

  layout <- "
  AC
  AC
  BC
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

ggsave("simulations/sim_OR_zone_size_paper.jpg", p_oc, width = 15, height = 21)




# DpR endpoint ------------------------------------------------------------

bayes_results_all <- readRDS("simulations/sim_DpR_unconstrained_bayes_results_df_5_31_sigma40.rds")

# + Generate summary table --------------------------------------

post_inference_all <-  bayes_results_all$post_inference_all %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         gate = ifelse(is.na(control.delta_gate), "Gated Control: No", "Gated Control: Yes"),
         control.mu.diff = control.mu - control_h.mu) %>%
  select("nsim", "borrowing", "gate",
         "true_value.compare_true",
         "control.delta_SAM", "control.ess_h",
         "control.mu",
         "est2_lalonde") %>%
    distinct() %>%
  subset(gate == "Gated Control: Yes") %>%
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
  subset(gate == "Gated Control: Yes" & control.delta_SAM != 0.05) %>%
  pivot_wider(names_from = c(borrowing, control.delta_SAM, control.ess_h),
              values_from = value, names_sort = T) %>%
  mutate(name = factor(name, levels = c("FGR", "FSR", "CGR", "CSR", "mean_pmd", "sd_pmd"))) %>%
  select(c(3,2,10:15)) %>%
  arrange(name, control.mu) %>%
  subset(name %in% c("mean_pmd", "sd_pmd"))


kbl(metrics_tab_1,
    escape = F,
    format = "latex",
    row.names = FALSE,
    digits = 4,
    col.names = c("Metric", "control.mu",
                  rep(c("45", "90", "180"), 2))) %>%
  kable_styling(bootstrap_options = c("condensed"), full_width = FALSE) %>%
  kable_classic(full_width = FALSE) %>%
  add_header_above(c(" " = 1, " " = 1,
                     "Borrowing with delta = 0.1" = 3,
                     "Borrowing with delta = 0.15" = 3)) %>%
  collapse_rows(columns = 1, valign = "middle")


# + Visualize risks -------------------------------------------------------

risk_df1 <- risk_df %>%
  pivot_longer(c(FSR, FGR)) %>%
  subset(gate == "Gated Control: Yes" & control.delta_SAM !=0.05) %>%
  mutate(control.mu = factor(control.mu,
                             levels = c("-0.4", "-0.35", "-0.3", "-0.25", "-0.2", "-0.15", "-0.1", "-0.05", "0")),
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
  control.delta_SAM = factor(c(0.1, 0.15),
                             levels = c("0.1", "0.15"),
                             labels = c(expression(paste(delta, " = 0.1")),
                                        expression(paste(delta, " = 0.15")))),
  xmin = c(3, 2),
  xmax = c(7, 8),
  ymin = -Inf,
  ymax = Inf)

p_risk_borrow <- ggplot(risk_df1 %>% filter(borrowing == "Borrowing: Yes"),
                        aes(x = control.mu, y = value, color = name)) +
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
  xlab(~ paste("Concurrent Control DpR Mean (", theta[c], ")")) +
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
  subset(borrowing == "Borrowing: No") %>%
  select(-c(control.delta_SAM, control.ess_h)) %>%
  distinct() %>%
  mutate(control.delta_SAM = NA,
         control.ess_h = factor(0, labels = expression(n[ch*","*e] == 0 * " (No Borrowing)"))) %>%
  ggplot(aes(x = control.mu, y = value, color = name)) +
  geom_point(size = 2) +
  facet_grid( ~ control.ess_h,
              labeller = labeller(control.ess_h = label_parsed)) +
  labs(y = "Risk Value",
       color = "") +
  xlab(~ paste("Concurrent Control DpR Mean (", theta[c], ")")) +
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

# p_risk <- top / (p_risk_borrow + theme(legend.position = "none"))  +
#   plot_layout(heights = c(1, 2), guides = 'collect', axis_titles = 'collect')

p_risk <- (p_risk_borrow + theme(legend.position = "none")) / top  +
  plot_layout(heights = c(2, 1), guides = 'collect', axis_titles = 'collect')


ggsave( "simulations/sim_DpR_unconstrained_conflict_vs_risk_paper.jpg", p_risk, width = 12, height = 8)


# + Generate plots for OC -------------------------------------------------

oc_df <- bayes_results_all$oc_all %>%
  subset(control.delta_gate %in% c(0.1, 0.15)) %>%
  mutate(borrowing = ifelse(is.na(control.w), "Borrowing: Yes", "Borrowing: No"),
         control.mu.diff = control.mu - control_h.mu) %>%
  select("i", "borrowing",
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
  tmp[tmp$decision_pr == "Go", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] * 0.995
  tmp[tmp$decision_pr == "Consider", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] + tmp[tmp$decision_pr == "Consider", "proportion_pr"] + 0.02
  tmp[tmp$decision_pr == "No-Go", "label_ypos"] <- 1 + 0.012
  oc_new <- rbind(oc_new, tmp)
}

oc_new <- oc_new[, c("true_value.compare_true", "control.delta_SAM", "control.ess_h",
                     "borrowing",
                     "control.mu.diff",
                     "decision_pr", "proportion_pr", "label_ypos")]

oc_borrow <- oc_new %>% subset(borrowing == "Borrowing: Yes")

oc_noborrow <- oc_new %>%
  subset(borrowing == "Borrowing: No") %>%
  select(-c(control.delta_SAM, control.ess_h)) %>%
  distinct() %>%
  mutate(control.delta_SAM = NA,
         control.ess_h = 0)

oc_new <- bind_rows(oc_borrow, oc_noborrow)
oc_new$control.ess_h <- factor(as.character(oc_new$control.ess_h),
                               levels = c("0", "45", "90", "180"),
                               labels = c(
                                 expression(n[ch*","*e] == 0 * " (No Borrowing)"),
                                 expression(n[ch*","*e] == 45),
                                 expression(n[ch*","*e] == 90),
                                 expression(n[ch*","*e] == 180)
                               ))

oc_new$control.mu.diff <- factor(oc_new$control.mu.diff,
                                    levels = c("-0.2", "-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15", "0.2"))
oc_new$true_value.compare_true <- factor(oc_new$true_value.compare_true)
oc_new$control.delta_SAM <- factor(oc_new$control.delta_SAM,
                                   levels = c("0.1", "0.15"),
                                   labels = c(expression(paste(delta, " = 0.1")),
                                              expression(paste(delta, " = 0.15"))))

myColors <- c("red","#F0E442","#009E73")
names(myColors) <- levels(oc_borrow$decision_pr)

oc_borrow <- subset(oc_new, borrowing == "Borrowing: Yes")
oc_noborrow <- subset(oc_new, borrowing == "Borrowing: No")

p_list <- list()
for (i in levels(oc_borrow$true_value.compare_true)) {

  oc_tmp_noborrow <- oc_noborrow %>% filter(true_value.compare_true == i)
  oc_tmp_borrow <- oc_borrow %>% filter(true_value.compare_true == i)

  oc_tmp_borrow1 <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[2])
  oc_tmp_borrow2a <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[3])
  oc_tmp_borrow2b <- oc_tmp_borrow %>% filter(control.ess_h == levels(oc_new$control.ess_h)[4])

  p_noborrow <- ggplot(oc_tmp_noborrow,
                       aes(x = control.mu.diff, y = proportion_pr, fill = decision_pr)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
              vjust = 1.6, fontface = "bold", size = 4.2) +
    scale_fill_manual(name = "Decision", values = myColors) +
    facet_grid( ~ control.ess_h,
                labeller = labeller(control.ess_h = label_parsed)) +
    labs(x = expression(Delta[c] == theta[c] - theta[ch]), y = "Probability") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15))

  make_plot <- function(data, legend_pos = "none") {
    ggplot(data, aes(x = control.mu.diff, y = proportion_pr, fill = decision_pr)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
                vjust = 1.6, fontface = "bold", size = 4.2) +
      scale_fill_manual(name = "Decision", values = myColors) +
      facet_grid(control.delta_SAM ~ control.ess_h,
                 labeller = labeller(control.ess_h = label_parsed,
                                     control.delta_SAM = label_parsed)) +
      labs(x = expression(Delta[c] == theta[c] - theta[ch]), y = "Probability") +
      theme_bw(base_size = 14) +
      theme(legend.position = legend_pos,
            strip.text = element_text(size = 16),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15),
            axis.title = element_text(size = 17),
            axis.text = element_text(size = 15))
  }

  p_borrow1 <- make_plot(oc_tmp_borrow1, legend_pos = "bottom")
  p_borrow2a <- make_plot(oc_tmp_borrow2a)
  p_borrow2b <- make_plot(oc_tmp_borrow2b)

  layout <- "
  AC
  AC
  BC
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

ggsave("simulations/sim_DpR_unconstrained_zone_size_paper.jpg", p_oc, width = 15, height = 21)







