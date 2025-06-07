#' Plot Operating Characteristics (OC)
#'
#' @description
#' Creates a stacked bar plot of the operating characteristics.
#'
#' @param oc_data A data frame, typically from `run_simulation()`.
#' @param x_var A character string for the x-axis variable.
#' @param facet_formula A formula for faceting the plot (e.g., `~ trt_n + true_value.compare_true`).
#' @param ... Additional arguments passed to `ggplot2` functions.
#'
#' @return A ggplot object.
#' @export
plot_oc <- function(oc_data, x_var, facet_formula = NULL, colors = c("No-Go" = "red", "Consider" = "#F0E442", "Go" = "#009E73")) {

  oc_data$decision_pr <- factor(oc_data$decision_pr, levels = c("No-Go", "Consider", "Go"))

  p <- ggplot(oc_data, aes_string(x = x_var, y = "proportion_pr", fill = "decision_pr")) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = x_var, y = "Probability", fill = "Decision") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(legend.position = "bottom")

  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula)
  }

  # Add any extra ggplot arguments
  p <- p + list(...)

  return(p)
}

#' Plot Risk Profile (Type I and Type II Errors)
#'
#' @description
#' Creates a plot of the False Go Rate (FGR) and False Stop Rate (FSR).
#'
#' @param oc_data A data frame, typically from `run_simulation()`.
#' @param x_var A character string for the x-axis variable.
#' @param facet_formula A formula for faceting the plot (e.g., `~ trt_n`).
#' @param ... Additional arguments passed to `ggplot2` functions.
#'
#' @return A ggplot object.
#' @export
plot_risk <- function(oc_data, x_var, facet_formula = NULL, ...) {
  lrv_val <- min(abs(oc_data$true_value.compare_true))
  tv_val <- max(abs(oc_data$true_value.compare_true))

  risk_data <- oc_data %>%
    filter(
      (abs(true_value.compare_true) == lrv_val & decision_pr == "Go") |
        (abs(true_value.compare_true) == tv_val & decision_pr == "No-Go")
    ) %>%
    mutate(risk_type = ifelse(decision_pr == "Go", "Type I Error (FGR)", "Type II Error (FSR)"))

  p <- ggplot(risk_data, aes_string(x = x_var, y = "proportion_pr", color = "risk_type")) +
    geom_point(size = 3) +
    geom_hline(aes(yintercept = 0.2, color = "Type I Error (FGR)"), linetype = "dashed") +
    geom_hline(aes(yintercept = 0.1, color = "Type II Error (FSR)"), linetype = "dashed") +
    labs(x = x_var, y = "Risk Value", color = "") +
    theme_bw() +
    theme(legend.position = "bottom")

  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula)
  }

  p <- p + list(...)
  return(p)
}


#' Plot Borrowing Weights from SAM Prior
#'
#' @description
#' Creates boxplots of the prior and posterior borrowing weights.
#'
#' @param post_params A data frame, typically from `run_simulation()`.
#' @param x_var A character string for the x-axis variable.
#' @param facet_formula A formula for faceting the plots (e.g., `~ trt_n`).
#' @param ... Additional arguments passed to `ggplot2` functions.
#'
#' @return A patchwork object combining the two plots.
#' @export
plot_borrowing <- function(post_params, x_var, facet_formula = NULL, ...) {
  p1 <- ggplot(post_params, aes_string(x = x_var, y = "control.w_prior")) +
    geom_boxplot() +
    labs(y = "Prior Weight", title = "Prior Weight of Informative Component") +
    theme_bw()

  p2 <- ggplot(post_params, aes_string(x = x_var, y = "control.w_post")) +
    geom_boxplot() +
    labs(y = "Posterior Weight", title = "Posterior Weight of Informative Component") +
    theme_bw()

  if (!is.null(facet_formula)) {
    p1 <- p1 + facet_grid(facet_formula)
    p2 <- p2 + facet_grid(facet_formula)
  }

  p1 <- p1 + list(...)
  p2 <- p2 + list(...)

  p1 / p2 # Using patchwork to combine plots
}
