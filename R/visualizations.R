#' Plot Operating Characteristics (OC)
#'
#' @description
#' Creates a stacked bar plot of the operating characteristics.
#'
#' @param oc_data A data frame, typically `processed_results$oc_all`.
#' @param x_var A character string for the x-axis variable.
#' @param facet_formula A formula for faceting the plot (e.g., `~ trt_n + borrowing`).
#' @param ... Additional arguments passed to `ggplot2` functions.
#'
#' @return A ggplot object.
#' @export
plot_oc <- function(oc_data, x_var, facet_formula = NULL, ...) {
  oc_data$decision_pr <- factor(oc_data$decision_pr, levels = c("no-go", "consider", "go"))

  p <- ggplot(oc_data, aes_string(x = x_var, y = "proportion_pr", fill = "decision_pr")) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = x_var, y = "Probability", fill = "Decision") +
    scale_fill_manual(values = c("no-go" = "red", "consider" = "#F0E442", "go" = "#009E73")) +
    theme_bw() +
    theme(legend.position = "bottom")

  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula)
  }

  p <- p + list(...)
  return(p)
}


#' Plot Risk Profile (Type I and Type II Errors)
#'
#' @description
#' Creates a plot of the False Go Rate (FGR) and False Stop Rate (FSR).
#'
#' @param oc_data A data frame, typically `processed_results$oc_all`.
#' @param x_var A character string for the x-axis variable.
#' @param facet_formula A formula for faceting the plot.
#' @param ... Additional arguments passed to `ggplot2` functions.
#'
#' @return A ggplot object.
#' @export
plot_risk <- function(oc_data, x_var, facet_formula = NULL, ...) {
  lrv_val <- min(abs(oc_data$true_value.compare_true))
  tv_val <- max(abs(oc_data$true_value.compare_true))

  risk_data <- oc_data %>%
    filter(
      (abs(true_value.compare_true) == lrv_val & decision_pr == "go") |
        (abs(true_value.compare_true) == tv_val & decision_pr == "no-go")
    ) %>%
    mutate(risk_type = ifelse(decision_pr == "go", "Type I Error (FGR)", "Type II Error (FSR)"))

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
