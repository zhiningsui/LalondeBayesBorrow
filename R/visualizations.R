#' Plot Operating Characteristics (OC)
#'
#' @description
#' Creates a stacked bar plot of the operating characteristics, showing the
#' probabilities of "Go", "No-Go", and "Consider" decisions for each
#' simulation scenario.
#'
#' @param oc_data A data frame containing the operating characteristics,
#'   typically from `process_sim_results()`. Must contain columns for the
#'   decision proportions (e.g., `proportion_pr`) and the decision types
#'   (e.g., `decision_pr`).
#' @param x_var A character string specifying the variable for the x-axis
#'   (e.g., `"control.n"` or `"control.p.diff"`).
#' @param facet_formula A formula for faceting the plot, to split it into
#'   panels based on other setting variables (e.g., `~ borrowing` or
#'   `true_value.compare_true ~ borrowing`).
#' @param ... Additional arguments passed to `ggplot2` functions for further
#'   customization (e.g., `labs()`, `theme()`).
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar labs scale_fill_manual theme_bw theme facet_grid
#' @export
#'
#' @examples
#' # Assuming 'processed_results' is the output from process_sim_results()
#' # and contains an 'oc_all' data frame.
#'
#' # plot_oc(
#' #   oc_data = processed_results$oc_all,
#' #   x_var = "control.n",
#' #   facet_formula = true_value.compare_true ~ borrowing,
#' #   labs(
#' #     title = "Operating Characteristics by Sample Size and Borrowing",
#' #     x = "Sample Size per Arm",
#' #     y = "Probability of Decision"
#' #   )
#' # )
plot_oc <- function(oc_data, x_var, facet_formula = NULL, plot_type = "stepwise", ...) {
  
  if ("decision_pr" %in% names(oc_data)) {
    oc_data$decision_pr <- factor(oc_data$decision_pr, levels = c("No-Go", "Consider", "Go"))
    y_var <- "proportion_pr"
    fill_var <- "decision_pr"
  } else {
    stop("A decision proportion column (e.g., 'proportion_pr') was not found in the data.")
  }

  oc_new <- data.frame()
  oc_tmp <- oc_data
  for (i in seq(1, nrow(oc_tmp), 3)) {
    tmp <- oc_tmp[i:(i+2),]
    tmp[tmp$decision_pr == "Go", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"]
    tmp[tmp$decision_pr == "Consider", "label_ypos"] <- tmp[tmp$decision_pr == "Go", "proportion_pr"] + tmp[tmp$decision_pr == "Consider", "proportion_pr"]
    tmp[tmp$decision_pr == "No-Go", "label_ypos"] <- 1
    oc_new <- rbind(oc_new, tmp)
  }


  if (plot_type == "stepwise") {
    p <- ggplot(oc_new, aes_string(x = x_var, y = y_var, fill = fill_var)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(y = label_ypos, label = scales::percent(proportion_pr, 0.01)),
                vjust = 1.6, fontface = "bold", size = 4.2) 
  } else if (plot_type == "smooth") {
    p <- ggplot(oc_new, aes_string(x = x_var, y = y_var, fill = fill_var)) +
      geom_area(color = "black", linewidth = 0.1, position = "stack") 
  }
  p <- p +
    scale_fill_manual(name = "Decision",
                      values = c("No-Go" = "red", "Consider" = "#F0E442", "Go" = "#009E73")) +
    labs(y = "Probability") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 16),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15))

  # Apply faceting if a formula is provided
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, labeller = label_parsed)
  }

  # Apply any additional ggplot2 customizations
  p <- p + list(...)

  return(p)
}

