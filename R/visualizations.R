#' Plot Operating Characteristics (OC)
#'
#' @description
#' Creates a stacked plot of the operating characteristics, showing the
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
#' @param plot_type A character string specifying the type of plot.
#'   Either "stepwise" for a bar plot with text or "smooth" for a stacked area plot.
#'   Defaults to "stepwise".
#' @param ... Additional arguments passed to `ggplot2` functions for further
#'   customization (e.g., `labs()`, `theme()`).
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar labs scale_fill_manual theme_bw theme facet_grid geom_text element_text label_parsed geom_area
#' @importFrom scales percent
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


#' Plot Posterior Distributions for a Single Run
#'
#' @description
#' This function visualizes and compares the posterior distributions for the
#' treatment arm, the control arm (without borrowing), and the hybrid control
#' arm (with borrowing) from a single simulation run.
#'
#' @param post_params A data frame with exactly two rows: one for a "borrowing"
#'   scenario and one for a "no borrowing" scenario from the same simulation run.
#' @param post_inference A data frame with the corresponding posterior inference
#'   results for the "borrowing" scenario. This should be filtered to a single row.
#' @param endpoint A character string specifying the endpoint type: "binary" or "continuous".
#' @param title_format A character string that defines the format for the plot title.
#'   You can use `{column_name}` syntax to insert values from the `post_params`
#'   data frame (e.g., `"Control: {control.count}/{control.n}"`). If `NULL`, a
#'   default title will be used.
#' @param x_range A numeric vector of length two specifying the range for the x-axis.
#'   If `NULL`, a reasonable default will be used.
#' @param ... Additional arguments passed to `ggplot2` functions for further
#'   customization (e.g., `labs()`, `theme()`).
#'
#' @return A ggplot object visualizing the posterior distributions.
#'
#' @importFrom dplyr %>% mutate select filter recode
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_area geom_line scale_color_manual scale_fill_manual labs annotate theme_minimal theme element_text
#' @importFrom stats dbeta dnorm
#' @importFrom glue glue
#' @importFrom grDevices rgb
#' @importFrom rlang `%||%`
#' @export
plot_posterior <- function(post_params, post_inference, endpoint, title_format = NULL, x_range = NULL, ...) {

  # --- Input Validation ---
  if (!is.data.frame(post_params) || nrow(post_params) != 2) {
    stop("'post_params' must be a data frame with exactly two rows: one for a 'borrowing' scenario and one for the 'no borrowing' scenario.")
  }
  if (!"borrow" %in% names(post_params)) {
    stop("'borrow' column not found in 'post_params'. It should contain 'Yes' and 'No' values or similar.")
  }
  if (sum(c("Yes", "No") %in% post_params$borrow) != 2) {
    stop("'post_params' must contain one row where borrow == 'Yes' and one row where borrow == 'No'.")
  }
  if (nrow(post_inference) != 1) {
    stop("Please provide the single row from 'post_inference_all' that corresponds to the borrowing scenario (e.g., borrow == 'Yes').")
  }
  if (!endpoint %in% c("binary", "continuous")) {
    stop("'endpoint' must be either 'binary' or 'continuous'.")
  }

  # --- Extract parameter values ---
  post_params_borrow <- post_params %>% filter(borrow == "Yes")
  post_params_noborrow <- post_params %>% filter(borrow == "No")

  # --- Determine plot range and compute densities ---
  if (endpoint == "binary") {
    x_vals <- seq(x_range[1] %||% 0, x_range[2] %||% 1, length.out = 1000)
    df <- tibble::tibble(x = x_vals) %>%
      mutate(
        control_noborrow = dbeta(x, post_params_noborrow$control.a2_post, post_params_noborrow$control.b2_post),
        treatment = dbeta(x, post_params_noborrow$treatment.a2_post, post_params_noborrow$treatment.b2_post),
        control_borrow1 = dbeta(x, post_params_borrow$control.a1_post, post_params_borrow$control.b1_post),
        control_borrow2 = dbeta(x, post_params_borrow$control.a2_post, post_params_borrow$control.b2_post),
        control_mixture = post_params_borrow$control.w_post * control_borrow1 +
          (1 - post_params_borrow$control.w_post) * control_borrow2
      )
    x_lab <- "Response Rate"
    response_text <- paste0("Concurrent Control: ", post_params_noborrow$control.count, "/", post_params_noborrow$control.n)
  } else { # continuous
    if(is.null(x_range)) {
      all_means <- c(post_params_borrow$control.mu1_post, post_params_borrow$control.mu2_post, post_params_noborrow$treatment.mu2_post)
      all_sds <- c(post_params_borrow$control.sigma1_post, post_params_borrow$control.sigma2_post, post_params_noborrow$treatment.sigma2_post)
      min_x <- min(all_means - 3 * all_sds, na.rm = TRUE)
      max_x <- max(all_means + 3 * all_sds, na.rm = TRUE)
      x_range <- c(min_x, max_x)
    }
    x_vals <- seq(x_range[1], x_range[2], length.out = 1000)
    df <- tibble::tibble(x = x_vals) %>%
      mutate(
        control_noborrow = dnorm(x, post_params_noborrow$control.mu2_post, post_params_noborrow$control.sigma2_post),
        treatment = dnorm(x, post_params_noborrow$treatment.mu2_post, post_params_noborrow$treatment.sigma2_post),
        control_borrow1 = dnorm(x, post_params_borrow$control.mu1_post, post_params_borrow$control.sigma1_post),
        control_borrow2 = dnorm(x, post_params_borrow$control.mu2_post, post_params_borrow$control.sigma2_post),
        control_mixture = post_params_borrow$control.w_post * control_borrow1 +
          (1 - post_params_borrow$control.w_post) * control_borrow2
      )
    x_lab <- "Value"
    response_text <- paste0("Concurrent Control: n=", post_params_noborrow$control.n)
  }

  # --- Prepare data for plotting ---
  df_long <- df %>%
    select(x, control_noborrow, control_mixture, treatment) %>%
    pivot_longer(-x, names_to = "curve", values_to = "density") %>%
    mutate(curve = recode(curve,
                          control_noborrow = "Concurrent Control",
                          control_mixture = "Hybrid Control (SAM Prior)",
                          treatment = "Experimental Arm"))

  # --- Define plot aesthetics ---
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

  # --- Create annotation text ---
  prob_text <- paste0(
    "P(Minimal) = ", sprintf("%.3f", post_inference$pr_m),"\n",
    "P(Lower) = ", sprintf("%.3f", post_inference$pr_l), "\n",
    "P(Target) = ", sprintf("%.3f", post_inference$pr_t), "\n\n",
    "Decision: ", stringr::str_to_title(post_inference$decision_pr)
  )

  # --- Create the plot ---
  p <- ggplot(df_long, aes(x = x, y = density, fill = curve, color = curve)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line(size = 1.2) +
    scale_color_manual(values = line_colors) +
    scale_fill_manual(values = fill_colors) +
    labs(
      title = response_text,
      x = x_lab,
      y = "Density",
      color = "",
      fill = ""
    ) +
    annotate("text", x = Inf, y = Inf, label = prob_text, hjust = 1.1, vjust = 1.1,
             size = 4.5, fontface = "bold", color = "black") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      legend.position = "bottom"
    )

  # Apply any additional ggplot2 customizations
  p <- p + list(...)

  return(p)
}

