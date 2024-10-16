#' Obtain Operating Characteristics (OC)
#'
#' This function calculates the operating characteristics (OC) by computing the proportion of decisions made based on the posterior probability and the credible interval from a set of decision results.
#' It combines these proportions with the provided settings data to return a comprehensive summary.
#'
#' @param settings A data frame. Contains the settings or parameters used in the simulations or trials.
#' @param decisions A data frame. Contains the decisions made in each simulation or trial, with columns for posterior probability-based decisions (`decision_pr`) and credible interval-based decisions (`decision_ci`).
#'
#' @return A data frame that includes the original settings along with the calculated proportions of decisions based on the posterior probability and credible interval. The output includes:
#'   - `proportion_pr`: Proportion of decisions made based on the posterior probability.
#'   - `proportion_ci`: Proportion of decisions made based on the credible interval.
#' @export
#'
#' @examples
#' # Example usage:
#' decisions <- data.frame(decision_pr = c("go", "no-go", "consider"),
#'                         decision_ci = c("go", "no-go", "consider"))
#' obtain_oc(decisions)
#'
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data
obtain_oc <- function(decisions){
  # Proportion of decisions made based on the posterior probability
  proportion_pr <- decisions %>%
    group_by(.data$decision_pr) %>%
    summarise(count_pr = n(), .groups = 'drop') %>%
    complete(decision_pr = factor(levels(.data$decision_pr)),
             fill = list(count_pr = 0)) %>%
    mutate(proportion_pr = .data$count_pr / sum(.data$count_pr, na.rm = TRUE))
  # Proportion of decisions made based on the credible interval
  proportion_ci <- decisions %>%
    group_by(.data$decision_ci) %>%
    summarise(count_ci = n(), .groups = 'drop') %>%
    complete(decision_ci = factor(levels(.data$decision_ci)),
             fill = list(count_ci = 0)) %>%
    mutate(proportion_ci = .data$count_ci / sum(.data$count_ci, na.rm = TRUE))

  OC <- full_join(proportion_pr, proportion_ci, by = c("decision_pr" = "decision_ci"))
  return(OC)
}
