#' Process and Consolidate Simulation Results
#'
#' @description
#' This function takes a list of results from multiple simulation scenarios
#' (each being the output of `bayesian_lalonde_decision`) and consolidates them
#' into a final, tidy list of data frames ready for analysis and visualization.
#'
#' @param bayes_results A list where each element is a list of results from one
#'   call to `bayesian_lalonde_decision`.
#'
#' @return A list of consolidated data frames: `post_params_all`, `post_est_ci_all`,
#'   `post_inference_all`, and `oc_all`.
#'
#' @importFrom dplyr bind_rows group_by do ungroup
#' @export
process_sim_results <- function(bayes_results) {
  if (!is.list(bayes_results) || length(bayes_results) == 0) {
    warning("The 'bayes_results' list is empty. Returning NULL.")
    return(NULL)
  }

  final_results <- list(
    post_params_all = dplyr::bind_rows(lapply(bayes_results, `[[`, "post_params")),
    post_est_ci_all = dplyr::bind_rows(lapply(bayes_results, `[[`, "post_est_ci")),
    post_inference_all = dplyr::bind_rows(lapply(bayes_results, `[[`, "post_inference"))
  )

  setting_cols <- names(final_results$post_inference_all)[!grepl("nsim|pr_|decision_|ci_|est_|sd_", names(final_results$post_inference_all))]

  oc_all <- final_results$post_inference_all %>%
    group_by(across(all_of(setting_cols))) %>%
    do(obtain_oc(.)) %>%
    ungroup()

  final_results$oc_all <- oc_all

  return(final_results)
}
