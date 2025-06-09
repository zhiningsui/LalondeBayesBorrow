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

  final_results$post_params_all <- final_results$post_params_all %>%
    mutate(borrow = ifelse(is.na(control.w), "Yes",
                           ifelse(control.w == 0, "No", paste0("Fixed weight = ", control.w) ))) %>%
    select(-control.w)

  final_results$post_est_ci_all <- final_results$post_est_ci_all %>%
    mutate(borrow = ifelse(is.na(control.w), "Yes",
                           ifelse(control.w == 0, "No", paste0("Fixed weight = ", control.w) ))) %>%
    select(-control.w)

  final_results$post_inference_all <- final_results$post_inference_all %>%
    mutate(borrow = ifelse(is.na(control.w), "Yes",
                           ifelse(control.w == 0, "No", paste0("Fixed weight = ", control.w) ))) %>%
    select(-control.w)

  start_idx <- which(names(final_results$post_params_all) == "nsim")
  end_idx <- which(names(final_results$post_params_all) == "control.w_prior")
  all_cols <- 1:ncol(final_results$post_params_all)
  setting_cols <- all_cols[!all_cols %in% c(start_idx: end_idx)]
  final_results$post_params_all <- final_results$post_params_all[, c(setting_cols, c(start_idx: end_idx))]

  start_idx <- which(names(final_results$post_est_ci_all) == "nsim")
  end_idx <- which(names(final_results$post_est_ci_all) == "compare_ci_u_95ci")
  all_cols <- 1:ncol(final_results$post_est_ci_all)
  setting_cols <- all_cols[!all_cols %in% c(start_idx: end_idx)]
  final_results$post_est_ci_all <- final_results$post_est_ci_all[, c(setting_cols, c(start_idx: end_idx))]

  start_idx <- which(names(final_results$post_inference_all) == "nsim")
  end_idx <- which(names(final_results$post_inference_all) == "decision_ci")
  all_cols <- 1:ncol(final_results$post_inference_all)
  setting_cols <- all_cols[!all_cols %in% c(start_idx: end_idx)]
  final_results$post_inference_all <- final_results$post_inference_all[, c(setting_cols, c(start_idx: end_idx))]

  idx <- which(names(final_results$post_inference_all) == "nsim")
  oc_all <- final_results$post_inference_all %>%
    group_by(across(all_of(1: (idx-1) ))) %>%
    do(obtain_oc(.)) %>%
    ungroup()

  final_results$oc_all <- oc_all

  return(final_results)
}


#' Create a Data Frame of Posterior Mean Differences (PMD)
#'
#' @description
#' This function calculates the Posterior Mean Difference (PMD) by comparing any
#' two specified borrowing scenarios from a simulation run. It correctly handles
#' both one-to-one comparisons (e.g., 'SAM Prior' vs. 'Fixed Weight') and
#' one-to-many comparisons against a common baseline (e.g., comparing multiple
#' borrowing scenarios to a single 'No Borrowing' run).
#'
#' @param post_inference_all A data frame containing the posterior inference results
#'   from a simulation run, typically from `process_sim_results()`.
#' @param borrow_to_compare A character vector of length two specifying the two borrowing
#'   conditions from the 'borrow' column to compare (e.g., c("Yes", "No")). The
#'   difference will be calculated as `borrow_to_compare[1] - borrow_to_compare[2]`.
#' @param pmd_var The name of the single posterior estimate variable to use for
#'   calculating the PMD. Defaults to `"est2_lalonde"`.
#'
#' @return A data frame summarizing the mean and standard deviation of the PMD for
#'   each relevant combination of settings.
#'
#' @importFrom dplyr %>% group_by summarize filter select mutate inner_join left_join distinct n_distinct across all_of
#' @importFrom tidyr pivot_wider
#' @export
create_pmd_summary <- function(post_inference_all,
                               borrow_to_compare = c("Yes", "No"),
                               pmd_var = "est2_lalonde",
                               group_vars = NULL) {

  # --- Argument Validation ---
  if (length(borrow_to_compare) != 2 || !is.character(borrow_to_compare)) {
    stop("'borrow_to_compare' must be a character vector of length two.")
  }
  if (!"borrow" %in% names(post_inference_all)) {
    stop("'borrow' column not found in the input data frame.")
  }

  # --- Automatically determine group_vars if not provided ---
  if (is.null(group_vars)) {
    all_names <- names(post_inference_all)
    group_cols <- 1 : ( which(names(post_inference_all) == "nsim") - 1)
    group_vars <- all_names[group_cols]
    group_vars <- group_vars[group_vars != "borrow"]
    group_vars <- group_vars[!grepl("^treatment|^true_value", group_vars)]

    if (length(group_vars) == 0) {
      stop("Could not automatically determine grouping variables. Please specify `group_vars` manually.")
    }
  }

  if ("No" %in% borrow_to_compare) {
    # --- Case 1: Comparing TO the "No" borrowing baseline ---
    baseline_level <- "No"
    compare_level <- setdiff(borrow_to_compare, "No")

    baseline_data <- post_inference_all %>%
      filter(borrow == baseline_level) %>%
      select(nsim, all_of(group_vars), "baseline_est" = all_of(pmd_var)) %>%
      distinct()

    compare_data <- post_inference_all %>%
      filter(borrow == compare_level) %>%
      select(nsim, all_of(group_vars), "compare_est" = all_of(pmd_var)) %>%
      distinct()

    if (nrow(baseline_data) == 0) stop("No baseline (No borrowing) results found.")
    if (nrow(compare_data) == 0) stop(paste0("No results found for comparison level: '", compare_level, "'."))

    matching_keys <- c("nsim", group_vars)[
      apply(baseline_data[, c("nsim", group_vars)], 2, function(x) length(unique(x)) != 1)
    ]

    baseline_data_renamed <- baseline_data %>%
      rename_with(~ paste0(.x, "_baseline"), .cols = setdiff(names(baseline_data), c(matching_keys, "baseline_est")))

    pmd_data <- compare_data %>%
      left_join(baseline_data_renamed, by = matching_keys) %>%
      select(nsim, all_of(group_vars), compare_est, baseline_est) %>%
      mutate(pmd = compare_est - baseline_est)

    # varying_cols <- c("nsim", group_vars)[
    #   apply(baseline_data[, c("nsim", group_vars)], 2, function(x) length(unique(x)) != 1)
    # ]
    #
    # for (i in 1: nrow(baseline_data)) {
    #   baseline <- baseline_data[i, ]
    #   match_vals <- baseline[, vary_cols, drop = FALSE]
    #
    #   matches <- purrr::reduce(
    #     vary_cols,
    #     .init = rep(TRUE, nrow(compare_data)),
    #     .f = function(cond, col) cond & (compare_data[[col]] == match_vals[[col]])
    #   )
    #
    #   ## select rows of compare_data that have the values in vary_cols matching the values in baseline, and add baseline_est to this subset baseline["baseline_est"]
    #   compare_data[matches, "baseline_est"] <- baseline$baseline_est
    # }

  } else {
    # --- Case 2: Comparing two different borrowing scenarios ---
    #   comp1_data <- post_inference_all %>%
    #     filter(borrow == borrow_to_compare[1]) %>%
    #     select(all_of(c(group_vars, "nsim", pmd_var)))
    #
    #   comp2_data <- post_inference_all %>%
    #     filter(borrow == borrow_to_compare[2]) %>%
    #     select(all_of(c(group_vars, "nsim", pmd_var)))
    #
    #   if (nrow(comp1_data) == 0 || nrow(comp2_data) == 0) {
    #     stop("One or both of the specified 'borrow_to_compare' levels were not found.")
    #   }
    #
    #   # Use inner_join to only compare scenarios where BOTH levels were run
    #   pmd_data <- inner_join(
    #     comp1_data,
    #     comp2_data,
    #     by = c(group_vars[!group_vars %in% "borrow"], "nsim"),
    #     suffix = c(".c1", ".c2")
    #   ) %>%
    #     mutate(pmd = .data[[paste0(pmd_var, ".c1")]] - .data[[paste0(pmd_var, ".c2")]])
    # }

    # --- Identify Complete Groups for PMD Calculation ---
    pmd_data_raw <- post_inference_all %>%
      filter(borrow %in% borrow_to_compare) %>%
      select(all_of(group_vars), nsim, borrow, all_of(pmd_var))

    complete_groups <- pmd_data_raw %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(n_borrow_levels = n_distinct(borrow), .groups = 'drop') %>%
      filter(n_borrow_levels == 2) %>%
      select(-n_borrow_levels)

    if (nrow(complete_groups) < n_distinct(select(pmd_data_raw, all_of(group_vars)))) {
      warning("Some setting groups were excluded from the PMD calculation because they did not have results for both specified borrowing conditions.")
    }

    if (nrow(complete_groups) == 0) {
      warning("No complete groups found for PMD calculation. Returning an empty data frame.")
      return(data.frame())
    }

    # Use inner_join to only compare scenarios where BOTH levels were run
    pmd_data <- inner_join(
      comp1_data,
      comp2_data,
      by = c(group_vars[!group_vars %in% "borrow"], "nsim"),
      suffix = c(".c1", ".c2")
    ) %>%
      mutate(pmd = .data[[paste0(pmd_var, ".c1")]] - .data[[paste0(pmd_var, ".c2")]])
  }

  # --- Summarize the PMD ---
  pmd_summary <- pmd_data %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(
      mean_pmd = mean(pmd, na.rm = TRUE),
      sd_pmd = sd(pmd, na.rm = TRUE),
      .groups = "drop"
    )

  return(pmd_summary)
}



#' Create a Data Frame of Operating Characteristic Risks
#'
#' @description
#' This function takes the operating characteristics from a simulation run and
#' creates a data frame summarizing the False Go Rate (FGR) and False
#' Stop Rate (FSR) for each simulation scenario.
#'
#' @param oc_all A data frame from `process_sim_results()`.
#' @param lrv The Lower Reference Value used in the simulation.
#' @param tv The Target Value used in the simulation.
#' @param group_vars A character vector of variable names to group by. If `NULL` (default),
#'   the function will automatically use all setting columns.
#'
#' @return A data frame with the FGR and FSR for each simulation scenario.
#'
#' @importFrom dplyr %>% group_by summarize filter select mutate case_when
#' @importFrom tidyr pivot_wider
#' @export
create_risk_df <- function(oc_all, lrv, tv, group_vars = NULL) {
  # --- Argument Validation ---
  if (missing(lrv) || missing(tv)) {
    stop("Arguments 'lrv' and 'tv' must be provided to calculate risks correctly.")
  }

  # --- Automatically determine group_vars if not provided ---
  if (is.null(group_vars)) {
    all_names <- names(oc_all)
    result_patterns <- "^proportion_|^decision_pr|^count_"
    group_vars <- all_names[!grepl(result_patterns, all_names)]
    exclude_patterns <- "^true_value"
    group_vars <- group_vars[!grepl(exclude_patterns, group_vars) & group_vars != "treatment.p"]

    if (length(group_vars) == 0) {
      stop("Could not automatically determine grouping variables. Please specify `group_vars` manually.")
    }
  }

  # --- Calculate Risks (FGR and FSR) ---
  risk_data <- oc_all %>%
    filter(
      (true_value.compare_true == lrv & decision_pr == "Go") |
        (true_value.compare_true == tv & decision_pr == "No-Go")
    ) %>%
    mutate(
      risk_type = case_when(
        decision_pr == "Go" ~ "FGR",
        decision_pr == "No-Go" ~ "FSR",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(across(all_of(c(group_vars, "risk_type")))) %>%
    summarize(proportion = mean(proportion_pr, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = risk_type, values_from = proportion)

  return(risk_data)
}
