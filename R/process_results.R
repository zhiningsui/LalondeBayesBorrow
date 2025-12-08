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

  res_names <- names(bayes_results)
  # Derive design_scenario and prior_name from the names
  meta_df <- tibble(
    name = res_names,
    design_scenario = str_extract(res_names, "(?<=design)\\d+"),
    prior_name = str_extract(res_names, "(?<=__)\\S+")
  )

  # Helper function to bind rows with new metadata columns
  bind_with_meta <- function(component) {
    dplyr::bind_rows(
      lapply(seq_along(bayes_results), function(i) {
        df <- bayes_results[[i]][[component]]
        if (is.null(df)) return(NULL)
        df %>%
          mutate(
            design_scenario = meta_df$design_scenario[i],
            prior_name = meta_df$prior_name[i]
          )
      })
    )
  }

  final_results <- list(
    post_params_all = bind_with_meta("post_params"),
    post_est_ci_all = bind_with_meta("post_est_ci"),
    post_inference_all = bind_with_meta("post_inference")
  )


  # final_results <- list(
  #   post_params_all = dplyr::bind_rows(lapply(bayes_results, `[[`, "post_params")),
  #   post_est_ci_all = dplyr::bind_rows(lapply(bayes_results, `[[`, "post_est_ci")),
  #   post_inference_all = dplyr::bind_rows(lapply(bayes_results, `[[`, "post_inference"))
  # )

  # scn <- final_results$post_params_all %>%
  #   distinct(pick(all_of(scn_var))) %>%
  #   mutate(dat_gen_scenario = 1:n())

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

#' Create a Posterior Mean Difference (PMD) summary
#'
#' @description
#' Compute the Posterior Mean Difference (PMD) between two borrowing
#' conditions (e.g., `"Yes"` vs `"No"`) across simulations and user-specified
#' grouping variables. Works for one-to-one comparisons and the common
#' "vs. No borrowing" baseline case by the same unified logic.
#'
#' @details
#' The function filters `post_inference_all` to the two levels in
#' `borrow_to_compare`, pivots wide on `borrow`, and computes
#' \eqn{\text{PMD} = \hat{\theta}_{\text{level1}} - \hat{\theta}_{\text{level2}}}
#' within each simulation (`nsim`) and group. It then summarizes
#' `mean_pmd` and `sd_pmd` across simulations for each group.
#'
#' If `group_vars = NULL`, the function will:
#' 1) use any present columns from `scn_var` and `prior_var`;
#' 2) optionally construct a compact control-only scenario ID
#'    (`design_scenario.control`) when `control.only = TRUE` and there are
#'    multiple distinct "control" columns to define the grouping.
#'
#' Returned metadata tibbles `priors` and `scenarios` are included **only** when
#' `group_vars = NULL` and the corresponding inputs exist.
#'
#' @param post_inference_all A data frame of posterior results (one row per simulation),
#'   typically from `bayesian_lalonde_decision(...)[["post_inference"]]`.
#'   Must contain columns `nsim`, `borrow`, and the PMD variable (see `pmd_var`).
#' @param borrow_to_compare Character(2). The two levels of `borrow` to contrast
#'   (difference computed as `borrow_to_compare[1] - borrow_to_compare[2]`).
#' @param scn_var Character vector of scenario column names to include if present
#'   (default `"design_scenario"`).
#' @param prior_var Character vector of prior/borrrowing column names to include
#'   if present (default `"prior_name"`).
#' @param pmd_var Single column name with the posterior estimate used for PMD
#'   (default `"est2_lalonde"`).
#' @param group_vars Character vector of grouping variables. If `NULL`, they are
#'   auto-detected from `scn_var`/`prior_var` (and possibly a compact control ID).
#' @param control.only Logical. If `TRUE`, and control-related columns exist,
#'   create a compact `design_scenario.control` grouping ID (default `TRUE`).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pmd_summary}: tibble with one row per group and columns
#'         `mean_pmd`, `sd_pmd`.
#'   \item \code{priors}: optional tibble of prior-related grouping columns
#'         (returned only when `group_vars = NULL` and `prior_var` found).
#'   \item \code{scenarios}: optional tibble of scenario-related grouping columns
#'         (returned only when `group_vars = NULL` and `scn_var` found).
#' }
#'
#' @importFrom dplyr %>% group_by summarise filter select mutate inner_join left_join distinct n_distinct across all_of relocate everything
#' @importFrom tidyr pivot_wider
#' @export
create_pmd_summary <- function(post_inference_all,
                               borrow_to_compare = c("Yes", "No"),
                               scn_var = "design_scenario",
                               prior_var = "prior_name",
                               pmd_var = "est2_lalonde",
                               group_vars = NULL,
                               control.only = TRUE) {

  # --- Argument Validation ---
  if (!is.data.frame(post_inference_all)) {
    stop("'post_inference_all' must be a data frame.")
  }
  if (length(borrow_to_compare) != 2 || !is.character(borrow_to_compare)) {
    stop("'borrow_to_compare' must be a character vector of length two.")
  }
  if (!"borrow" %in% names(post_inference_all)) {
    stop("'borrow' column not found in the input data frame.")
  }

  # --- Automatically determine group_vars if not provided ---
  if (is.null(group_vars)) {
    cols_group <- intersect(scn_var, names(post_inference_all))
    if (length(cols_group) == 0) stop("None of scn_var columns exist in post_inference_all.")
    if (!setequal(scn_var, cols_group)) {
      warning("Some scn_var names are missing in post_inference_all: ",
              paste(setdiff(scn_var, cols_group), collapse = ", "))
    }

    constancy_check <- post_inference_all %>%
      group_by(across(all_of(scn_var))) %>%
      summarise(across(everything(), n_distinct), .groups = "drop")
    # Keep columns that are constant (n_distinct == 1) across all scenarios
    constant_cols <- constancy_check %>%
      select(where(~ all(.x == 1))) %>%
      names()

    scn_params_all <- post_inference_all %>%
      distinct(across(all_of(c(scn_var, constant_cols)))) %>%
      select(!where(~ n_distinct(.) == 1))

    scn_vars_all <- colnames(scn_params_all)

    if (control.only) {
      scn_vars <- scn_vars_all[grepl("control", scn_vars_all)]
      scn_control <- post_inference_all %>%
        select(all_of(scn_vars)) %>%
        distinct() %>%
        mutate(design_scenario.control = 1:n())

      post_inference_all <- post_inference_all %>%
        left_join(scn_control)

      scn_vars <- c(scn_var, "design_scenario.control", scn_vars)
      scn_params <- scn_control

      scn_vars_all <- unique(c(scn_vars_all, scn_vars))
      scn_params_all <- post_inference_all %>%
        distinct(across(all_of(scn_vars_all))) %>%
        select(!where(~ n_distinct(.) == 1))

    } else {
      scn_params <- scn_params_all
    }

    cols_group <- intersect(prior_var, names(post_inference_all))
    if (length(cols_group) == 0) stop("None of prior_var columns exist in post_inference_all.")
    if (!setequal(prior_var, cols_group)) {
      warning("Some prior_var names are missing in post_inference_all: ",
              paste(setdiff(prior_var, cols_group), collapse = ", "))
    }
    constancy_check <- post_inference_all %>%
      group_by(across(all_of(prior_var))) %>%
      summarise(across(everything(), n_distinct), .groups = "drop")

    constant_cols <- constancy_check %>%
      select(where(~ all(.x == 1))) %>%
      names()
    prior_params <- post_inference_all %>%
      distinct(across(all_of(c(prior_var, constant_cols)))) %>%
      select(!where(~ n_distinct(.) == 1))

    prior_vars <- colnames(prior_params)
    group_vars <- c(scn_vars, prior_vars)
  }

  if ("No" %in% borrow_to_compare) {
    # --- Case 1: Comparing TO the "No" borrowing baseline ---
    baseline_level <- "No"
    compare_level <- setdiff(borrow_to_compare, "No")

    baseline_data <- post_inference_all %>%
      filter(borrow == baseline_level) %>%
      select(nsim, all_of(group_vars), "baseline_est" = all_of(pmd_var))

    compare_data <- post_inference_all %>%
      filter(borrow == compare_level) %>%
      select(nsim, all_of(group_vars), "compare_est" = all_of(pmd_var))

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

  } else {
    # --- Case 2: Comparing two different borrowing scenarios ---
    comp1_data <- post_inference_all %>%
      filter(borrow == borrow_to_compare[1]) %>%
      select(all_of(c(group_vars, "nsim", pmd_var)))

    comp2_data <- post_inference_all %>%
      filter(borrow == borrow_to_compare[2]) %>%
      select(all_of(c(group_vars, "nsim", pmd_var)))

    if (nrow(comp1_data) == 0 || nrow(comp2_data) == 0) {
      stop("One or both of the specified 'borrow_to_compare' levels were not found.")
    }

    # --- Identify Complete Groups for PMD Calculation ---
    pmd_data_raw <- post_inference_all %>%
      filter(borrow %in% borrow_to_compare) %>%
      select(all_of(group_vars), nsim, borrow, all_of(pmd_var))

    complete_groups <- pmd_data_raw %>%
      group_by(across(all_of(colnames(scn_params)))) %>%
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
  if (control.only) {
    group_vars <- group_vars[group_vars != "design_scenario"]
  }
  pmd_summary <- pmd_data %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(
      mean_pmd = mean(pmd, na.rm = TRUE),
      sd_pmd = sd(pmd, na.rm = TRUE),
      .groups = "drop"
    )

  return(list(
    pmd_summary = pmd_summary,
    priors = prior_params,
    scenarios = scn_params_all)
  )
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
    mutate(
      risk_type = case_when(
        (true_value.compare_true == lrv & decision_pr == "Go") ~ "FGR",
        (true_value.compare_true == tv & decision_pr == "No-Go") ~ "FSR",
        (true_value.compare_true == lrv & decision_pr == "No-Go") ~ "CSR",
        (true_value.compare_true == tv & decision_pr == "Go") ~ "CGR",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(risk_type)) %>%
    # filter(
    #   (true_value.compare_true == lrv & decision_pr == "Go") |
    #     (true_value.compare_true == tv & decision_pr == "No-Go") |
    #     (true_value.compare_true == lrv & decision_pr == "No-Go") |
    #     (true_value.compare_true == tv & decision_pr == "Go")
    # ) %>%
    # mutate(
    #   risk_type = case_when(
    #     decision_pr == "Go" ~ "FGR",
    #     decision_pr == "No-Go" ~ "FSR",
    #     TRUE ~ NA_character_
    #   )
    # ) %>%
    group_by(across(all_of(c(group_vars, "risk_type")))) %>%
    summarize(proportion = mean(proportion_pr, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = risk_type, values_from = proportion)

  return(risk_data)
}
