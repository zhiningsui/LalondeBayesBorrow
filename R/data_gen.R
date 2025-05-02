#' Generate g-scores
#'
#' This function generates g-scores for multiple arms of a study across multiple simulation repetitions.
#' It simulates individual-level g-scores for each arm based on specified probabilities, means,
#' and standard deviations for the g-score distributions. The function outputs true median values
#' and adjusted median estimates along with their variance.
#'
#' @param n_arms A scalar. Number of arms to be generated.
#' @param nsim A scalar. Number of repetitions for the simulation.
#' @param n_list A named vector. Each element represents the number of patients in each arm.
#' @param prob_list A named vector. Each element represents the probability of a g-score being zero in each arm.
#' @param mu_list A named vector. Each element represents the mean of the normal distribution for generating the log of positive g-scores in each arm.
#' @param sigma_list A named vector. Each element represents the standard deviation of the normal distribution for generating the log of positive g-scores in each arm.
#' @param arm_names A named vector. Each element represents the name of each arm.
#'
#' @return A list of data frames:
#'   - \code{data}: A data frame containing individual-level g-scores for all arms and all repetitions of the simulation.
#'   - \code{true_value}: A data frame with the true medians for all arms and the true median ratio for treatment vs. control (if both are present).
#'   - \code{median_est}: A data frame with adjusted median estimates for each arm, along with the variance estimators for the median estimates.
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply detectCores
#' @importFrom data.table data.table rbindlist setcolorder
#' @importFrom rlang := .data
data_gen_gscore <- function(n_arms = 2, nsim = 10, n_list = NULL, prob_list = NULL, mu_list = NULL, sigma_list = NULL, arm_names = NULL) {
  true_value <- data.frame(arm = arm_names,
                           true_value = mapply(get_true_gscore_median,
                                               prob_list, mu_list, sigma_list)) %>%
    pivot_wider(names_from = .data$arm, values_from = .data$true_value, names_glue = "{arm}.true_value")

  if (sum(arm_names %in% c("treatment", "control")) >= 2) {
    true_value <- true_value %>%
      mutate(compare_true = !!sym(paste0(arm_names[arm_names == "treatment"], ".true_value")) / !!sym(paste0(arm_names[arm_names == "control"], ".true_value")))
  } else {
    true_value <- true_value %>%
      mutate(compare_true = NA)
  }

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncore <- 2
  } else {
    # use all cores in devtools::test()
    ncore <- detectCores() - 1
  }
  cl <- makeCluster(ncore)
  required_vars <- c("n_list", "mu_list", "sigma_list", "prob_list", "arm_names",
                     "simulate_gscore", "sd_gscore_median")

  clusterExport(cl, required_vars, envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(data.table)
  })

  results <- parLapply(cl, 1:nsim, function(nrep) {
    data_temp <- data.table::data.table(
      nsim = rep(nrep, each = sum(n_list)),
      id = unlist(lapply(n_list, seq)),
      arm = unlist(as.list(mapply(rep, arm_names, n_list))),
      g_score = unlist(as.list(mapply(simulate_gscore, n_list, mu_list, sigma_list, prob_list)))
    )

    data_temp$arm <- factor(data_temp$arm, levels = arm_names)
    rownames(data_temp) <- NULL

    med_temp <- data_temp %>%
      group_by(.data$arm) %>%
      summarize(
        median_unadjusted = stats::median(.data$g_score),
        median_sd = sd_gscore_median(.data$g_score)
      ) %>%
      mutate(
        n = n_list,
        median_adjusted = ifelse(.data$median_unadjusted == 0, min(data_temp$g_score[data_temp$g_score > 0]) / sqrt(n), .data$median_unadjusted),
        nsim = nrep
      )

    list(data = data_temp, median = med_temp)
  })

  stopCluster(cl)

  data <- rbindlist(lapply(results, "[[", "data"))
  g_median <- rbindlist(lapply(results, "[[", "median"))
  setcolorder(g_median, c("nsim", "arm", "n", "median_unadjusted", "median_adjusted", "median_sd"))
  return(list(data = data, true_value = true_value, median_est = g_median))
}



#' Generate Binary Response Data
#'
#' This function simulates binary endpoint data (0/1) for multiple study arms
#' across several simulation repetitions. Each arm's data is generated based on
#' the specified probabilities of response using a binomial distribution. The
#' function returns individual-level data, true response rates, and response
#' frequencies per arm per simulation.
#'
#' @param nsim A scalar. The number of repetitions of the simulation to be performed.
#' @param n_list A named list. Each element represents the number of patients in
#'   the corresponding arm (names must match arm_names). Must contain positive integers.
#' @param prob_list A named list. Each element represents the probability of response (binary outcome = 1)
#'   for the corresponding arm (names must match arm_names). Must contain numeric values between 0 and 1.
#' @param arm_names A character vector of arm names.
#' @param n_arms Number of arms (optional, inferred from arm_names if not provided).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: A data.table with columns nsim, id, arm, resp, containing
#'     all simulated binary outcomes (0/1) for all arms and all repetitions.
#'   \item \code{true_value}: A data frame of true response rates for each arm,
#'     along with the difference in response rates between the treatment and control
#'     arms (if both are present).
#'   \item \code{freq}: A data.table of response frequencies per arm per simulation,
#'     with columns for nsim, arm, count (number of responses), and n (total patients).
#' }
#'
#' @importFrom dplyr %>% group_by summarise mutate ungroup select
#' @importFrom tidyr pivot_wider
#' @importFrom data.table data.table rbindlist setcolorder
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply detectCores
#' @importFrom stats rbinom
#' @importFrom rlang := .data !! sym
#'
#' @export
data_gen_binary <- function(n_arms = 2, nsim = 10, n_list = NULL, prob_list = NULL, arm_names = NULL) {

  # --- Input Validation ---
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 ) {
    stop("nsim must be a single positive integer.")
  }
  if (!is.character(arm_names) || length(arm_names) == 0) {
    stop("arm_names must be a character vector with at least one arm name.")
  }

  expected_list_names <- sort(arm_names)
  if (!all(sort(names(n_list)) == expected_list_names) ||
      !all(sort(names(prob_list)) == expected_list_names)) {
    stop("Names of n_list and prob_list must match arm_names.")
  }

  # Validate contents of the lists
  for (arm in arm_names) {
    if (!is.numeric(n_list[[arm]]) || length(n_list[[arm]]) != 1 || n_list[[arm]] <= 0 ) {
      stop(paste("Invalid sample size for arm", arm, ": must be a single positive integer."))
    }
    if (!is.numeric(prob_list[[arm]]) || length(prob_list[[arm]]) != 1 || prob_list[[arm]] < 0 || prob_list[[arm]] > 1) {
      stop(paste("Invalid probability for arm", arm, ": must be a single numeric value between 0 and 1."))
    }
  }

  # --- True rate values and comparison ---
  rate_true_df <- data.frame(arm = arm_names,
                             rate_true = unlist(prob_list)) %>%
    pivot_wider(names_from = .data$arm, values_from = .data$rate_true,
                names_glue = "{arm}.rate_true")

  # Calculate true difference if both 'treatment' and 'control' arms exist
  if ("treatment" %in% arm_names && "control" %in% arm_names) {
    rate_true <- rate_true_df %>%
      mutate(compare_true = !!sym("treatment.rate_true") - !!sym("control.rate_true"))
  } else {
    rate_true <- rate_true_df %>%
      mutate(compare_true = NA) # Comparison not applicable
  }

  # --- Parallel Processing Setup ---
  # Determine number of cores to use
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncore <- 2
  } else {
    # use all cores available minus 1, or at least 1
    ncore <- parallel::detectCores() - 1
  }
  if (ncore < 1) ncore <- 1

  cl <- makeCluster(ncore)
  # Ensure cluster is stopped when function exits (either normally or due to error)
  on.exit(stopCluster(cl))

  # Export necessary objects and load libraries to cluster environment
  clusterExport(cl, list("n_list", "prob_list", "arm_names", "n_arms"),
                envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(data.table)
    library(stats) # For rbinom
    library(rlang) # For !! sym, .data
  })

  # --- Parallel Simulation Execution ---
  results <- parLapply(cl, 1:nsim, function(nrep) {
    # Using mapply to apply rbinom with arm-specific parameters
    simulated_values <- unlist(as.list(
      mapply(stats::rbinom,
             n = n_list,
             size = as.list(rep(1, n_arms)),
             prob = prob_list,
             SIMPLIFY = FALSE)))

    # Create data.table for the current replicate
    data_temp <- data.table::data.table(
      nsim = rep(nrep, each = sum(unlist(n_list))),
      id = unlist(lapply(n_list, seq_len)), # Use seq_len for safety
      arm = unlist(mapply(rep, names(n_list), unlist(n_list), SIMPLIFY = FALSE)), # Use names(n_list) for arm names
      resp = simulated_values
    )

    # Calculate summary statistics for the current replicate
    summary_stats <- data_temp %>%
      group_by(arm) %>%
      summarise(
        nsim = unique(nsim),
        n = n(),
        count = sum(resp == 1),
        .groups = "drop"
      )

    list(data = data_temp, summary = summary_stats)
  })

  # --- Combine Results from all Simulations ---
  data_all <- data.table::rbindlist(lapply(results, "[[", "data"))
  summary_all <- data.table::rbindlist(lapply(results, "[[", "summary")) %>%
    pivot_wider(names_from = arm, values_from = c(n, count),
                names_glue = "{arm}.{.value}")

  # Return the combined results and true values
  return(list(
    data = data_all,
    true_value = rate_true,
    summary = summary_all
  ))
}

#' Generate continuous outcome data for multiple arms across multiple simulations
#'
#' Generates continuous outcome data by sampling from Normal distributions for each arm
#' in each simulation replicate.
#'
#' @param nsim Number of simulations
#' @param n_list Named list of sample sizes per arm (names should match arm_names). Must contain positive integers.
#' @param mu_list Named list of means per arm (names should match arm_names). Must contain numeric values.
#' @param sigma_list Named list of standard deviations per arm (names should match arm_names). Must contain non-negative numeric values.
#' @param arm_names Character vector of arm names.
#' @param n_arms Number of arms (optional, inferred from arm_names if not provided).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: A data.table with columns nsim, id, arm, value, containing all simulated continuous values.
#'   \item \code{summary}: A data.table with columns nsim and arm-specific summary statistics (n, mu_hat, sd, s).
#'   \item \code{true_value}: A data frame with true mean values per arm and potentially a true mean difference (treatment - control).
#' }
#' @importFrom dplyr %>% group_by summarise mutate
#' @importFrom tidyr pivot_wider
#' @importFrom data.table data.table rbindlist
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply detectCores
#' @importFrom stats rnorm sd
#' @importFrom rlang .data !! sym
#' @export
data_gen_continuous <- function(nsim = 10, n_list, mu_list, sigma_list, arm_names, n_arms = length(arm_names)) {

  # --- Input Validation ---
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0) {
    stop("nsim must be a single positive integer.")
  }
  if (!is.character(arm_names) || length(arm_names) == 0) {
    stop("arm_names must be a character vector with at least one arm name.")
  }

  expected_list_names <- sort(arm_names)
  if (!all(sort(names(n_list)) == expected_list_names) ||
      !all(sort(names(mu_list)) == expected_list_names) ||
      !all(sort(names(sigma_list)) == expected_list_names)) {
    stop("Names of n_list, mu_list, and sigma_list must match arm_names.")
  }

  # Validate contents of the lists
  for (arm in arm_names) {
    if (!is.numeric(n_list[[arm]]) || length(n_list[[arm]]) != 1 || n_list[[arm]] <= 0 ) {
      stop(paste("Invalid sample size for arm", arm, ": must be a single positive integer."))
    }
    if (!is.numeric(mu_list[[arm]]) || length(mu_list[[arm]]) != 1) {
      stop(paste("Invalid mean for arm", arm, ": must be a single numeric value."))
    }
    if (!is.numeric(sigma_list[[arm]]) || length(sigma_list[[arm]]) != 1 || sigma_list[[arm]] < 0) {
      stop(paste("Invalid standard deviation for arm", arm, ": must be a single non-negative numeric value."))
    }
  }

  # --- True mean values and comparison ---
  mean_true_df <- data.frame(arm = arm_names,
                             mean_true = unlist(mu_list)) %>%
    pivot_wider(names_from = .data$arm, values_from = .data$mean_true,
                names_glue = "{arm}.mean_true")

  # Calculate true difference if both 'treatment' and 'control' arms exist
  if ("treatment" %in% arm_names && "control" %in% arm_names) {
    mean_true <- mean_true_df %>%
      mutate(compare_true = !!sym(paste0("treatment", ".mean_true")) - !!sym(paste0("control", ".mean_true")))
  } else {
    mean_true <- mean_true_df %>%
      mutate(compare_true = NA) # Comparison not applicable
  }


  # --- Parallel Processing Setup ---
  # Determine number of cores to use
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncore <- 2
  } else {
    # use all cores available minus 1, or at least 1
    ncore <- parallel::detectCores() - 1
  }
  if (ncore < 1) ncore <- 1

  cl <- makeCluster(ncore)
  # Ensure cluster is stopped when function exits (either normally or due to error)
  on.exit(stopCluster(cl))

  # Export necessary objects and load libraries to cluster environment
  clusterExport(cl, list("n_list", "mu_list", "sigma_list", "arm_names"),
                envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(data.table)
    library(stats) # For rnorm, sd
    library(rlang) # For !! sym
  })

  # --- Parallel Simulation Execution ---
  results <- parLapply(cl, 1:nsim, function(nrep) {

    # Using mapply to apply rnorm with arm-specific parameters
    simulated_values <- unlist(mapply(stats::rnorm,
                                      n = n_list,
                                      mean = mu_list,
                                      sd = sigma_list,
                                      SIMPLIFY = FALSE))

    # Create data.table for the current replicate
    data_temp <- data.table::data.table(
      nsim = rep(nrep, each = sum(unlist(n_list))),
      id = unlist(lapply(n_list, seq_len)), # Use seq_len for safety
      arm = unlist(mapply(rep, names(n_list), unlist(n_list), SIMPLIFY = FALSE)), # Use names(n_list) for arm names
      value = simulated_values
    )

    # Calculate summary statistics for the current replicate
    summary_stats <- data_temp %>%
      group_by(arm) %>%
      summarise(
        nsim = unique(nsim),
        n = n(),
        mu_hat = mean(value),
        sd = stats::sd(value),
        .groups = "drop"
      ) %>%
      mutate(s = .data$sd / sqrt(.data$n)) # Standard error of the mean

    list(data = data_temp, summary = summary_stats)
  })

  stopCluster(cl)

  # --- Combine Results from all Simulations ---
  data_all <- data.table::rbindlist(lapply(results, "[[", "data"))
  summary_all <- data.table::rbindlist(lapply(results, "[[", "summary")) %>%
    pivot_wider(names_from = arm, values_from = c(n, mu_hat, sd, s),
                names_glue = "{arm}.{.value}")

  return(list(
    data = data_all,
    summary = summary_all,
    true_value = mean_true
  ))
}


#' Simulate Constrained Depth of Response (DpR) values for a single arm
#'
#' Generates DpR values based on a two-part mixture model with a point mass at -1
#' and a shifted Gamma distribution for the continuous part.
#'
#' @param n Sample size for the arm
#' @param p Proportion of -1 values (probability of complete response)
#' @param mu Target mean for the arm
#' @param sigma Target standard deviation for the arm
#'
#' @return A numeric vector of simulated DpR values for the arm.
#' @importFrom stats rbinom rgamma
#' @keywords internal
simulate_DpR <- function(n, p, mu, sigma) {

  # --- Parameter Validation for the current arm ---
  if (!is.numeric(n) || length(n) != 1 || n <= 0 ) {
    stop("Invalid sample size: must be a single positive integer.")
  }
  if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1) {
    stop("Invalid p value: must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(mu) || length(mu) != 1 || mu < -1) {
    stop("Invalid mean: must be a single numeric value greater than or equal to -1.")
  }
  # Sigma can be 0 if p is 1 and mu is -1
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma < 0) {
    stop("Invalid standard deviation: must be a single non-negative numeric value.")
  }

  sigma2 <- sigma^2

  # --- Simulate DpR values based on the mixture model ---
  if (p == 1) { # If p is 1
    # All values are -1. Mean must be -1 and variance 0.
    if (mu != -1.0 || sigma2 != 0.0) {
      stop("Invalid parameters: When p is 1, mu must be -1 and sigma2 must be 0.")
    }
    # Simulate n values of -1
    simulated_dpr <- rep(-1, n)

  } else { # If p < 1
    # Check the variance validity condition for p < 1
    required_sigma2_lower_bound <- (p / (1 - p)) * (mu + 1)^2
    if (sigma2 < required_sigma2_lower_bound - 1e-9) { # Using a small tolerance for floating point comparison for the strict inequality
      stop(paste0("Invalid variance: sigma2 (", sigma2, ") must be greater than or equal to ", required_sigma2_lower_bound, " for p < 1."))
    }

    # Calculate the denominator for alpha and beta formulas
    denominator <- (1 - p) * sigma2 - p * (mu + 1)^2
    if (denominator <= 1e-9) { # Using a small tolerance to check denominator is positive (should be if sigma2 is valid and p < 1)
      stop("Denominator for alpha/beta calculation is not sufficiently positive.")
    }

    # Calculate alpha (shape) and beta (rate)
    alpha <- (mu + 1)^2 / denominator
    beta <- (1 - p) * (mu + 1) / denominator

    # Ensure alpha and beta are positive (required for Gamma distribution)
    if (alpha <= 0 || beta <= 0) {
      stop(paste("Calculated alpha (", alpha, ") or beta (", beta, ") is not positive, which is required for a Gamma distribution."))
    }

    # Simulate DpR values for p < 1 using the mixture model
    is_complete_response <- rbinom(n, 1, p)
    gamma_values <- rgamma(n, shape = alpha, rate = beta)

    simulated_dpr <- numeric(n)
    simulated_dpr[is_complete_response == 1] <- -1
    simulated_dpr[is_complete_response == 0] <- gamma_values[is_complete_response == 0] - 1
  }

  return(simulated_dpr)
}

#' Generate Constrained Depth of Response (DpR) data for multiple arms across multiple simulations
#'
#' @param n_arms Number of arms (can be inferred from arm_names, but kept for consistency with template)
#' @param nsim Number of simulations
#' @param n_list Named list of sample sizes per arm (names should match arm_names)
#' @param p_list Named list of proportion of -1 per arm (names should match arm_names)
#' @param mu_list Named list of means per arm (theta_x) (names should match arm_names)
#' @param sigma_list Named list of standard deviations per arm (sigma_x) (names should match arm_names)
#' @param arm_names Character vector of arm names
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: A data.table with columns nsim, id, arm, value, containing all simulated DpR values.
#'   \item \code{summary}: A data.table with columns nsim and arm-specific summary statistics (n, mu_hat, sd, s).
#'   \item \code{true_value}: A data frame with true mean differences (currently only 'treatment' vs 'control').
#' }
#' @importFrom dplyr %>% group_by summarise mutate
#' @importFrom tidyr pivot_wider
#' @importFrom data.table data.table rbindlist
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply
#' @importFrom stats rbinom rgamma sd
#' @importFrom rlang .data !! sym
#' @export
data_gen_DpR <- function(n_arms, nsim, n_list, p_list, arm_names, mu_list, sigma_list) {

  # --- Input Validation ---
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0) {
    stop("nsim must be a single positive integer.")
  }
  if (!is.character(arm_names) || length(arm_names) == 0) {
    stop("arm_names must be a character vector with at least one arm name.")
  }
  if (!all(sort(names(n_list)) == sort(arm_names)) ||
      !all(sort(names(p_list)) == sort(arm_names)) ||
      !all(sort(names(mu_list)) == sort(arm_names)) ||
      !all(sort(names(sigma_list)) == sort(arm_names))) {
    stop("Names of n_list, p_list, mu_list, and sigma_list must match arm_names.")
  }

  # True mean values and comparison (adapted from continuous version)
  mean_true <- data.frame(arm = arm_names,
                          mean_true = unlist(mu_list)) %>%
    pivot_wider(names_from = .data$arm, values_from = .data$mean_true,
                names_glue = "{arm}.mean_true")

  # Check if both 'treatment' and 'control' arms exist for comparison
  if ("treatment" %in% arm_names && "control" %in% arm_names) {
    mean_true <- mean_true %>%
      mutate(compare_true = !!sym(paste0("treatment", ".mean_true")) - !!sym(paste0("control", ".mean_true")))
  } else {
    mean_true <- mean_true %>%
      mutate(compare_true = NA) # Comparison not applicable
  }

  # --- Parallel Processing Setup ---
  # Determine number of cores to use
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncore <- 2
  } else {
    # use all cores in devtools::test() - 1
    ncore <- parallel::detectCores() - 1
  }
  if (ncore < 1) ncore <- 1 # Ensure at least one core

  cl <- makeCluster(ncore)
  on.exit(stopCluster(cl)) # Ensure cluster is stopped on exit

  # Export necessary objects and load libraries to cluster
  clusterExport(cl, list("n_list", "p_list", "mu_list", "sigma_list", "arm_names", "simulate_DpR"),
                envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(data.table)
    library(stats) # For rbinom, rgamma, sd
    library(rlang) # For !! sym (if needed within cluster)
  })

  # --- Parallel Simulation Execution ---
  results <- parLapply(cl, 1:nsim, function(nrep) {

    data_list_replicate <- list()

    for (arm_name in arm_names) {
      n_arm <- n_list[[arm_name]]
      p_arm <- p_list[[arm_name]]
      mu_arm <- mu_list[[arm_name]]
      sigma_arm <- sigma_list[[arm_name]]

      simulated_dpr <- tryCatch({
        simulate_DpR(n = n_arm, p = p_arm, mu = mu_arm, sigma = sigma_arm)
      }, error = function(e) {
        stop(paste("Error simulating DpR for arm", arm_name, "in simulation replicate", nrep, ":", e$message))
      })

      # Store the simulated values for the current arm in the current simulation replicate's list
      data_list_replicate[[arm_name]] <- data.table(
        nsim = nrep,
        id = 1:n_arm,
        arm = arm_name,
        value = simulated_dpr
      )
    }
    data_temp <- data.table::rbindlist(data_list_replicate)

    # Calculate summary statistics for the current replicate
    summary_stats <- data_temp %>%
      group_by(arm) %>%
      summarise(
        nsim = unique(nsim),
        n = n(),
        mu_hat = mean(value),
        sd = sd(value),
        .groups = "drop"
      ) %>%
      mutate(s = sd / sqrt(n)) # Standard error of the mean

    list(data = data_temp, summary = summary_stats)
  })

  # --- Combine Results from all Simulations ---
  data_all <- data.table::rbindlist(lapply(results, "[[", "data"))
  summary_all <- data.table::rbindlist(lapply(results, "[[", "summary")) %>%
    pivot_wider(names_from = arm, values_from = c(n, mu_hat, sd, s),
                names_glue = "{arm}.{.value}")

  # Return the combined results and true values
  return(list(
    data = data_all,
    summary = summary_all,
    true_value = mean_true
  ))
}
