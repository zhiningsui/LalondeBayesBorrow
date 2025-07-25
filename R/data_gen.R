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
#' @param seed An integer or `NULL`. If an integer, a seed is set for
#'   reproducibility of the random number generation across parallel workers.
#'   Defaults to `NULL`, which means no explicit seed is set beyond R's default behavior.
#'
#' @return A list containing:
#'   * `data`: A data.table with columns `nsim`, `id`, `arm`, `resp`, containing
#'     all simulated binary outcomes (0/1) for all arms and all repetitions.
#'   * `true_value`: A data frame of true response rates for each arm,
#'     along with the difference in response rates between the treatment and control
#'     arms (if both are present).
#'   * `summary`: A data.table of response frequencies per arm per simulation, # Corrected from 'freq' to 'summary' based on code
#'     with columns for `nsim`, `arm`, `count` (number of responses), and `n` (total patients).
#'
#' @importFrom dplyr %>% group_by summarise mutate ungroup select
#' @importFrom tidyr pivot_wider
#' @importFrom data.table data.table rbindlist setcolorder
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply detectCores
#' @importFrom stats rbinom
#' @importFrom rlang := .data !! sym
#'
#' @export
#' @seealso `bayesian_lalonde_decision`
data_gen_binary <- function(n_arms = 2, nsim = 10, n_list = NULL, prob_list = NULL, arm_names = NULL, seed = NULL) {

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

  if (!is.null(seed)) {
    set.seed(seed)
    clusterSetRNGStream(cl, seed)
  } else {
    clusterSetRNGStream(cl)
  }

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
#' in each simulation replicate. It produces both the individual-level data and
#' arm-specific summary statistics (sample size, mean, standard deviation, standard error).
#' The function supports parallel processing for faster execution across simulations
#' and includes options for reproducibility via setting a random seed.
#'
#' @param nsim A single positive integer. Number of simulations to run.
#' @param n_list A named list of positive integers. Sample sizes for each arm.
#'   Names must exactly match `arm_names`.
#' @param mu_list A named list of numeric values. True means for the Normal
#'   distribution for each arm. Names must exactly match `arm_names`.
#' @param sigma_list A named list of non-negative numeric values. True standard
#'   deviations for the Normal distribution for each arm. Names must exactly
#'   match `arm_names`.
#' @param arm_names A character vector. Names of the arms in the simulation. Used
#'   for naming outputs and matching list elements.
#' @param seed An integer or `NULL`. If an integer, a seed is set for
#'   reproducibility of the random number generation across parallel workers.
#'   Defaults to `NULL`, which means no explicit seed is set beyond R's default behavior.
#'
#' @return A list containing:
#'   * `data`: A `data.table` with columns `nsim` (simulation replicate ID),
#'     `id` (patient ID within the replicate and arm), `arm` (arm name),
#'     `value` (simulated continuous outcome). Contains all individual simulated values.
#'   * `summary`: A `data.table` with columns `nsim` and arm-specific summary
#'     statistics calculated from the simulated data for each arm and simulation
#'     replicate. Columns are in wide format, named `[arm_name].n`, `[arm_name].mu_hat`
#'     (sample mean), `[arm_name].sd` (sample standard deviation), and `[arm_name].s`
#'     (sample standard error of the mean, `sd/sqrt(n)`).
#'   * `true_value`: A `data.frame` with the true input mean values per arm
#'     (columns `[arm_name].mean_true`) and a true mean difference between 'treatment'
#'     and 'control' arms if both exist (column `compare_true`).
#'
#' @importFrom dplyr %>% group_by summarise mutate select
#' @importFrom tidyr pivot_wider
#' @importFrom data.table data.table rbindlist
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply detectCores clusterSetRNGStream
#' @importFrom stats rnorm sd
#' @importFrom rlang .data !! sym
#' @export
#'
#' @examples
#' # Example 1: Generate data for 2 arms over 5 simulations
#' n_list <- list(treatment = 50, control = 60)
#' mu_list <- list(treatment = 10, control = 8)
#' sigma_list <- list(treatment = 3, control = 3.5)
#' arm_names <- c("treatment", "control")
#'
#' sim_data <- data_gen_continuous(nsim = 5, n_list, mu_list, sigma_list, arm_names)
#'
#' # Display structure of results
#' str(sim_data)
#'
#' # Display first few rows of individual data
#' print(head(sim_data$data))
#'
#' # Display the summary data
#' print(sim_data$summary)
#'
#' # Display the true values
#' print(sim_data$true_value)
#'
#' # Example 2: Generate data with a seed for reproducibility
#' sim_data_reproducible <- data_gen_continuous(nsim = 3, n_list, mu_list,
#'                                             sigma_list, arm_names, seed = 123)
#' print(sim_data_reproducible$summary) # Note the summary values
#'
#' sim_data_reproducible_again <- data_gen_continuous(nsim = 3, n_list, mu_list,
#'                                                  sigma_list, arm_names, seed = 123)
#' print(sim_data_reproducible_again$summary) # Values should match the previous run
#'
#' @seealso \code{\link{bayesian_lalonde_decision}}
data_gen_continuous <- function(nsim = 10, n_list, mu_list, sigma_list, arm_names, seed = NULL) {

  # --- Input Validation ---
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0) {
    stop("nsim must be a single positive integer.")
  }
  if (!is.character(arm_names) || length(arm_names) == 0) {
    stop("arm_names must be a character vector with at least one arm name.")
  }

  if (!is.list(n_list) || is.null(names(n_list)) || !is.list(mu_list) || is.null(names(mu_list)) || !is.list(sigma_list) || is.null(names(sigma_list))) {
    stop("n_list, mu_list, and sigma_list must be named lists.")
  }

  expected_list_names <- sort(arm_names)
  if (!all(sort(names(n_list)) == expected_list_names) ||
      !all(sort(names(mu_list)) == expected_list_names) ||
      !all(sort(names(sigma_list)) == expected_list_names)) {
    stop("Names of n_list, mu_list, and sigma_list must match arm_names.")
  }


  # Validate contents of the lists
  for (arm in arm_names) {
    if (!is.numeric(n_list[[arm]]) || length(n_list[[arm]]) != 1 || !is.finite(n_list[[arm]]) || n_list[[arm]] <= 0 || n_list[[arm]] != floor(n_list[[arm]]) ) {
      stop(paste("Invalid sample size for arm '", arm, "': must be a single positive integer.", sep=""))
    }
    if (!is.numeric(mu_list[[arm]]) || length(mu_list[[arm]]) != 1 || !is.finite(mu_list[[arm]])) {
      stop(paste("Invalid mean for arm '", arm, "': must be a single finite numeric value.", sep=""))
    }
    if (!is.numeric(sigma_list[[arm]]) || length(sigma_list[[arm]]) != 1 || !is.finite(sigma_list[[arm]]) || sigma_list[[arm]] < 0) {
      stop(paste("Invalid standard deviation for arm '", arm, "': must be a single non-negative finite numeric value.", sep=""))
    }
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed))) {
    stop("seed must be a single integer or NULL.")
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
  ncore <- if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "") {
    2 # Use 2 cores in R CMD check
  } else {
    max(1, detectCores() - 1) # Use max available cores minus 1, or at least 1
  }

  cl <- tryCatch({
    makeCluster(ncore)
  }, error = function(e) {
    warning(paste("Could not create parallel cluster with", ncore, "cores. Falling back to sequential processing:", e$message))
    return(NULL) # Indicate fallback
  })

  parallel_enabled <- !is.null(cl)

  if (parallel_enabled) {

    on.exit(stopCluster(cl), add = TRUE)

    # Export necessary objects and load libraries to cluster environment
    clusterExport(cl, list("n_list", "mu_list", "sigma_list", "arm_names", ".data", "!!", "sym"),
                  envir = environment())

    if (!is.null(seed)) {
      set.seed(seed)
      clusterSetRNGStream(cl, seed)
    } else {
      clusterSetRNGStream(cl)
    }

    clusterEvalQ(cl, {
      # Load necessary libraries on each worker
      library(dplyr)
      library(tidyr) # Needed for pivot_wider within worker if done there (it's done later)
      library(data.table)
      library(stats) # For rnorm, sd
      library(rlang) # For !! sym
    })

    cat(sprintf("Cluster setup complete with %d cores. Starting parallel computation...\n", ncore))

    # --- Parallel Simulation Execution ---
    results <- parLapply(cl, 1:nsim, function(nrep) {

      # Generate data for each arm using mapply
      simulated_values_list <- mapply(stats::rnorm,
                                      n = n_list,
                                      mean = mu_list,
                                      sd = sigma_list,
                                      SIMPLIFY = FALSE)

      data_temp <- data.table::rbindlist(lapply(names(simulated_values_list), function(arm_name) {
        data.table::data.table(
          arm = arm_name,
          value = simulated_values_list[[arm_name]]
        )
      }))

      data_temp[, nsim := nrep] # Use := for data.table assignment
      data_temp[, id := seq_len(.N), by = arm] # .N is number of rows in data.table, by=arm assigns id per arm

      summary_stats <- data_temp[, list(
        nsim = unique(nsim), # Should be just one unique value per replicate
        n = .N,
        mu_hat = mean(value),
        sd = stats::sd(value) # Use stats::sd explicitly
      ), by = arm]

      summary_stats[, s := sd / sqrt(n)] # Use := for data.table assignment

      list(data = data_temp, summary = summary_stats)
    })

  } else { # Sequential Processing Fallback
    cat("Sequential computation enabled.\n")
    # --- Sequential Simulation Execution ---
    # If cluster creation failed or ncore=1, run sequentially
    results <- vector("list", nsim)
    # Set seed for sequential run if provided
    if (!is.null(seed)) {
      set.seed(seed)
    }

    for (nrep in 1:nsim) {
      # Generate data for each arm using mapply
      simulated_values_list <- mapply(stats::rnorm,
                                      n = n_list,
                                      mean = mu_list,
                                      sd = sigma_list,
                                      SIMPLIFY = FALSE)

      # Create individual data.table for the current replicate
      data_temp <- data.table::rbindlist(lapply(names(simulated_values_list), function(arm_name) {
        data.table::data.table(
          arm = arm_name,
          value = simulated_values_list[[arm_name]]
        )
      }))

      # Add nsim and id columns after combining arms
      data_temp[, nsim := nrep] # Use := for data.table assignment
      data_temp[, id := seq_len(.N), by = arm] # .N is number of rows in data.table, by=arm assigns id per arm


      # Calculate summary statistics for the current replicate using data.table
      summary_stats <- data_temp[, list(
        nsim = unique(nsim), # Should be just one unique value per replicate
        n = .N,
        mu_hat = mean(value),
        sd = stats::sd(value) # Use stats::sd explicitly
      ), by = arm]

      # Calculate standard error of the mean (s = sd / sqrt(n))
      summary_stats[, s := sd / sqrt(n)] # Use := for data.table assignment

      results[[nrep]] <- list(data = data_temp, summary = summary_stats)
    }
  }

  # --- Combine Results from all Simulations ---
  data_all <- data.table::rbindlist(lapply(results, "[[", "data"))
  summary_all_long <- data.table::rbindlist(lapply(results, "[[", "summary")) # Combine summaries (long format)

  summary_all_wide <- summary_all_long %>%
    pivot_wider(names_from = arm, values_from = c(n, mu_hat, sd, s),
                names_glue = "{arm}.{.value}")

  summary_cols <- names(summary_all_wide)
  summary_all_wide <- summary_all_wide %>% select(nsim, everything())

  return(list(
    data = data_all, # Returns data.table
    summary = data.table::as.data.table(summary_all_wide),
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
#'   * `data`: A data.table with columns `nsim`, `id`, `arm`, `value`, containing all simulated DpR values.
#'   * `summary`: A data.table with columns `nsim` and arm-specific summary statistics (n, mu_hat, sd, s).
#'   * `true_value`: A data frame with true mean differences (currently only 'treatment' vs 'control').
#'
#' @importFrom dplyr %>% group_by summarise mutate
#' @importFrom tidyr pivot_wider
#' @importFrom data.table data.table rbindlist
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply
#' @importFrom stats rbinom rgamma sd
#' @importFrom rlang .data !! sym
#' @export
#' @seealso `bayesian_lalonde_decision`
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
