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
#' @examples
#' # Example of generating g-scores for a study with 2 arms, 100 patients in each arm,
#' # and 10 simulation repetitions:
#' n_list <- c(treatment = 100, control = 100)
#' prob_list <- c(treatment = 0.2, control = 0.3)
#' mu_list <- c(treatment = 0.5, control = 0.6)
#' sigma_list <- c(treatment = 0.1, control = 0.2)
#' arm_names <- c(treatment = "treatment", control = "control")
#' sim_data_gscore <- data_gen_gscore(n_arms = 2, nsim = 10, n_list = n_list,
#'                                    prob_list = prob_list, mu_list = mu_list,
#'                                    sigma_list = sigma_list, arm_names = arm_names)
#'
#' @import dplyr
#' @import tidyr
#' @import parallel
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
#' This function simulates binary endpoint data (0/1) for multiple study arms across several simulation repetitions.
#' Each arm's data is generated based on the specified probabilities of response, and the function returns
#' individual-level data as well as summary statistics for the true response rates and response frequencies.
#'
#' @param n_arms A scalar. The number of arms to be generated in the simulation.
#' @param nsim A scalar. The number of repetitions of the simulation to be performed.
#' @param n_list A list. Each element represents the number of patients in each arm.
#' @param prob_list A list. Each element represents the probability of response (binary outcome = 1) for each arm.
#' @param arm_names A list. Each element represents the name of the corresponding arm.
#'
#' @return A list of data frames:
#'   - \code{data}: A data frame containing individual-level binary outcomes (0/1) for all arms and all repetitions of the simulation.
#'   - \code{true_value}: A data frame of true response rates for each arm, along with the difference in response rates between the treatment and control arms (if both are present).
#'   - \code{freq}: A data frame of response frequencies, where each row represents the number of responses and the total number of patients in each arm for each simulation repetition.
#' @export
#'
#' @examples
#' # Example of generating binary response data for 2 arms, 100 patients per arm,
#' # and 10 simulation repetitions:
#' n_list <- c(treatment = 50, control = 50, control_h = 200)
#' prob_list <- c(treatment = 0.4, control = 0.2, control_h = 0.3)
#' arm_names <- c(treatment = "treatment", control = "control", control_h = "control_h")
#' sim_data_bin <- data_gen_binary(n_arms = 3, nsim = 10, n_list = n_list,
#'                                 prob_list = prob_list, arm_names = arm_names)
#'
#' @import dplyr
#' @import tidyr
#' @import parallel
#' @importFrom data.table data.table rbindlist setcolorder
#' @importFrom rlang := .data
data_gen_binary <- function(n_arms = 2, nsim = 10, n_list = NULL, prob_list = NULL, arm_names = NULL) {
  rate_true <- data.frame(arm = arm_names,
                          rate_true = unlist(prob_list)) %>%
      pivot_wider(names_from = .data$arm, values_from = .data$rate_true,
                  names_glue = "{arm}.rate_true")

  if (sum(arm_names %in% c("treatment", "control")) >= 2) {
    rate_true <- rate_true %>%
      mutate(compare_true = !!sym(paste0(arm_names[arm_names == "treatment"], ".rate_true")) - !!sym(paste0(arm_names[arm_names == "control"], ".rate_true")))
  } else {
    rate_true <- rate_true %>%
      mutate(compare_true = NA)
  }

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncore <- 2
  } else {
    # use all cores in devtools::test()
    ncore <- parallel::detectCores() - 1
  }
  cl <- makeCluster(ncore)
  clusterExport(cl, list("n_list", "prob_list", "arm_names", "n_arms"), envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(data.table)
  })

  results <- parLapply(cl, 1:nsim, function(nrep) {
   data_temp <- data.table(
      nsim = rep(nrep, each = sum(n_list)),
      id = unlist(lapply(n_list, seq)),
      arm = unlist(as.list(mapply(rep, arm_names, n_list))),
      resp = unlist(as.list(mapply(stats::rbinom, n_list, as.list(rep(1, n_arms)), prob_list)))
    )
   data_temp$arm <- factor(data_temp$arm, levels = arm_names)
   data_temp
  })

  stopCluster(cl)

  data <- rbindlist(results)
  freq <- freq_binary(data)

  return(list(data = data, true_value = rate_true, freq = freq))

}

data_gen_continuous <- function(n_arms = 2, nsim = 10, n_list = NULL, mu_list = NULL, sigma_list = NULL,arm_names = NULL) {
  mean_true <- data.frame(arm = arm_names,
                          mean_true = unlist(mu_list)) %>%
    pivot_wider(names_from = .data$arm, values_from = .data$mean_true,
                names_glue = "{arm}.mean_true")

  if (sum(arm_names %in% c("treatment", "control")) >= 2) {
    mean_true <- mean_true %>%
      mutate(compare_true = !!sym(paste0(arm_names[arm_names == "treatment"], ".mean_true")) - !!sym(paste0(arm_names[arm_names == "control"], ".mean_true")))
  } else {
    mean_true <- mean_true %>%
      mutate(compare_true = NA)
  }

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncore <- 2
  } else {
    # use all cores in devtools::test()
    ncore <- parallel::detectCores() - 1
  }
  cl <- makeCluster(ncore)
  clusterExport(cl, list("n_list", "mu_list", "sigma_list", "arm_names", "n_arms"),
                envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(data.table)
  })

  results <- parLapply(cl, 1:nsim, function(nrep) {

    data_temp <- data.table(
      nsim = rep(nrep, each = sum(n_list)),
      id = unlist(lapply(n_list, seq)),
      arm = unlist(as.list(mapply(rep, arm_names, n_list))),
      value = unlist(mapply(rnorm, n = n_list, mean = mu_list, sd = sigma_list, SIMPLIFY = FALSE))
    )

    summary_stats <- data_temp %>%
      group_by(arm) %>%
      summarise(
        nsim = unique(nsim),
        n = n(),
        mu_hat = mean(value),
        sd = sd(value),
        .groups = "drop"
      ) %>%
      mutate(s = sd / sqrt(n))

    list(data = data_temp, summary = summary_stats)
  })

  stopCluster(cl)

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


#' Generate PTSS data across arms and simulations using transformed Beta
#'
#' @param n_arms Number of arms
#' @param nsim Number of simulations
#' @param n_list Named list of sample sizes per arm
#' @param mu_list Named list of means per arm
#' @param sigma_list Named list of standard deviations per arm
#' @param arm_names Character vector of arm names
#'
#' @return A list with elements `data`, `summary`, and `true_value`
#'
data_gen_ptss <- function(n_arms = 2, nsim = 10,
                          n_list = NULL, mu_list = NULL, sigma_list = NULL,
                          arm_names = NULL) {
  a <- -1
  b <- 1

  # Calculate true values
  # true_value <- data.frame(arm = arm_names) %>%
  #   rowwise() %>%
  #   mutate(true_value = mu_list[arm]) %>%
  #   pivot_wider(names_from = .data$arm, values_from = .data$true_value, names_glue = "{arm}.true_value")

  true_value <- data.frame(arm = arm_names) %>%
    rowwise() %>%
    mutate(true_value = a + (b - a) * ((mu_list[arm] - a) / (b - a))) %>%
    pivot_wider(names_from = .data$arm, values_from = .data$true_value, names_glue = "{arm}.true_value")

  if (all(c("treatment", "control") %in% arm_names)) {
    true_value <- true_value %>%
      mutate(compare_true = treatment.true_value - control.true_value)
  } else {
    true_value <- true_value %>% mutate(compare_true = NA)
  }

  ncore <- max(parallel::detectCores() - 1, 1)
  cl <- parallel::makeCluster(ncore)
  parallel::clusterExport(cl, varlist = c("simulate_ptss",
                                          "n_list", "mu_list", "sigma_list",
                                          "arm_names", "a", "b"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(data.table)
  })

  results <- parallel::parLapply(cl, 1:nsim, function(nrep) {
    data_temp <- data.table::data.table(
      nsim = rep(nrep, sum(n_list)),
      id = unlist(lapply(n_list, seq)),
      arm = factor(unlist(mapply(rep, arm_names, n_list, SIMPLIFY = TRUE, USE.NAMES = FALSE)),
                   levels = arm_names),
      ptss = unlist(as.list(mapply(simulate_ptss, n_list, mu_list, sigma_list)))
    )

    summary_stats <- data_temp %>%
      group_by(arm) %>%
      summarise(
        nsim = unique(nsim),
        n = n(),
        mean_ptss = mean(ptss),
        sd_ptss = sd(ptss),
        .groups = "drop"
      )

    list(data = data_temp, summary = summary_stats)
  })

  parallel::stopCluster(cl)

  data_all <- data.table::rbindlist(lapply(results, "[[", "data"))
  summary_all <- data.table::rbindlist(lapply(results, "[[", "summary"))

  return(list(
    data = data_all,
    summary = summary_all,
    true_value = true_value
  ))
}
