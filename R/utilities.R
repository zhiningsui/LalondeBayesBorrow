#' Calculate Response Frequencies for Binary Endpoints
#'
#' This function calculates the frequency of binary endpoint responses (0/1) across arms for multiple simulation repetitions.
#' It groups the data by arm and simulation repetition to count the number of responses and the total number of patients in each arm.
#'
#' @param data A data frame. Individual-level binary outcomes (0/1) for all arms and repetitions of the simulation.
#'
#' @return A data frame. Each row represents the number of responses and the total number of patients in each arm for each simulation repetition.
#' The data frame includes the following columns:
#'   - `nsim`: Simulation repetition number.
#'   - `arm.n`: Total number of patients in each arm.
#'   - `arm.count`: Total number of responses in each arm.
#' @export
#'
#' @examples
#' data <- data.frame(nsim = rep(1:10, each = 200),
#'                    arm = rep(c("A","B"), each = 100, times = 10),
#'                    resp = stats::rbinom(2000, 1, 0.5))
#' freq_binary(data = data)
#'
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data
freq_binary <- function(data) {
  freq <- data %>%
    group_by(.data$nsim, .data$arm, .data$resp) %>%
    summarise(count = n(), .groups = 'drop') %>%
    complete(.data$nsim, .data$arm, .data$resp, fill = list(count = 0)) %>%
    group_by(.data$nsim, .data$arm) %>%
    mutate(n = sum(count)) %>%
    ungroup() %>%
    dplyr::filter(.data$resp == 1) %>%
    group_by(.data$nsim, .data$arm) %>%
    pivot_wider(names_from = .data$arm, values_from = c(.data$n, .data$count),
                names_glue = "{arm}.{.value}") %>%
    select(.data$nsim, ends_with(".n"), ends_with(".count"))

  return(freq = freq)
}

#' Simulate PTSS values from a transformed Beta distribution
#'
#' @param n Integer. Number of observations to simulate.
#' @param mu Numeric. Target mean PTSS on the [-100, 100] scale.
#' @param sigma Numeric. Target standard deviation PTSS on the [-100, 100] scale.
#' @param a Numeric. Lower bound of PTSS scale (default = -100).
#' @param b Numeric. Upper bound of PTSS scale (default = 100).
#'
#' @return A numeric vector of simulated PTSS values.
#' @export
#'
#' @examples
#' set.seed(123)
#' simulate_ptss_transformed_beta(n = 100, mu = 20, sigma = 20)
simulate_ptss <- function(n, mu, sigma, a = -1, b = 1) {
  mu_star <- (mu - a) / (b - a)
  sigma_star <- sigma / (b - a)

  # Method-of-moments
  nu <- (mu_star * (1 - mu_star)) / (sigma_star^2) - 1

  if (nu <= 0 || is.nan(nu)) {
    warning(sprintf("Invalid combination of mu = %.3f and sigma = %.3f: nu = %.3f is not positive.", mu, sigma, nu))
    return(rep(NA, n))
  }

  alpha <- mu_star * nu
  beta <- (1 - mu_star) * nu

  # Simulate
  x <- rbeta(n, shape1 = alpha, shape2 = beta)
  y <- a + (b - a) * x
  return(y)
}

#' Validate PTSS parameters for transformed Beta distribution on [-1, 1]
#'
#' @param mu A numeric vector of means on [-1, 1]
#' @param sigma A numeric vector of standard deviations (must be positive)
#' @param a Lower bound of the interval (default: -1)
#' @param b Upper bound of the interval (default: 1)
#'
#' @return A logical vector indicating whether each (mu, sigma) pair is valid.
#'         Also throws warnings for invalid cases.
#'
#' @export
validate_ptss_params <- function(mu, sigma, a = -1, b = 1) {
  if (length(mu) != length(sigma)) stop("mu and sigma must be the same length")
  mu_star <- (mu - a) / (b - a)
  sigma_star <- sigma / (b - a)

  sigma_max <- sqrt(mu_star * (1 - mu_star))  # from Beta variance constraint
  is_valid <- (sigma_star < sigma_max) & (sigma_star > 0) & (mu_star > 0) & (mu_star < 1)

  for (i in which(!is_valid)) {
    warning(sprintf("Invalid PTSS parameters: mu = %.3f, sigma = %.3f; sigma must be < %.3f",
                    mu[i], sigma[i], sigma_max[i]))
  }

  return(is_valid)
}

#' Retrieve All Functions from an Environment
#'
#' This function returns the names of all functions within a specified environment.
#'
#' @param env An environment. The environment to search for functions. The default is the current environment.
#'
#' @return A character vector. Names of all functions found in the specified environment.
#' @export
get_all_functions <- function(env = environment()) {
  all_objects <- ls(env)
  functions <- sapply(all_objects, function(obj_name) is.function(get(obj_name, envir = env)))
  return(names(functions)[functions])
}


#' Convert to RBesT Mixture Distribution
#'
#' This function generates a mixture distribution in the format required by the RBesT package from a given set of parameters.
#' It supports both continuous and binary endpoints, constructing either normal or beta mixture distributions based on the endpoint type.
#'
#' @param post A data frame. Parameters for the mixture distribution, including the weight (`w`) of each component, and the distribution parameters for each component:
#'   - For continuous endpoints: `w`, `mu1`, `sigma1`, `mu2`, and `sigma2`.
#'   - For binary endpoints: `w`, `a1`, `b1`, `a2`, and `b2`.
#' @param endpoint A string. Specifies the type of endpoint ("continuous" for normal distributions, "binary" for beta distributions).
#'
#' @return An RBesT mixture distribution object. The object will be of class `mixnorm` for continuous endpoints or `mixbeta` for binary endpoints.
#' @export
#'
#' @examples
#' # Example for continuous endpoint
#' post <- data.frame(w = 0.5, mu1 = 0, sigma1 = 1, mu2 = 0.5, sigma2 = 1)
#' mix_dist <- convert_RBesT_mix(post, endpoint = "continuous")
#'
#' # Example for binary endpoint
#' post <- data.frame(w = 0.5, a1 = 1, b1 = 1, a2 = 2, b2 = 2)
#' mix_dist <- convert_RBesT_mix(post, "binary")
#' @import RBesT
convert_RBesT_mix <- function(post, endpoint){
  if (endpoint == "continuous") {
    if(!is.data.frame(post)) post <- as.data.frame(t(post))
    # Define the normal mixture distributions using RBesT
    if (post$w == 0) {
      post = mixnorm(c(1-unname(post$w), unname(post$mu2), unname(post$sigma2)))
    } else if (post$w == 1){
      post = mixnorm(c(unname(post$w), unname(post$mu1), unname(post$sigma1)))
    }
    else {
      post = mixnorm(c(unname(post$w), unname(post$mu1), unname(post$sigma1)),
                     c(1-unname(post$w), unname(post$mu2), unname(post$sigma2)))
    }
  } else if (endpoint == "binary") {
    # Define the beta mixture distributions using RBesT
    if (post$w == 0) {
      post = mixbeta(c(1-unname(post$w), unname(post$a2), unname(post$b2)))
    } else if (post$w == 1){
      post = mixbeta(c(unname(post$w), unname(post$a1), unname(post$b1)))
    }else {
      post = mixbeta(c(unname(post$w), unname(post$a1), unname(post$b1)),
                     c(1-unname(post$w), unname(post$a2), unname(post$b2)))
    }
  }
  return(post)
}

#' Create Data Generation Parameters
#'
#' This function processes and validates a set of parameters for generating data
#' related to treatment and control groups, along with their historical counterparts.
#' It checks for missing values and replaces them with `NA` while issuing warnings.
#' The function supports different endpoints, and the interpretation of the `prob` parameter
#' varies depending on the specified endpoint:
#'   - For the "g-score" endpoint, `prob` refers to the probability of zero g-scores.
#'   - For the "OR" endpoint, `prob` refers to the probability of an event occurring.
#'
#' @param params A list. Parameters for data generation.
#'   - `trt_n`: A scalar. Sample size for the treatment group.
#'   - `ctrl_n`: A scalar. Sample size for the control group. If not provided, it will be set to match `trt_n`.
#'   - `trt_prob`: A scalar. Probability for the treatment group, interpretation based on the `endpoint`.
#'   - `ctrl_prob`: A scalar. Probability for the control group, interpretation based on the `endpoint`.
#'   - `trt_mu`: A scalar. Mean value for the treatment group (for non-zero values in g-scores or other continuous endpoints).
#'   - `ctrl_mu`: A scalar. Mean value for the control group.
#'   - `trt_sigma`: A scalar. Standard deviation for the treatment group (for non-zero values in g-scores or other continuous endpoints).
#'   - `ctrl_sigma`: A scalar. Standard deviation for the control group.
#'   - Additional parameters for historical data (e.g., `trt_h_n`, `ctrl_h_prob`) may also be included.
#' @param endpoint A string. Specifies the type of endpoint ("g-score" or "OR").
#'
#' @return A nested list. Processed parameters for each group, including `treatment`, `control`, `treatment_h`, and `control_h` if historical data is provided.
#' Missing values are replaced with `NA`, and a warning is issued for each missing parameter.
#' @export
#'
#' @examples
#' params <- list(trt_n = 100, trt_prob = 0.6, trt_mu = 1.5, trt_sigma = 0.3,
#'                ctrl_n = 100, ctrl_prob = 0.4, ctrl_mu = 1.2, ctrl_sigma = 0.2)
#' create_data_gen_params(params, "g-score")
create_data_gen_params <- function(params, endpoint) {
  check_param <- function(param, param_name) {
    if (!is.null(param)) {
      return(param)
    } else {
      warning(paste("Warning: Parameter", param_name, "is missing. Replacing with NA."))
      return(NA)
    }
  }

  safe_log <- function(param, param_name) {
    if (!is.null(param)) {
      return(log(param))
    } else {
      warning(paste("Warning: Parameter", param_name, "is missing. Replacing with NA."))
      return(NA)
    }
  }

  beta_moments <- function(mu, sigma, a = -1, b = 1) {
    mu_star <- (mu - a) / (b - a)
    sigma_star <- sigma / (b - a)
    nu <- mu_star * (1 - mu_star) / sigma_star^2 - 1
    alpha <- mu_star * nu
    beta <- (1 - mu_star) * nu
    return(c(alpha = alpha, beta = beta))
  }

  if (endpoint == "ptss") {
    a <- -1
    b <- 1

    get_ptss_entry <- function(mu, sigma, n, name) {
      ab <- beta_moments(mu, sigma, a, b)
      list(
        n = n,
        mu = mu,
        sigma = sigma,
        alpha = unname(ab["alpha"]),
        beta = unname(ab["beta"]),
        name = name
      )
    }

    param_list <- list(
      treatment = get_ptss_entry(
        mu = check_param(params$trt_mu, "trt_mu"),
        sigma = check_param(params$trt_sigma, "trt_sigma"),
        n = check_param(params$trt_n, "trt_n"),
        name = "treatment"
      ),
      control = get_ptss_entry(
        mu = check_param(params$ctrl_mu, "ctrl_mu"),
        sigma = check_param(params$ctrl_sigma, "ctrl_sigma"),
        n = check_param(params$ctrl_n, "ctrl_n"),
        name = "control"
      ),
      treatment_h = get_ptss_entry(
        mu = check_param(params$trt_h_mu, "trt_h_mu"),
        sigma = check_param(params$trt_h_sigma, "trt_h_sigma"),
        n = check_param(params$trt_h_n, "trt_h_n"),
        name = "treatment_h"
      ),
      control_h = get_ptss_entry(
        mu = check_param(params$ctrl_h_mu, "ctrl_h_mu"),
        sigma = check_param(params$ctrl_h_sigma, "ctrl_h_sigma"),
        n = check_param(params$ctrl_h_n, "ctrl_h_n"),
        name = "control_h"
      )
    )

  } else if (endpoint == "g-score" || endpoint == "OR") {
    get_entry <- function(prefix, name) {
      list(
        n = check_param(params[[paste0(prefix, "_n")]], paste0(prefix, "_n")),
        prob = check_param(params[[paste0(prefix, "_prob")]], paste0(prefix, "_prob")),
        mu = safe_log(params[[paste0(prefix, "_mu")]], paste0(prefix, "_mu")),
        sigma = check_param(params[[paste0(prefix, "_sigma")]], paste0(prefix, "_sigma")),
        name = name
      )
    }

    param_list <- list(
      treatment = get_entry("trt", "treatment"),
      control = get_entry("ctrl", "control"),
      treatment_h = get_entry("trt_h", "treatment_h"),
      control_h = get_entry("ctrl_h", "control_h")
    )
  }

  # Filter out arms with missing required values
  if (endpoint == "g-score") {
    required_keys <- c("n", "prob", "mu", "sigma")
  } else if (endpoint == "OR") {
    required_keys <- c("n", "prob")
  } else if (endpoint == "ptss") {
    required_keys <- c("n", "mu", "sigma", "alpha", "beta")
  } else {
    stop("Unsupported endpoint type.")
  }

  param_list_filtered <- Filter(function(x) !any(sapply(x[required_keys], is.na)), param_list)


  # Filter out entries with missing n or prob parameters
  if(endpoint == "g-score"){
    param_list_filtered <- Filter(function(x) !any(sapply(x, is.na)), param_list)
  } else if (endpoint == "OR"){
    param_list_filtered <- Filter(function(x) !any(sapply(x[c("n", "prob")], is.na)), param_list)
  }

  param_list_filtered <- Filter(function(x) !any(sapply(x[required_keys], is.na)), param_list)

  return(param_list_filtered)
}


