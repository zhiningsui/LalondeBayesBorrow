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
#' @import stats
freq_binary <- function(data) {
  freq <- data %>%
    group_by(nsim, arm, resp) %>%
    summarise(count = n(), .groups = 'drop') %>%
    complete(nsim, arm, resp, fill = list(count = 0)) %>%
    group_by(nsim, arm) %>%
    mutate(n = sum(count)) %>%
    ungroup() %>%
    dplyr::filter(resp == 1) %>%
    group_by(nsim, arm) %>%
    pivot_wider(names_from = arm, values_from = c(n, count),names_glue = "{arm}.{.value}") %>%
    select(nsim, ends_with(".n"), ends_with(".count"))

  return(freq = freq)
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
  param_list <- list(
    treatment = list(
      n = check_param(params$trt_n, "trt_n"),
      prob = check_param(params$trt_prob, "trt_prob"),
      mu = safe_log(params$trt_mu, "trt_mu"),
      sigma = check_param(params$trt_sigma, "trt_sigma"),
      name = "treatment"
    ),
    control = list(
      n = check_param(params$ctrl_n, "ctrl_n"),
      prob = check_param(params$ctrl_prob, "ctrl_prob"),
      mu = safe_log(params$ctrl_mu, "ctrl_mu"),
      sigma = check_param(params$ctrl_sigma, "ctrl_sigma"),
      name = "control"
    ),
    treatment_h = list(
      n = check_param(params$trt_h_n, "trt_h_n"),
      prob = check_param(params$trt_h_prob, "trt_h_prob"),
      mu = safe_log(params$trt_h_mu, "trt_h_mu"),
      sigma = check_param(params$trt_h_sigma, "trt_h_sigma"),
      name = "treatment_h"
    ),
    control_h = list(
      n = check_param(params$ctrl_h_n, "ctrl_h_n"),
      prob = check_param(params$ctrl_h_prob, "ctrl_h_prob"),
      mu = safe_log(params$ctrl_h_mu, "ctrl_h_mu"),
      sigma = check_param(params$ctrl_h_sigma, "ctrl_h_sigma"),
      name = "control_h"
    )
  )
  # Filter out entries with missing n or prob parameters
  if(endpoint == "g-score"){
    param_list_filtered <- Filter(function(x) !any(sapply(x, is.na)), param_list)
  } else if (endpoint == "OR"){
    param_list_filtered <- Filter(function(x) !any(sapply(x[c("n", "prob")], is.na)), param_list)
  }

  return(param_list_filtered)
}



