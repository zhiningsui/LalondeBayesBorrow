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
#' This function generates a mixture distribution object compatible with the
#' RBesT package from a data frame containing component parameters.
#' It supports the creation of both normal (`mixnorm`) and beta (`mixbeta`)
#' mixture distributions based on the specified endpoint type. The function
#' is designed to handle mixture models with up to two components based on the
#' provided structure of the input data frame.
#'
#' @param post A data frame or a list of parameters. This object should contain
#'   the parameters defining the mixture components. The expected structure
#'   depends on the `endpoint` type:
#'   * For `endpoint = "continuous"` (Normal mixture): Expected columns/elements are
#'     `w` (weight of the first component), `mu1` (mean of the first component),
#'     `sigma1` (standard deviation of the first component), `mu2` (mean of the
#'     second component), and `sigma2` (standard deviation of the second component).
#'     The weight of the second component is implicitly `1 - w`.
#'   * For `endpoint = "binary"` (Beta mixture): Expected columns/elements are
#'     `w` (weight of the first component), `a1` (alpha parameter of the first component),
#'     `b1` (beta parameter of the first component), `a2` (alpha parameter of the
#'     second component), and `b2` (beta parameter of the second component).
#'     The weight of the second component is implicitly `1 - w`.
#'   The function currently assumes a single row/set of parameters for the mixture,
#'   representing either a one- or two-component mixture.
#' @param endpoint A character string. Specifies the type of endpoint and thus
#'   the type of mixture distribution to create. Must be either `"continuous"`
#'   for a normal mixture (`mixnorm`) or `"binary"` for a beta mixture (`mixbeta`).
#'
#' @return An RBesT mixture distribution object, specifically of class `mixnorm`
#'   if `endpoint` is `"continuous"` or `mixbeta` if `endpoint` is `"binary"`.
#'   Returns an error if the `endpoint` is not specified correctly or if required
#'   parameters are missing or invalid.
#' @export
#'
#' @import RBesT
#'
#' @examples
#' # Example for a continuous endpoint (Normal mixture)
#' post_continuous <- data.frame(w = 0.6, mu1 = 10, sigma1 = 2, mu2 = 15, sigma2 = 3)
#' mix_norm <- convert_RBesT_mix(post_continuous, endpoint = "continuous")
#' print(mix_norm)
#'
#' # Example for a binary endpoint (Beta mixture)
#' post_binary <- list(w = 0.8, a1 = 5, b1 = 15, a2 = 10, b2 = 10)
#' mix_beta <- convert_RBesT_mix(post_binary, endpoint = "binary")
#' print(mix_beta)
#'
#' # Example for a single-component mixture (weight = 1)
#' post_single <- data.frame(w = 1, mu1 = 12, sigma1 = 1.5, mu2 = NA, sigma2 = NA)
#' mix_single <- convert_RBesT_mix(post_single, endpoint = "continuous")
#' print(mix_single)
#'
#' # Example for a single-component mixture (weight = 0)
#' post_single_comp2 <- list(w = 0, mu1 = NA, sigma1 = NA, mu2 = 8, sigma2 = 2.5)
#' mix_single_comp2 <- convert_RBesT_mix(post_single_comp2, endpoint = "continuous")
#' print(mix_single_comp2)
#'
#' @seealso \code{\link[RBesT]{mixnorm}}, @seealso \code{\link[RBesT]{mixbeta}}, @seealso \code{\link[RBesT]{mix}}
convert_RBesT_mix <- function(post, endpoint){
  # --- Input Validation ---
  if (!is.data.frame(post) && !is.list(post)) {
    stop("Input 'post' must be a data frame or a list.")
  }

  if (is.data.frame(post) && nrow(post) != 1) {
    if (nrow(post) > 1) {
      warning("Input 'post' has more than one row. Only the first row will be used.")
    }
    post <- post[1, , drop = FALSE]
  } else if (is.list(post) && length(post) > 0 && is.list(post[[1]])) {
    if (length(post) > 1) {
      warning("Input 'post' is a list with multiple elements. Only the first element will be used.")
    }
    post <- post[[1]]
  }

  if (!is.character(endpoint) || length(endpoint) != 1 || !(endpoint %in% c("continuous", "binary"))) {
    stop("Input 'endpoint' must be a single string, either 'continuous' or 'binary'.")
  }

  required_params <- if (endpoint == "continuous") {
    c("w", "mu1", "sigma1", "mu2", "sigma2")
  } else {
    c("w", "a1", "b1", "a2", "b2")
  }

  if (!all(required_params %in% names(post))) {
    missing_params <- setdiff(required_params, names(post))
    stop(paste("Input 'post' is missing required parameters for endpoint '", endpoint, "': ",
               paste(missing_params, collapse = ", "), sep = ""))
  }

  for (param in required_params) {
    if (!is.numeric(post[[param]]) && !is.na(post[[param]])) {
      stop(paste("Parameter '", param, "' in 'post' must be numeric.", sep = ""))
    }
  }

  # Validate weights
  if (is.na(post$w) || post$w < 0 || post$w > 1) {
    stop("Parameter 'w' in 'post' must be a numeric value between 0 and 1.")
  }

  # --- Create RBesT Mixture Object ---
  w <- unname(post$w)

  if (endpoint == "continuous") {
    mu1 <- unname(post$mu1)
    sigma1 <- unname(post$sigma1)
    mu2 <- unname(post$mu2)
    sigma2 <- unname(post$sigma2)

    # Validate standard deviations
    if ((w > 0 && (is.na(sigma1) || sigma1 <= 0)) || (w < 1 && (is.na(sigma2) || sigma2 <= 0))) {
      stop("Standard deviations (sigma1, sigma2) must be positive for components with non-zero weights.")
    }

    if (w == 0) {
      if (is.na(mu2) || is.na(sigma2)) {
        stop("Parameters for the second component (mu2, sigma2) are missing when w = 0.")
      }
      mixture <- mixnorm(c(1, mu2, sigma2))
    } else if (w == 1) {
      if (is.na(mu1) || is.na(sigma1)) {
        stop("Parameters for the first component (mu1, sigma1) are missing when w = 1.")
      }
      mixture <- mixnorm(c(1, mu1, sigma1))
    } else {
      if (is.na(mu1) || is.na(sigma1) || is.na(mu2) || is.na(sigma2)) {
        stop("Parameters for both components (mu1, sigma1, mu2, sigma2) are required for a two-component mixture (0 < w < 1).")
      }
      mixture <- mixnorm(c(w, mu1, sigma1), c(1 - w, mu2, sigma2))
    }
  } else if (endpoint == "binary") {
    a1 <- unname(post$a1)
    b1 <- unname(post$b1)
    a2 <- unname(post$a2)
    b2 <- unname(post$b2)

    # Validate beta parameters (must be positive)
    if ((w > 0 && (is.na(a1) || a1 <= 0 || is.na(b1) || b1 <= 0)) || (w < 1 && (is.na(a2) || a2 <= 0 || is.na(b2) || b2 <= 0))) {
      stop("Beta parameters (a, b) must be positive for components with non-zero weights.")
    }

    if (w == 0) {
      if (is.na(a2) || is.na(b2)) {
        stop("Parameters for the second component (a2, b2) are missing when w = 0.")
      }
      mixture <- mixbeta(c(1, a2, b2))
    } else if (w == 1) {
      if (is.na(a1) || is.na(b1)) {
        stop("Parameters for the first component (a1, b1) are missing when w = 1.")
      }
      mixture <- mixbeta(c(1, a1, b1))
    } else {
      if (is.na(a1) || is.na(b1) || is.na(a2) || is.na(b2)) {
        stop("Parameters for both components (a1, b1, a2, b2) are required for a two-component mixture (0 < w < 1).")
      }
      mixture <- mixbeta(c(w, a1, b1), c(1 - w, a2, b2))
    }
  }
  return(mixture)
}

#' Create Data Generation Parameters
#'
#' This function processes and validates a set of parameters for generating data
#' for different study arms (treatment, control, and their historical counterparts).
#' It checks for missing parameters, replaces them with `NA` while issuing warnings,
#' and structures the valid parameters into a nested list format suitable for
#' data simulation functions. The function supports various endpoint types,
#' each requiring a specific set of parameters and potentially applying
#' transformations (like logarithm for means in "g-score").
#'
#' Arms (treatment, control, historical treatment, historical control) are
#' included in the output list only if all their required parameters for the
#' specified `endpoint` are present (not `NULL` and not resulting in `NA` after
#' initial checks), and if the sample size `n` is positive.
#'
#' The interpretation of parameters depends on the `endpoint`:
#' * **continuous**: Parameters for a normal continuous outcome include a mean (`mu`) and a standard deviation (`sigma`). Required parameters per arm: `n`, `mu`, `sigma`.
#' * **binary**: Parameters for a binary outcome include a probability of the event occurring (`p`). Required parameters per arm: `n`, `p`.
#' * **DpR** (Depth of Response): Parameters include a probability of complete tumor shrinkage (`p`), a mean (`mu`), and a standard deviation (`sigma`) for the continuous outcome part. Required parameters per arm: `n`, `p`, `mu`, `sigma`.
#' * **g-score**: Parameters for zero-inflated continuous g-scores include a probability of observing a zero value (`p`), a mean (`mu`, log-transformed internally) and a standard deviation (`sigma`) for the non-zero part. Required parameters per arm: `n`, `p`, `mu`, `sigma`.
#'
#' @param params A list. Contains the parameters for data generation. Expected
#'   parameters follow a naming convention: `[arm_prefix]_[parameter_name]`,
#'   where `[arm_prefix]` can be `trt` (treatment), `ctrl` (control),
#'   `trt_h` (historical treatment), or `ctrl_h` (historical control).
#'   Required `[parameter_name]`s depend on the `endpoint`. Examples include
#'   `n` (sample size), `p` (probability), `mu` (mean), `sigma` (standard deviation).
#' @param endpoint A character string. Specifies the type of endpoint. Must be one
#'   of `"continuous"`, `"binary"`, `"DpR"`, or `"g-score"`.
#'
#' @return A named list where each element corresponds to a study arm
#'   (`treatment`, `control`, `treatment_h`, `control_h`). Each arm's element is
#'   a list containing its data generating parameters (`n`, `p`, `mu`, `sigma`, etc.,
#'   depending on the endpoint) and its `name`. Arms for which required parameters
#'   were missing (`NULL` or resulted in `NA`) or `n` was non-positive are excluded
#'   from the list. Missing parameters for included arms are reported via warnings
#'   and set to `NA`.
#' @export
#'
#' @examples
#' # Example for 'continuous' endpoint with missing historical data
#' params_continuous <- list(
#'   trt_n = 150, trt_mu = 75, trt_sigma = 10,
#'   ctrl_n = 150, ctrl_mu = 70, ctrl_sigma = 11,
#'   trt_h_n = 80, trt_h_mu = 76 # Missing sigma, this arm will be excluded
#' )
#' data_params_continuous <- create_data_gen_params(params_continuous, "continuous")
#' str(data_params_continuous)
#'
#' # Example for 'binary' endpoint
#' params_binary <- list(
#'   trt_n = 200, trt_p = 0.7,
#'   ctrl_n = 200, ctrl_p = 0.5,
#'   ctrl_h_n = 100, ctrl_h_p = 0.55
#' )
#' data_params_binary <- create_data_gen_params(params_binary, "binary")
#' str(data_params_binary)
#'
#' # Example for 'DpR' endpoint
#' params_DpR <- list(
#'   trt_n = 120, trt_p = 0.9, trt_mu = 5, trt_sigma = 1.2, # p=prob of complete shrinkage
#'   ctrl_n = 120, ctrl_p = 0.8, ctrl_mu = 4.5, ctrl_sigma = 1.5
#' )
#' data_params_DpR <- create_data_gen_params(params_DpR, "DpR")
#' str(data_params_DpR)
#'
#' # Example for 'g-score' endpoint
#' params_gscore <- list(
#'   trt_n = 100, trt_p = 0.2, trt_mu = 1.5, trt_sigma = 0.5, # p=prob of zero g-score
#'   ctrl_n = 100, ctrl_p = 0.4, ctrl_mu = 1.2, ctrl_sigma = 0.6
#' )
#' data_params_gscore <- create_data_gen_params(params_gscore, "g-score")
#' str(data_params_gscore)
#'
create_data_gen_params <- function(params, endpoint) {

  check_param <- function(param, param_name) {
    if (is.null(param)) {
      warning(paste("Parameter '", param_name, "' is missing. Replacing with NA.", sep = ""))
      return(NA)
    }
    if (is.numeric(param) && !is.finite(param)) {
      warning(paste("Parameter '", param_name, "' has non-finite numeric value. Replacing with NA.", sep = ""))
      return(NA)
    }
    return(param)
  }

  safe_log <- function(param, param_name) {
    param_checked <- check_param(param, param_name)
    if (is.na(param_checked)) {
      return(NA)
    }
    if (!is.numeric(param_checked)) {
      warning(paste("Parameter '", param_name, "' is not numeric, cannot apply log transformation. Replacing with NA.", sep = ""))
      return(NA)
    }
    if (param_checked <= 0) {
      warning(paste("Parameter '", param_name, "' is not positive, cannot apply log transformation. Replacing with NA.", sep = ""))
      return(NA)
    }
    return(log(param_checked))
  }

  get_arm_params <- function(prefix, name, params, endpoint) {
    arm_list <- list(name = name)

    if (endpoint == "continuous") {
      arm_list$n <- check_param(params[[paste0(prefix, "_n")]], paste0(prefix, "_n"))
      arm_list$mu <- check_param(params[[paste0(prefix, "_mu")]], paste0(prefix, "_mu"))
      arm_list$sigma <- check_param(params[[paste0(prefix, "_sigma")]], paste0(prefix, "_sigma"))
    } else if (endpoint == "binary") {
      arm_list$n <- check_param(params[[paste0(prefix, "_n")]], paste0(prefix, "_n"))
      arm_list$p <- check_param(params[[paste0(prefix, "_p")]], paste0(prefix, "_p"))
    } else if (endpoint == "DpR") {
      arm_list$n <- check_param(params[[paste0(prefix, "_n")]], paste0(prefix, "_n"))
      arm_list$p <- check_param(params[[paste0(prefix, "_p")]], paste0(prefix, "_p"))
      arm_list$mu <- check_param(params[[paste0(prefix, "_mu")]], paste0(prefix, "_mu"))
      arm_list$sigma <- check_param(params[[paste0(prefix, "_sigma")]], paste0(prefix, "_sigma"))
    } else if (endpoint == "g-score") {
      arm_list$n <- check_param(params[[paste0(prefix, "_n")]], paste0(prefix, "_n"))
      arm_list$p <- check_param(params[[paste0(prefix, "_p")]], paste0(prefix, "_p"))
      # Apply log transform to mu for g-score
      arm_list$mu <- safe_log(params[[paste0(prefix, "_mu")]], paste0(prefix, "_mu"))
      arm_list$sigma <- check_param(params[[paste0(prefix, "_sigma")]], paste0(prefix, "_sigma"))
    } else {
      stop("Internal error: Unsupported endpoint type encountered during parameter extraction.")
    }

    return(arm_list)
  }

  # --- Validate endpoint ---
  supported_endpoints <- c("continuous", "binary", "DpR", "g-score")
  if (!is.character(endpoint) || length(endpoint) != 1 || !(endpoint %in% supported_endpoints)) {
    stop(paste("Input 'endpoint' must be a single string from:", paste(supported_endpoints, collapse = ", ")))
  }

  # --- Extract parameters for each potential arm ---
  param_list_raw <- list(
    treatment = get_arm_params("trt", "treatment", params, endpoint),
    control = get_arm_params("ctrl", "control", params, endpoint),
    treatment_h = get_arm_params("trt_h", "treatment_h", params, endpoint),
    control_h = get_arm_params("ctrl_h", "control_h", params, endpoint)
  )

  # --- Define required keys for filtering based on endpoint ---
  required_keys <- switch(endpoint,
                          "continuous" = c("n", "mu", "sigma"),
                          "binary" = c("n", "p"),
                          "DpR" = c("n", "p", "mu", "sigma"),
                          "g-score" = c("n", "p", "mu", "sigma"),
                          stop("Unsupported endpoint type.")
  )

  # --- Filter out arms where any required parameter is NA or n is not positive ---
  param_list_filtered <- Filter(function(arm) {
    if (is.null(arm) || length(arm) == 0) return(FALSE)

    all_required_present_and_not_na <- all(required_keys %in% names(arm)) &&
      !any(sapply(arm[required_keys], is.na))

    n_val <- arm$n
    is_n_positive_numeric <- is.numeric(n_val) && length(n_val) == 1 && n_val > 0 && is.finite(n_val)

    return(all_required_present_and_not_na && is_n_positive_numeric)

  }, param_list_raw)

  return(param_list_filtered)
}


