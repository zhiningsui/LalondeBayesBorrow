#' Obtain Posterior Distribution for Single Arm
#'
#' This function calculates the posterior distribution for one-arm trial data
#' with either continuous or binary endpoints, incorporating current data
#' and optionally borrowing from historical data. It supports the SAM prior
#' (Self-Adaptive Mixture) for automatic conflict detection when borrowing
#' is enabled (`w = NULL`) and allows specification of a fixed mixture weight
#' for the informative component (`w` is numeric). A power prior mechanism is
#' included via the `ess_h` parameter to control the influence of historical data.
#' Additionally, a gating mechanism is applied: if the prior-data conflict exceeds
#' a threshold (`delta_gate`), borrowing is disabled regardless of the SAM weight.
#'
#' @param endpoint A string. The type of endpoint, either `'continuous'`
#'   (assumes Normal likelihood and Normal-Inverse-Gamma conjugacy on mean/variance)
#'   or `'binary'` (assumes Binomial likelihood and Beta conjugacy).
#' @param current A named list or vector containing current trial data:
#'   * For `'continuous'`: must include `n` (number of patients), `mu_hat` (sample mean), `s` (sample standard deviation the mean).
#'   * For `'binary'`: must include `n` (number of patients), `count` (number of responses/events).
#'   Sample size `n` must be positive.
#' @param historical A named list or vector containing historical trial data
#'   in the same format as the `current` parameter. If `NULL` (default),
#'   no historical data is used, and borrowing is disabled (`w` is effectively 0).
#'   If provided, historical sample size `n` must be non-negative.
#' @param delta_SAM A scalar positive value, or `NULL`. The clinical significant difference
#'   threshold used in the SAM prior calculation to detect prior-data conflict.
#'   Only used if `w` is `NULL`. The interpretation depends on the `endpoint`:
#'   for `'continuous'`, it's on the mean scale; for `'binary'`, it's on the
#'   probability scale. If `NULL` and `w` is `NULL`, default values (0.2 for
#'   continuous, 0.1 for binary) are used.
#' @param delta_gate A scalar positive value or NULL. Threshold for the gating mechanism: if prior-data conflict exceeds this, borrowing is disabled. If NULL, gate is not applied.
#' @param w A scalar value between 0 and 1, or `NULL`.
#'   * If a numeric value (0 to 1): Specifies a fixed initial weight for
#'       the informative (historical) component of the mixture prior. Borrowing occurs
#'       according to this fixed weight.
#'   * If `NULL` (default): The SAM prior is used. The initial weight for
#'       the informative component is automatically calculated based on the
#'       agreement between historical and current data, using `delta`.
#'   Note: If `historical` is `NULL` or `ess_h` is 0, `w` is effectively forced to 0,
#'   and no borrowing occurs regardless of the input `w`.
#' @param a A scalar positive value. The alpha parameter for the *non-informative*
#'   Beta prior component (used only for `endpoint = "binary"`). Defaults to 0.01.
#' @param b A scalar positive value. The beta parameter for the *non-informative*
#'   Beta prior component (used only for `endpoint = "binary"`). Defaults to 0.01.
#' @param a0 A scalar positive value. The base alpha parameter for the *informative*
#'   Beta prior component (used only for `endpoint = "binary"`). This parameter
#'   (along with `b0`) serves as a base prior which is updated by the discounted
#'   historical data to form the informative prior component parameters. Defaults to 0.01.
#' @param b0 A scalar positive value. The base beta parameter for the *informative*
#'   Beta prior component (used only for `endpoint = "binary"`). Defaults to 0.01.
#' @param theta0 A scalar. The mean parameter for the *non-informative* Normal prior
#'   component (used only for `endpoint = "continuous"`). Defaults to 0.
#' @param s0 A scalar positive value. The standard deviation parameter for the
#'   *non-informative* Normal prior component (used only for `endpoint = "continuous"`).
#'   Defaults to 100.
#' @param ess_h A scalar non-negative finite value, or `NULL`. The Effective Sample Size
#'   to use for the historical data. This implements a power prior approach, where
#'   the historical data likelihood is raised to the power of `ess_h / nh`.
#'   If `NULL` (default), the power is 1 (full historical data influence).
#'   If 0, historical data is fully discounted, equivalent to no borrowing.
#'
#' @return A list containing the results of the posterior calculation:
#'   * `w_prior`: The initial mixture weight used for the informative prior
#'       component *before* updating with current data. This is either the input `w`
#'       or the weight calculated by the SAM prior if `w` was `NULL`. Returns 0
#'       if no historical data was provided, historical `n` was 0, or `ess_h` was 0.
#'   * `post`: A named list containing the parameters of the resulting
#'       two-component mixture posterior distribution, suitable for use with RBesT:
#'       * `w`: The *posterior* weight of the informative component. This is calculated
#'           by updating the `w_prior` with the current data via Likelihood Ratio Test.
#'       * For continuous (`mixnorm` format): `mu1`, `sigma1` (mean and SD of
#'           the informative posterior component), `mu2`, `sigma2` (mean and SD of
#'           the non-informative posterior component). These will be `NA_real_` if the
#'           corresponding posterior component has zero weight (`post$w` is 0 or 1).
#'       * For binary (`mixbeta` format): `a1`, `b1` (parameters of the
#'           informative Beta posterior component), `a2`, `b2` (parameters of the
#'           non-informative Beta posterior component). These will be `NA_real_` if the
#'           corresponding posterior component has zero weight (`post$w` is 0 or 1).
#'   Returns an error if inputs are invalid or required parameters are missing.
#' @importFrom stats dnorm dbinom
#' @importFrom utils hasName
#' @export
#'
#' @examples
#' # Continuous Endpoint - SAM prior
#' current_cont <- list(n = 50, mu_hat = 1.8, s = 0.5)
#' historical_cont <- list(n = 100, mu_hat = 1.5, s = 0.4)
#' # delta = 0.2 for continuous is on the mean scale
#' post_sam_cont <- posterior_distribution(endpoint = "continuous",
#'                                         current = current_cont,
#'                                         historical = historical_cont,
#'                                         delta_SAM = 0.2, delta_gate = 0.2,
#'                                         theta0 = 1.6, s0 = 2)
#' str(post_sam_cont)
#'
#' # Continuous Endpoint - Fixed weight (no SAM), with ESS discounting
#' post_fixed_cont_ess <- posterior_distribution(endpoint = "continuous",
#'                                              current = current_cont,
#'                                              historical = historical_cont,
#'                                              w = 0.6, ess_h = 40,
#'                                              theta0 = 1.6, s0 = 2)
#' str(post_fixed_cont_ess)
#'
#' # Continuous Endpoint - No historical data (w is ignored)
#' post_nocorr_cont <- posterior_distribution(endpoint = "continuous",
#'                                            current = current_cont,
#'                                            historical = NULL,
#'                                            theta0 = 1.6, s0 = 2)
#' str(post_nocorr_cont)
#'
#' # Binary Endpoint - SAM prior
#' current_bin <- list(n = 60, count = 45)
#' historical_bin <- list(n = 100, count = 60)
#' # delta = 0.1 for binary is on the probability scale (0 to 1)
#' post_sam_bin <- posterior_distribution(endpoint = "binary",
#'                                        current = current_bin,
#'                                        historical = historical_bin,
#'                                        delta_SAM = 0.1, delta_gate = 0.1,
#'                                        a = 0.5, b = 0.5,
#'                                        a0 = 1, b0 = 1)
#' str(post_sam_bin)
#'
#' # Binary Endpoint - Fixed weight with ESS discounting
#' post_fixed_ess_bin <- posterior_distribution(endpoint = "binary",
#'                                              current = current_bin,
#'                                              historical = historical_bin,
#'                                              w = 0.8, ess_h = 30,
#'                                              a = 0.5, b = 0.5, a0 = 1, b0 = 1)
#' str(post_fixed_ess_bin)
#'
#' # Binary Endpoint - No borrowing (w=0 explicitly)
#' post_noborrow_bin <- posterior_distribution(endpoint = "binary",
#'                                             current = current_bin,
#'                                             historical = historical_bin,
#'                                             w = 0, a = 0.5, b = 0.5)
#' str(post_noborrow_bin)
#'
#' @seealso \code{\link{posterior_inference}}, \code{\link{posterior_prob}},
#'   \code{\link[RBesT]{mix}}, \code{\link[RBesT]{mixnorm}}, \code{\link[RBesT]{mixbeta}}
posterior_distribution <- function(endpoint, current, historical = NULL, delta_SAM = NULL,
                                   delta_gate = NULL,
                                   w = NULL, a = 0.01, b = 0.01, a0 = 0.01, b0 = 0.01,
                                   theta0 = 0, s0 = 100, ess_h = NULL) {
  # --- Input Validation ---
  supported_endpoints <- c("continuous", "binary")
  if (!is.character(endpoint) || length(endpoint) != 1 || !(endpoint %in% supported_endpoints)) {
    stop(paste("Input 'endpoint' must be a single string, either 'continuous' or 'binary'. Supported endpoints are:", paste(supported_endpoints, collapse = ", ")))
  }
  if (!is.list(current) && !is.vector(current)) {
    stop("Input 'current' must be a named list or vector.")
  }
  if (!is.null(historical) && !is.list(historical) && !is.vector(historical)) {
    stop("Input 'historical' must be NULL or a named list/vector.")
  }

  # Check required parameters for current and historical data structure
  if (endpoint == "continuous") {
    req_current_names <- req_historical_names <- c("n", "mu_hat", "s")
  } else { # binary
    req_current_names <- req_historical_names <- c("n", "count")
  }

  if (!all(hasName(current, req_current_names))) {
    stop(paste("Current data is missing required elements for", endpoint, "endpoint:", paste(setdiff(req_current_names, names(current)), collapse = ", ")))
  }
  if (!is.null(historical) && !all(hasName(historical, req_historical_names))) {
    stop(paste("Historical data is provided but is missing required elements for", endpoint, "endpoint:", paste(setdiff(req_historical_names, names(historical)), collapse = ", ")))
  }

  # --- Extract and Validate Data Values ---
  extract_and_validate <- function(data, name, param_name, is_positive = FALSE, is_non_negative = FALSE, is_proportion = FALSE, is_integer = FALSE) {
    val <- data[[param_name]]
    if (is.null(val)) {
      stop(paste("Internal error: Parameter '", param_name, "' missing from ", name, " data after initial check.", sep=""))
    }
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
      stop(paste("Parameter '", param_name, "' in ", name, " data must be a single finite numeric value.", sep=""))
    }
    if (is_positive && val <= 0) {
      stop(paste("Parameter '", param_name, "' in ", name, " data must be positive.", sep=""))
    }
    if (is_non_negative && val < 0) {
      stop(paste("Parameter '", param_name, "' in ", name, " data must be non-negative.", sep=""))
    }
    if (is_proportion && (val < 0 || val > 1)) {
      stop(paste("Parameter '", param_name, "' in ", name, " data must be between 0 and 1.", sep=""))
    }
    if (is_integer && val >= 0 && val != floor(val)) {
      stop(paste("Parameter '", param_name, "' in ", name, " data must be an integer.", sep=""))
    }
    return(as.numeric(val))
  }

  # Extract current data
  n <- extract_and_validate(current, "current", "n", is_positive = TRUE, is_integer = TRUE)
  if (endpoint == "continuous") {
    y <- extract_and_validate(current, "current", "mu_hat")
    s <- extract_and_validate(current, "current", "s", is_positive = TRUE)
  } else { # binary
    x <- extract_and_validate(current, "current", "count", is_non_negative = TRUE, is_integer = TRUE)
    if (x > n) {
      stop("Current count cannot be greater than current sample size 'n'.")
    }
  }

  # Extract historical data (or set defaults if NULL)
  if (!is.null(historical)) {
    nh <- extract_and_validate(historical, "historical", "n", is_non_negative = TRUE, is_integer = TRUE)
    if (endpoint == "continuous") {
      yh <- extract_and_validate(historical, "historical", "mu_hat")
      sh <- extract_and_validate(historical, "historical", "s", is_positive = TRUE)
    } else { # binary
      xh <- extract_and_validate(historical, "historical", "count", is_non_negative = TRUE, is_integer = TRUE)
      if (!is.na(xh) && !is.na(nh) && xh > nh) { # Check only if both are not NA
        stop("Historical count cannot be greater than historical sample size 'n'.")
      }
    }
  } else {
    nh <- 0
    if (endpoint == "continuous") {
      yh <- NA_real_
      sh <- NA_real_
    } else { # binary
      xh <- NA_real_
    }
  }

  # Validate and extract prior parameters (using default if necessary)
  extract_prior_param <- function(param, name, is_positive = FALSE, is_non_negative = FALSE) {
    if (is.null(param)) {
      stop(paste("Internal error: Prior parameter '", name, "' is NULL.", sep=""))
    }
    if (!is.numeric(param) || length(param) != 1 || !is.finite(param)) {
      stop(paste("Prior parameter '", name, "' must be a single finite numeric value.", sep=""))
    }
    if (is_positive && param <= 0) {
      stop(paste("Prior parameter '", name, "' must be positive.", sep=""))
    }
    if (is_non_negative && param < 0) {
      stop(paste("Prior parameter '", name, "' must be non-negative.", sep=""))
    }
    return(as.numeric(param))
  }

  if (endpoint == "continuous") {
    theta0 <- extract_prior_param(theta0, "theta0")
    s0 <- extract_prior_param(s0, "s0", is_positive = TRUE)
  } else { # binary
    a <- extract_prior_param(a, "a", is_positive = TRUE)
    b <- extract_prior_param(b, "b", is_positive = TRUE)
    a0 <- extract_prior_param(a0, "a0", is_positive = TRUE)
    b0 <- extract_prior_param(b0, "b0", is_positive = TRUE)
  }

  # Determine effective historical sample size (ess_h)
  if (is.null(ess_h)) {
    ess_h_used <- nh # Default ess_h is historical n
  } else {
    if (!is.numeric(ess_h) || length(ess_h) != 1 || !is.finite(ess_h) || ess_h < 0) {
      stop("Parameter 'ess_h' must be a single non-negative finite numeric value when not NULL.")
    }
    ess_h_used <- ess_h
  }

  # If effective historical data size is 0, force no borrowing regardless of w
  if (ess_h_used < .Machine$double.eps) { # Treat near-zero ess_h as zero
    ess_h_used <- 0
    if (!is.null(w) && w > 0) {
      warning("Effective historical sample size (ess_h) is 0, forcing borrowing weight w to 0.")
    } else if (is.null(w)) {
      warning("Effective historical sample size (ess_h) is 0, forcing SAM prior weight w to 0.")
    }
    w_input_used <- 0
  } else {
    # Validate fixed w if provided
    if (!is.null(w)) {
      if (!is.numeric(w) || length(w) != 1 || !is.finite(w) || w < 0 || w > 1) {
        stop("Parameter 'w' must be a single numeric value between 0 and 1 when not NULL.")
      }
    }
    w_input_used <- w # Keep original w input if ess_h > 0
  }

  # Calculate power prior exponent (alpha)
  power_alpha <- if (nh > 0) ess_h_used / nh else 0

  # --- Determine Prior Mixture Weight (w_prior) ---
  w_prior <- w_input_used # Start with input w (could be NULL or fixed or forced to 0)

  if (is.null(w_input_used) && power_alpha > 0) {
    if (is.null(delta_SAM) || !is.numeric(delta_SAM) || length(delta_SAM) != 1 || !is.finite(delta_SAM) || delta_SAM <= 0) {
      delta_SAM <- if (endpoint == "continuous") 0.2 else 0.1
      warning(paste("Parameter 'delta_SAM' was NULL or invalid. Using default delta_SAM =", delta_SAM, "for", endpoint, "endpoint."))
    }

    if (endpoint == "continuous") {
      if (is.na(sh) || sh <= 0)  stop("Historical standard deviation 's' must be positive for SAM calculation with continuous endpoint.")

      # Pooled standard deviation
      pooled_n_total <- nh + n
      df_pooled <- if (pooled_n_total > 2) pooled_n_total - 2 else 1 # Avoid negative or zero df
      sigma2_pooled <- ( ifelse(nh > 1, nh - 1, 0) * nh * sh^2 + ifelse(n > 1, n - 1, 0) * n * s^2 ) / df_pooled
      s_pooled <- sqrt( sigma2_pooled / n )

      if (!is.null(delta_gate) && abs(y - yh) >= delta_gate) {
        message("Gate applied: |current - historical| >= delta_gate. Setting w_prior = 0.")
        w_prior <- 0
      } else {
        if (!is.null(delta_gate)) {
          message("Gate not triggered: |current - historical| < delta_gate. Proceeding with SAM prior.")
        }

        # SAM starts here

        if (s_pooled < .Machine$double.eps) { # Avoid division by zero if s_pooled is 0
          warning("Pooled SD ~0. Unstable SAM weight. Forcing w_prior = 0.")
          w_prior <- 0
        } else {
          # SAM conflict calculation for Normal means based on likelihood ratio of shifted means
          # Likelihood at observed difference (y-yh)
          loglik_obs_diff <- dnorm(y - yh, mean = 0, sd = s_pooled, log = TRUE)

          # Maximum likelihood at shifted difference (delta or -delta)
          loglik_shifted_plus <- dnorm(y - yh, mean = delta_SAM, sd = s_pooled, log = TRUE)
          loglik_shifted_minus <- dnorm(y - yh, mean = -delta_SAM, sd = s_pooled, log = TRUE)

          loglik_max_shifted <- max(loglik_shifted_plus, loglik_shifted_minus)
          log_R <- loglik_obs_diff - loglik_max_shifted

          R <- if (log_R > 700) 1e300 else if (log_R < -700) 1e-300 else exp(log_R)
          w_prior <- R / (1 + R)
        }
      }
    } else if (endpoint == "binary") {
      if (is.na(xh) || is.na(nh) || nh <= 0) {
        stop("Internal error: Historical data invalid for SAM calculation.")
      }

      thetah_prior_mean <- (a + xh) / (a + b + nh)
      thetac_current_mean <- x / n

      if (!is.null(delta_gate) && abs(thetac_current_mean - thetah_prior_mean) >= delta_gate) {
        message("Gate applied: |current - historical| >= delta_gate. Setting w_prior = 0.")
        w_prior <- 0
      } else {  # Otherwise, compute SAM prior weight as usual
        if (!is.null(delta_gate)) {
          message("Gate not triggered: |current - historical| < delta_gate. Proceeding with SAM prior.")
        }

        if (thetah_prior_mean <= .Machine$double.eps || thetah_prior_mean >= 1 - .Machine$double.eps) {
          warning("Historical prior mean near 0/1. Unstable SAM weight. Forcing w_prior = 0.")
          w_prior <- 0
        } else {
          loglik_obs <- dbinom(x, size = n, prob = thetah_prior_mean, log = TRUE)
          # Calculate likelihood of current data under shifted probabilities (+/- delta)
          thetah_plus_delta <- min(thetah_prior_mean + delta_SAM, 1 - .Machine$double.eps) # Subtract tiny epsilon to avoid prob=1 edge case in dbinom
          thetah_minus_delta <- max(thetah_prior_mean - delta_SAM, .Machine$double.eps) # Add tiny epsilon to avoid prob=0 edge case in dbinom

          loglik_shifted_plus <- dbinom(x, size = n, prob = thetah_plus_delta, log = TRUE)
          loglik_shifted_minus <- dbinom(x, size = n, prob = thetah_minus_delta, log = TRUE)
          loglik_max_shifted <- max(loglik_shifted_plus, loglik_shifted_minus)

          log_R <- loglik_obs - loglik_max_shifted
          R <- if (log_R > 700) 1e300 else if (log_R < -700) 1e-300 else exp(log_R)
          w_prior <- R / (1 + R)
        }
      }
    }
  } else if (is.null(w_input_used) && power_alpha == 0) {
    w_prior <- 0 # If w is NULL but effective historical data is zero, force w_prior to 0
  }
  w_prior <- max(0, min(1, w_prior))

  # --- Calculate Posterior Distribution Parameters ---
  post <- list()
  ws <- 0

  if (w_prior < .Machine$double.eps) { # Treat near-zero prior weight as zero
    w_prior <- 0
    ws <- 0

    if (endpoint == "continuous") {
      theta2 <- (s^2 * theta0 + s0^2 * y) / (s^2 + s0^2)
      sig2 <- sqrt((s0^2 * s^2) / (s0^2 + s^2))

      post = list(
        w = 0,
        mu1 = NA_real_, sigma1 = NA_real_,
        mu2 = unname(theta2), sigma2 = unname(sig2)
      )
    } else { # binary
      a2 <- a + x
      b2 <- b + n - x

      post = list(
        w = 0,
        a1 = NA_real_, b1 = NA_real_,
        a2 = unname(a2), b2 = unname(b2)
      )
    }

  } else if (w_prior > 1 - .Machine$double.eps) { # Treat near-one prior weight as one
    w_prior <- 1
    ws <- 1

    if (endpoint == "continuous") {
      sh_discount_var <- sh^2 / power_alpha
      theta1 <- (s^2 * yh + sh_discount_var * y ) / (s^2 + sh_discount_var )
      sig1 <- sqrt( (s^2 * sh_discount_var) / (s^2 + sh_discount_var) )

      post = list(
        w = 1,
        mu1 = unname(theta1), sigma1 = unname(sig1),
        mu2 = NA_real_, sigma2 = NA_real_
      )

    } else { # binary
      a1 <- a0 + power_alpha * xh + x
      b1 <- b0 + power_alpha * (nh - xh) + n - x

      post = list(
        w = 1,
        a1 = unname(a1), b1 = unname(b1),
        a2 = NA_real_, b2 = NA_real_
      )
    }

  }
  else {  # Case 3: Borrowing occurs with prior weight > 0 and < 1 (either fixed or SAM)

    if (endpoint == "continuous") {
      if (is.na(sh) || sh <= 0) {
        stop("Historical standard deviation 's' must be positive to calculate informative posterior component.")
      }

      sh_discount_var <- sh^2 / power_alpha
      theta1 <- (s^2 * yh + sh_discount_var * y ) / (s^2 + sh_discount_var )
      sig1 <- sqrt( (s^2 * sh_discount_var) / (s^2 + sh_discount_var) )

      theta2 <- (s^2 * theta0 + s0^2 * y) / (s^2 + s0^2)
      sig2 <- sqrt(s^2 * s0^2 / (s^2 + s0^2))

      log_part1_lik <- dnorm(y, mean = yh, sd = sqrt(sh_discount_var + s^2), log = TRUE)
      log_part2_lik <- dnorm(y, mean = theta0, sd = sqrt(s0^2 + s^2), log = TRUE)
      log_ratio <- (log(1 - w_prior) + log_part2_lik) - (log(w_prior) + log_part1_lik)

      if (log_ratio > 700) {
        ws <- 0
      } else if (log_ratio < -700) {
        ws <- 1
      } else {
        ws <- 1 / (1 + exp(log_ratio))
      }

      post = list(
        w = unname(ws),
        mu1 = unname(theta1), sigma1 = unname(sig1),
        mu2 = unname(theta2), sigma2 = unname(sig2)
      )

    } else { # binary
      if (is.na(xh) || is.na(nh) || nh <= 0) {
        stop("Internal error: Historical data invalid for informative posterior calculation.")
      }

      a1 <- a0 + power_alpha * xh + x
      b1 <- b0 + power_alpha * (nh - xh) + n - x

      a2 <- a + x
      b2 <- b + n - x

      a_h_prior <- a0 + power_alpha * xh
      b_h_prior <- b0 + power_alpha * (nh - xh)

      log_part1_lik <- lbeta(a1, b1) - lbeta(a_h_prior, b_h_prior)
      log_part2_lik <- lbeta(a2, b2) - lbeta(a, b)
      log_ratio <- (log(1 - w_prior) + log_part2_lik) - (log(w_prior) + log_part1_lik)

      if (log_ratio > 700) {
        ws <- 0
      } else if (log_ratio < -700) {
        ws <- 1
      } else {
        ws <- 1 / (1 + exp(log_ratio))
      }

      post = list(
        w = unname(ws),
        a1 = unname(a1), b1 = unname(b1),
        a2 = unname(a2), b2 = unname(b2)
      )
    }
  }

  # Ensure final posterior weight is within [0, 1]
  post$w <- max(0, min(1, post$w))

  # Set parameters of zero-weight component to NA_real_ for RBesT compatibility
  if (post$w < .Machine$double.eps) { # Treat near-zero posterior weight as zero
    post$w <- 0
    if (endpoint == "continuous") {
      post$mu1 <- NA_real_
      post$sigma1 <- NA_real_
    } else { # binary
      post$a1 <- NA_real_
      post$b1 <- NA_real_
    }
  } else if (post$w > 1 - .Machine$double.eps) { # Treat near-one posterior weight as one
    post$w <- 1
    if (endpoint == "continuous") {
      post$mu2 <- NA_real_
      post$sigma2 <- NA_real_
    } else { # binary
      post$a2 <- NA_real_
      post$b2 <- NA_real_
    }
  }

  return(list(w_prior = unname(w_prior), post = post))
}
