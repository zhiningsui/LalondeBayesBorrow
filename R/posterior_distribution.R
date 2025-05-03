#' Obtain Posterior Distribution for Continuous or Binary Endpoints
#'
#' This function calculates the posterior distribution for one-arm trial data with either continuous or binary endpoints, using current and optionally historical data.
#' It supports the SAM prior and allows for optional borrowing from historical data.
#'
#' @param endpoint A string. The type of endpoint, either "continuous" or "binary".
#' @param current A named vector. Current trial data:
#'   - For a continuous endpoint, the vector should include:
#'     - `n`: Number of patients.
#'     - `mu_hat`: Mean of the normal sampling distribution.
#'     - `s`: Standard deviation of the normal sampling distribution.
#'   - For a binary endpoint, the vector should include:
#'     - `n`: Number of patients.
#'     - `count`: Number of responses.
#' @param historical A named vector. Historical trial data in the same format as the `current` vector. If `NULL`, borrowing is not performed.
#' @param delta A scalar. The clinically significant difference for the SAM prior, used as the threshold for detecting prior-data conflict (default is `log(0.85)` for continuous endpoints).
#' @param w A scalar. Initial weight assigned to the informative component of the mixture prior. Set to `0` if no borrowing. If `NULL`, the SAM prior is used for automatic conflict detection.
#' @param a A scalar. The alpha parameter of the beta prior for binary endpoints (used in the mixture prior if borrowing).
#' @param b A scalar. The beta parameter of the beta prior for binary endpoints (used in the mixture prior if borrowing).
#'
#' @return A list containing the posterior distribution parameters:
#'   - For continuous endpoints:
#'     - `w_prior`: Weight of the prior after conflict detection (or user-defined weight).
#'     - `post`: A list with:
#'       - `w`: Posterior weight of the informative component.
#'       - `theta1`, `sigma1`: Posterior mean and standard deviation from informative prior.
#'       - `theta2`, `sigma2`: Posterior mean and standard deviation from non-informative prior.
#'   - For binary endpoints:
#'     - `w_prior`: Weight of the prior after conflict detection (or user-defined weight).
#'     - `post`: A list with:
#'       - `w`: Posterior weight of the informative component.
#'       - `a1`, `b1`: Parameters of the informative beta posterior.
#'       - `a2`, `b2`: Parameters of the non-informative beta posterior.
#' @export
posterior_distribution <- function(endpoint, current, historical = NULL, delta = log(0.85),
                                   w = NULL, a = 0.01, b = 0.01, a0 = 0.01, b0 = 0.01,
                                   theta0 = 0, s0 = 100, ess_h = NULL) {
  # w.input record the initial weight; NULL = SAM prior, 0 = No borrowing
  w.input <- w

  # if no historical data is provided, automatically set w.input to 0 to indicate no borrowing.
  if (is.null(historical)) {
    if (is.null(w.input) | w.input != 0 ) stop("posterior_distribution(): No historical data provided. Cannot do borrowing with w not equal to zero.")
    w.input <- 0
    historical <- current
  }

  if (endpoint == "continuous") {
    if(!("n" %in% names(historical) & "mu_hat" %in% names(historical) & "s" %in% names(historical) &
         "n" %in% names(current) & "mu_hat" %in% names(current) & "s" %in% names(current))) {
      stop("posterior_distribution(): Current data and/or historical data do not have n, mu_hat, or s column.")
    }

    nh <- historical["n"]
    yh <- historical["mu_hat"]
    sh <- historical["s"]

    n <- current["n"]
    y <- current["mu_hat"]
    s <- current["s"]

    # Obtain discounting parameter for the power prior
    if (!is.null(ess_h)) {
      power_alpha <- ess_h/nh
    } else {
      power_alpha <- 1
    }

    if (!is.null(w.input) && w.input == 0) { # if no borrowing (either because of no historical data, or because of large conflict)
      theta2 <- (s^2 * theta0 + s0^2 * y) / (s^2 + s0^2)
      sig2 <- sqrt((s0^2 * s^2) / (s0^2 + s^2))

      post = list(
        w = 0,
        mu1 = NA, sigma1 = NA,
        mu2 = unname(theta2), sigma2 = unname(sig2)
      )
    } else {
      if (is.null(w.input)) {
        s_pooled <- sqrt( ( (nh-1)*sh^2 + (n-1)*s^2 ) / (nh+n-2) )
        R <- exp(-((y - yh) / s_pooled)^2 / 2 - max(-((y - yh - delta) / s_pooled)^2 / 2, -((y - yh + delta) / s_pooled)^2 / 2))
        w <- R / (1 + R)
      }

      sh_discount <- sqrt( sh^2 / power_alpha )

      theta1 <- (s^2 * yh + sh_discount^2 * y ) / (s^2 + sh_discount^2 )
      sig1 <- sqrt( (s^2 * sh_discount^2) / (s^2 + sh_discount^2) )

      theta2 <- (s^2 * theta0 + s0^2 * y) / (s^2 + s0^2)
      sig2 <- sqrt(s^2 * s0^2 / (s^2 + s0^2))

      # theta2 <- (s^2 * yh + sh^2 * y * nh) / (s^2 + sh^2 * nh)
      # sig2 <- sqrt(s^2 * sh^2 * nh / (s^2 + sh^2 * nh))

      part1 <- 1 / sqrt(s^2 + sh_discount^2) * exp(-0.5 * (y - yh)^2 / (s^2 + sh_discount^2) )
      part2 <- 1 / sqrt(s^2 + s0^2) * exp(-0.5 * (y - theta0)^2 / (s^2 + s0^2) )

      ws <- w * part1 / (w * part1 + (1 - w) * part2)

      post = list(
        w = unname(ws),
        mu1 = unname(theta1), sigma1 = unname(sig1),
        mu2 = unname(theta2), sigma2 = unname(sig2)
      )
    }

    return(list(w_prior = unname(w), post = post))

  } else if (endpoint == "binary") {
    if(!("n" %in% names(historical) & "count" %in% names(historical) &
         "n" %in% names(current) & "count" %in% names(current))) {
      stop("posterior_distribution(): Current data and/or historical data do not have n or count column.")
    }

    # if no borrowing, assign historical to NULL, i.e., no historical data
    if(!is.null(w.input)){
      if(w.input == 0) {
        historical = NULL
      }
    }

    x <- as.numeric(current["count"])
    n <- as.numeric(current["n"])
    xh <- as.numeric(historical["count"])
    nh <- as.numeric(historical["n"])

    # Obtain discounting parameter for the power prior
    if (!is.null(ess_h)) {
      power_alpha <- ess_h/nh
    } else {
      power_alpha <- 1
    }

    a2 <- a + x
    b2 <- b + n - x

    if (!is.null(w.input) && w.input == 0) {
      post = list(
        w = 0,
        a1 = NA, b1 = NA,
        a2 = unname(a2), b2 = unname(b2)
      )
    } else {
      if (is.null(w.input)) {
        thetah <- (a + xh) / (a + b + nh)
        R <- thetah^x * (1 - thetah)^(n - x) / max((thetah + delta)^x * (1 - thetah - delta)^(n - x), (thetah - delta)^x * (1 - thetah + delta)^(n - x))
        w <- R / (1 + R)
      }
      a1 <- a0 + power_alpha * xh + x
      b1 <- b0 + power_alpha * (nh - xh) + n - x
      z0 <- beta(a2, b2) / beta(a, b)
      z1 <- beta(a1, b1) / beta(a0 + power_alpha * xh, b0 + power_alpha * (nh - xh))
      ws <- w * z1 / (w * z1 + (1 - w) * z0)

      post = list(
        w = unname(ws),
        a1 = unname(a1), b1 = unname(b1),
        a2 = unname(a2), b2 = unname(b2)
      )
    }
    return(list(w_prior = unname(w), post = post))
  }
}
