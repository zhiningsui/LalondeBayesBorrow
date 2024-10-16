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
#' @param gate A scalar.
#'
#' @return A list containing the posterior distribution parameters:
#'   - For continuous endpoints:
#'     - `w_prior`: Weight of the prior after conflict detection (or user-defined weight).
#'     - `post`: A list with:
#'       - `w`: Posterior weight of the informative component.
#'       - `mu1`, `sigma1`: Posterior mean and standard deviation from informative prior.
#'       - `mu2`, `sigma2`: Posterior mean and standard deviation from non-informative prior.
#'   - For binary endpoints:
#'     - `w_prior`: Weight of the prior after conflict detection (or user-defined weight).
#'     - `post`: A list with:
#'       - `w`: Posterior weight of the informative component.
#'       - `a1`, `b1`: Parameters of the informative beta posterior.
#'       - `a2`, `b2`: Parameters of the non-informative beta posterior.
#' @export
#'
#' @examples
#' # Example with a continuous endpoint:
#' current_data <- c(n = 100, mu_hat = 2.5, s = 0.5)
#' historical_data <- c(n = 80, mu_hat = 2.0, s = 0.6)
#' posterior_distribution(endpoint = "continuous", current = current_data,
#'                        historical = historical_data, delta = log(0.85),
#'                        w = NULL)
#'
#' # Example with a binary endpoint:
#' current_data <- c(n = 100, count = 45)
#' historical_data <- c(n = 90, count = 40)
#' posterior_distribution(endpoint = "binary", current = current_data,
#'                        historical = historical_data, a = 1, b = 1, w = NULL)
posterior_distribution <- function(endpoint, current, historical = NULL, delta = log(0.85), w = NULL, a = NULL, b = NULL, gate = NULL) {
  # w.input record the initial weight; NULL = SAM prior, 0 = No borrowing
  w.input <- w

  # if no historical data is provided, automatically set w.input to 0 to indicate no borrowing.
  if (is.null(historical)){
    if (is.null(w.input) | w.input != 0 ) stop("posterior_distribution(): No historical data provided. Cannot do borrowing with w not equal to zero.")
    w.input <- 0
    historical <- current
  }

  if (endpoint == "continuous") {

    if(!("n" %in% names(historical) & "mu_hat" %in% names(historical) & "s" %in% names(historical) &
         "n" %in% names(current) & "mu_hat" %in% names(current) & "s" %in% names(current))) {
      stop("posterior_distribution(): Current data and/or historical data do not have n, mu_hat, or s column.")
    }

    # if no borrowing, assign historical to NULL, i.e., no historical data
    if(!is.null(w.input)){
      if(w.input == 0) {
        historical = NULL
      }
    }

    nh <- historical["n"]
    y <- current["mu_hat"]
    s <- current["s"]
    yh <- historical["mu_hat"]
    sh <- historical["s"]

    # if gate is given, consider gate value before any borrowing
    if (!is.null(gate)) {
      if (is.null(historical)) {stop("posterior_distribution(): No historical data provided. Cannot compare prior-data conflict with the gate.")}
      if (y-yh > gate | yh-y > gate){
        w.input = 0
        w = 0
      }
    }

    if (!is.null(w.input) && w.input == 0) { # if no borrowing (either because of no historical data, or because of large conflict)
      if (is.null(historical)){ # if no historical data is provided, using N(0,100^2)
        mu2 <- (100^2 * y) / (s^2 + 100^2)
        sig2 <- sqrt((100^2 * s^2) / (100^2 + s^2))
      } else { # if historical data is provided
        mu2 <- (s^2 * yh + sh^2 * y * nh) / (s^2 + sh^2 * nh)
        sig2 <- sqrt(s^2 * sh^2 * nh / (s^2 + sh^2 * nh))
      }
      post = list(
        w = 0,
        mu1 = NA, sigma1 = NA,
        mu2 = unname(mu2), sigma2 = unname(sig2)
      )
    } else {
      if (is.null(w.input)) {
        R <- exp(-((y - yh) / s)^2 / 2 - max(-((y - yh - delta) / s)^2 / 2, -((y - yh + delta) / s)^2 / 2))
        w <- R / (1 + R)
      }
      mu1 <- (s^2 * yh + sh^2 * y) / (s^2 + sh^2)
      sig1 <- sqrt(s^2 * sh^2 / (s^2 + sh^2))
      mu2 <- (s^2 * yh + sh^2 * y * nh) / (s^2 + sh^2 * nh)
      sig2 <- sqrt(s^2 * sh^2 * nh / (s^2 + sh^2 * nh))
      part1 <- 1 / sqrt(s^2 + sh^2) * exp(-0.5 * (y - yh)^2 / (s^2 + sh^2) )
      part2 <- 1 / sqrt(s^2 + nh * sh^2) * exp(-0.5 * (y - yh)^2 / (s^2 + nh * sh^2) )
      # part1 <- sig1 / (s * sh) * exp(0.5 * (mu1^2 / sig1^2 - y^2 / s^2 - yh^2 / sh^2))
      # part2 <- sig2 / (s * sh / sqrt(1 / nh)) * exp(0.5 * (mu2^2 / sig2^2 - y^2 / s^2 - yh^2 / (sh^2 * nh)))
      ws <- w * part1 / (w * part1 + (1 - w) * part2)
      post = list(
        w = unname(ws),
        mu1 = unname(mu1), sigma1 = unname(sig1),
        mu2 = unname(mu2), sigma2 = unname(sig2)
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

    if (!is.null(gate)) {
      if (is.null(historical)) {stop("posterior_distribution(): No historical data provided. Cannot compare prior-data conflict with the gate.")}
      if (x/n-xh/nh > gate | xh/nh-x/n > gate){
        w.input = 0
        w = 0
      }
    }

    a2 <- a + x
    b2 <- b + n - x

    if (exists("w.input") && !is.null(w.input) && w.input == 0) {
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
      a1 <- a + xh + x
      b1 <- b + nh + n - xh - x
      z0 <- beta(a2, b2) / beta(a, b)
      z1 <- beta(a1, b1) / beta(a + xh, b + nh - xh)
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
