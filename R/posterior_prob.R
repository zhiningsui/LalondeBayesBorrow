#' Calculate Posterior Probability
#'
#' This function calculates the posterior probability of a parameter estimate (or the
#' difference between two parameter estimates) falling within a specified range,
#' based on provided RBesT mixture posterior distributions.
#' For a single posterior distribution, it computes the probability of the parameter
#' being greater than, less than, or between specified values. For two posterior
#' distributions, it computes the probability of the difference between the
#' first and second distributions falling within the range.
#'
#' @param range_type A character string. Specifies the type of comparison range.
#'   Must be one of `'greater'`, `'less'`, or `'between'`.
#'   Defaults to `'greater'`.
#' @param post1 An RBesT mixture distribution object (class `mix`, `mixnorm`, `mixbeta`, etc.).
#'   Represents the first posterior distribution.
#' @param post2 An RBesT mixture distribution object (class `mix`, `mixnorm`, `mixbeta`, etc.).
#'   Optional second posterior distribution object for pairwise comparison of `post1 - post2`.
#'   If `NULL` (default), the probability is calculated for `post1` itself.
#' @param value A scalar numeric value. The threshold used for comparison.
#'   If `range_type` is "greater" or "less", this is the single threshold.
#'   If `range_type` is "between", this is the lower bound of the interval.
#' @param value2 A scalar numeric value. The upper bound for the "between" comparison.
#'   This parameter is required if `range_type` is `"between"`.
#'
#' @return A scalar numeric value representing the calculated posterior probability,
#'   between 0 and 1. Returns `NA` if there is an error during the probability calculation
#'   (e.g., `pmixdiff` fails for two distributions).
#' @export
#' @import RBesT
#'
#' @examples
#' # Assuming RBesT is installed and loaded
#' library(RBesT) # Load RBesT for example
#' # Example 1: Single arm - Probability of being greater than a value
#' post_single_norm <- mixnorm(c(0.6, 10, 2), c(0.4, 15, 3))
#' prob_greater <- posterior_prob(post1 = post_single_norm, value = 12, range_type = "greater")
#' cat("P(post_single_norm > 12):", prob_greater, "\n")
#'
#' # Example 2: Single arm - Probability of being less than a value
#' prob_less <- posterior_prob(post1 = post_single_norm, value = 11, range_type = "less")
#' cat("P(post_single_norm < 11):", prob_less, "\n")
#'
#' # Example 3: Single arm - Probability of being between two values
#' prob_between <- posterior_prob(post1 = post_single_norm, value = 11,
#'                               range_type = "between", value2 = 14)
#' cat("P(11 < post_single_norm < 14):", prob_between, "\n")
#'
#' # Example 4: Two arms - Probability of difference > 0 (post1 > post2)
#' post_ctrl_norm <- mixnorm(c(1, 12, 2.5)) # Control arm posterior
#' prob_diff_greater_0 <- posterior_prob(post1 = post_single_norm, post2 = post_ctrl_norm,
#'                                      value = 0, range_type = "greater")
#' cat("P(post_single_norm - post_ctrl_norm > 0):", prob_diff_greater_0, "\n")
#'
#' # Example 5: Two arms - Probability of difference < a value
#' prob_diff_less_neg1 <- posterior_prob(post1 = post_single_norm, post2 = post_ctrl_norm,
#'                                      value = -1, range_type = "less")
#' cat("P(post_single_norm - post_ctrl_norm < -1):", prob_diff_less_neg1, "\n")
#'
#' # Example 6: Two arms - Probability of difference between two values
#' prob_diff_between <- posterior_prob(post1 = post_single_norm, post2 = post_ctrl_norm,
#'                                    value = -2, range_type = "between", value2 = 1)
#' cat("P(-2 < post_single_norm - post_ctrl_norm < 1):", prob_diff_between, "\n")
#'
#' @seealso \code{\link{posterior_distribution}}, \code{\link{posterior_inference}},
#'   \code{\link[RBesT]{pmix}}, \code{\link[RBesT]{pmixdiff}}
posterior_prob <- function(range_type = "greater", post1, post2 = NULL, value, value2 = NULL) {

  # --- Input Validation ---
  if (!inherits(post1, "mix")) {
    stop("Input 'post1' must be an RBesT mixture distribution object (class 'mix').")
  }
  if (!is.null(post2) && !inherits(post2, "mix")) {
    stop("Input 'post2' must be NULL or an RBesT mixture distribution object (class 'mix').")
  }
  if (!is.numeric(value) || length(value) != 1 || !is.finite(value)) {
    stop("Input 'value' must be a single finite numeric value.")
  }

  valid_range_types <- c("greater", "less", "between")
  if (!is.character(range_type) || length(range_type) != 1 || !(range_type %in% valid_range_types)) {
    stop(paste("Input 'range_type' must be one of:", paste(valid_range_types, collapse = ", ")))
  }

  if (range_type == "between") {
    if (is.null(value2)) {
      stop("Input 'value2' must be provided when range_type is 'between'.")
    }
    if (!is.numeric(value2) || length(value2) != 1 || !is.finite(value2)) {
      stop("Input 'value2' must be a single finite numeric value when range_type is 'between'.")
    }
    if (value >= value2) {
      stop("For range_type 'between', 'value' (lower bound) must be less than 'value2' (upper bound).")
    }
  } else {
    if (!is.null(value2)) {
      warning("Input 'value2' is ignored when range_type is not 'between'.")
      value2 <- NULL
    }
  }

  # --- Calculate Probability ---
  if (is.null(post2)) { # Single Arm Probability
    if (range_type == "greater") {
      return(pmix(post1, value, FALSE)) # Pr(post1 > value)
    } else if (range_type == "less") {
      return(pmix(post1, value))
    } else if (range_type == "between") {
      if (is.null(value2)) stop("value2 must be provided for range_type 'between'")
      return(pmix(post1, value2) - pmix(post1, value))
    }

    if (range_type == "greater") {
      prob <- pmix(post1, q = value, lower.tail = FALSE) # Pr(post1 > value)
    } else if (range_type == "less") {
      prob <- pmix(post1, q = value, lower.tail = TRUE) # P(post1 <= value)
    } else if (range_type == "between") { # P(value < post1 <= value2) = P(post1 <= value2) - P(post1 <= value)
      prob_le_value2 <- pmix(post1, q = value2, lower.tail = TRUE)
      prob_le_value <- pmix(post1, q = value, lower.tail = TRUE)
      prob <- prob_le_value2 - prob_le_value
    }
  } else { # Two Arms Probability
    if (range_type == "greater") { # Pr( post1 - post2 > value)
      prob <- tryCatch(
        pmixdiff(post1, post2, q = value, lower.tail = FALSE),
        error = function(e) {
          warning("Error calculating posterior probability for difference (range_type='greater'): ", e$message)
          NA_real_
        }
      )
    } else if (range_type == "less") {  # Pr( post1 - post2 <= value)
      prob <- tryCatch(
        pmixdiff(post1, post2, q = value, lower.tail = TRUE),
        error = function(e) {
          warning("Error calculating posterior probability for difference (range_type='less'): ", e$message)
          NA_real_
        }
      )
    } else if (range_type == "between") {  # Pr( value < post1 - post2 <= value2) = P(X <= value2) - P(X <= value)
      prob <- tryCatch(
        {
          prob_le_value2 <- pmixdiff(post1, post2, q = value2, lower.tail = TRUE)
          prob_le_value <- pmixdiff(post1, post2, q = value, lower.tail = TRUE)
          prob_le_value2 - prob_le_value
        },
        error = function(e) {
          warning("Error calculating posterior probability for difference (range_type='between'): ", e$message)
          NA_real_
        }
      )
    }
  }

  if (!is.na(prob)) {
    prob <- max(0, min(1, prob))
  }

  return(prob)
}
