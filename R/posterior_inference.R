#' Perform Posterior Inference
#'
#' This function calculates key summary statistics (posterior mean, standard deviation,
#' and credible intervals) for one or two RBesT mixture posterior distributions.
#' When two distributions are provided (e.g., for treatment and control arms),
#' it also computes these statistics for the difference or ratio between them.
#' The function supports exponentiating the results, useful when the original
#' parameter of interest is on a log scale (e.g., log-odds, log-mean ratio).
#'
#' @param post1 An RBesT mixture distribution object (class `mix`, `mixnorm`, `mixbeta`, etc.).
#'   Represents the posterior distribution for the first entity (e.g., treatment arm).
#' @param post2 An RBesT mixture distribution object (class `mix`, `mixnorm`, `mixbeta`, etc.), or `NULL`.
#'   Represents the posterior distribution for the second entity (e.g., control arm).
#'   If `NULL` (default), inference is performed only for `post1`.
#' @param quantiles A numeric vector of length 2 specifying the lower and upper
#'   quantiles for the credible interval (e.g., `c(0.025, 0.975)` for a 95% credible interval).
#'   Must be sorted and within (0, 1).
#' @param EXP_TRANSFORM A logical value. If `TRUE`, the calculated mean, credible interval bounds,
#'   and the comparison estimate are exponentiated (`exp(x)`). This is typically
#'   used when the mixture distribution models a parameter on the log scale (e.g., log-mean, log-odds).
#'   Note that the standard deviation (`sd1`, `sd2`, `sd_compare`) is always calculated
#'   on the original scale of the input mixture distributions.
#'
#' @return A data frame with exactly one row. Contains the following columns:
#'   * `est1`: Posterior mean for `post1`. Exponentiated if `EXP_TRANSFORM = TRUE`.
#'   * `sd1`: Posterior standard deviation for `post1`. Always on the original scale of `post1`.
#'   * `ci1_l`: Lower bound of the credible interval for `post1`. Exponentiated if `EXP_TRANSFORM = TRUE`.
#'   * `ci1_u`: Upper bound of the credible interval for `post1`. Exponentiated if `EXP_TRANSFORM = TRUE`.
#'   * `est2`: Posterior mean for `post2` (if `post2` is not `NULL`). Exponentiated if `EXP_TRANSFORM = TRUE`. `NA` if `post2` is `NULL`.
#'   * `sd2`: Posterior standard deviation for `post2` (if `post2` is not `NULL`). Always on the original scale of `post2`. `NA` if `post2` is `NULL`.
#'   * `ci2_l`: Lower bound of the credible interval for `post2` (if `post2` is not `NULL`). Exponentiated if `EXP_TRANSFORM = TRUE`. `NA` if `post2` is `NULL`.
#'   * `ci2_u`: Upper bound of the credible interval for `post2` (if `post2` is not `NULL`). Exponentiated if `EXP_TRANSFORM = TRUE`. `NA` if `post2` is `NULL`.
#'   * `est_compare`: Estimated treatment effect. If `EXP_TRANSFORM = FALSE`, this is `est1 - est2` (difference of means/parameters). If `EXP_TRANSFORM = TRUE`, this is `exp(est1) / exp(est2)` (ratio of exponentiated means/parameters). `NA` if `post2` is `NULL`.
#'   * `sd_compare`: Standard deviation of the *difference* (`post1 - post2`) on the original scale. Calculated as `sqrt(sd1^2 + sd2^2)`, assuming independence. This is *not* the standard deviation of the ratio if `EXP_TRANSFORM = TRUE`. `NA` if `post2` is `NULL`.
#'   * `compare_ci_l`: Lower bound of the credible interval for the treatment effect (`post1 - post2`). Exponentiated if `EXP_TRANSFORM = TRUE`. `NA` if `post2` is `NULL`.
#'   * `compare_ci_u`: Upper bound of the credible interval for the treatment effect (`post1 - post2`). Exponentiated if `EXP_TRANSFORM = TRUE`. `NA` if `post2` is `NULL`.
#'
#' @export
#' @import RBesT
#'
#' @examples
#' # Assuming RBesT is installed and loaded
#' library(RBesT) # Load RBesT for example
#' # Example 1: Single arm, continuous (Normal mixture)
#' post_single_norm <- mixnorm(c(0.6, 10, 2), c(0.4, 15, 3))
#' inference_single <- posterior_inference(post1 = post_single_norm, quantiles = c(0.025, 0.975))
#' print(inference_single)
#'
#' # Example 2: Two arms, continuous (Difference)
#' post_ctrl_norm <- mixnorm(c(1, 12, 2.5)) # Control arm posterior
#' inference_diff <- posterior_inference(post1 = post_single_norm, post2 = post_ctrl_norm,
#'                                     quantiles = c(0.025, 0.975))
#' print(inference_diff)
#'
#' # Example 3: Two arms, binary (Beta mixture), inference on probability scale
#' post_trt_beta <- mixbeta(c(0.7, 10, 5), c(0.3, 2, 8)) # Treatment posterior for probability
#' post_ctrl_beta <- mixbeta(c(1, 8, 12)) # Control posterior for probability
#' inference_prob_diff <- posterior_inference(post1 = post_trt_beta, post2 = post_ctrl_beta,
#'                                          quantiles = c(0.05, 0.95))
#' print(inference_prob_diff)
#'
#' # Example 4: Two arms, log-transformed data, inference on ratio scale
#' # Assume post_log_trt and post_log_ctrl are posteriors on the log scale
#' post_log_trt <- mixnorm(c(1, log(10), 0.5))
#' post_log_ctrl <- mixnorm(c(1, log(5), 0.6))
#' inference_log_ratio <- posterior_inference(post1 = post_log_trt, post2 = post_log_ctrl,
#'                                          quantiles = c(0.025, 0.975), EXP_TRANSFORM = TRUE)
#' # Note: est_compare is exp(log(10) - log(5)) = 10/5 = 2
#' # credible interval is exp(CI_log_diff) = exp(CI_log_trt - CI_log_ctrl) = CI_ratio
#' print(inference_log_ratio)
#'
#' @seealso \code{\link{posterior_distribution}}, \code{\link{posterior_prob}},
#'   \code{\link[RBesT]{summary.mix}}, \code{\link[RBesT]{qmix}}, \code{\link[RBesT]{qmixdiff}}
posterior_inference <- function(post1, post2 = NULL, quantiles, EXP_TRANSFORM = FALSE) {

  # --- Input Validation ---
  if (!inherits(post1, "mix")) {
    stop("Input 'post1' must be an RBesT mixture distribution object (class 'mix').")
  }
  if (!is.null(post2) && !inherits(post2, "mix")) {
    stop("Input 'post2' must be NULL or an RBesT mixture distribution object (class 'mix').")
  }
  if (!is.numeric(quantiles) || length(quantiles) != 2 || any(quantiles <= 0) || any(quantiles >= 1) || quantiles[1] >= quantiles[2]) {
    stop("Input 'quantiles' must be a numeric vector of length 2, sorted, and within (0, 1).")
  }
  if (!is.logical(EXP_TRANSFORM) || length(EXP_TRANSFORM) != 1) {
    stop("Input 'EXP_TRANSFORM' must be a single logical value (TRUE or FALSE).")
  }

  calculate_post_stats <- function(post) {
    summary_stats <- summary(post)
    est <- summary_stats["mean"]
    sd <- summary_stats["sd"]
    conf <- qmix(mix = post, quantiles)

    return(list(est = est, sd = sd, conf = conf))
  }

  # --- Calculate stats for post1 ---
  result1 <- calculate_post_stats(post1)

  SINGLE_ARM <- is.null(post2)

  if (!SINGLE_ARM){
    result2 <- calculate_post_stats(post2)

    # Calculate comparison statistics (Difference or Ratio)
    est_compare_original_scale <- result1$est - result2$est
    sd_compare_original_scale <- sqrt((result1$sd)^2 + (result2$sd)^2)
    conf_compare_original_scale <- tryCatch({
      qmixdiff(post1, post2, quantiles)
    }, error = function(e) {
      warning("Error calculating credible interval for the difference using qmixdiff: ", e$message)
      rep(NA_real_, 2)
    })

    # Apply exponent transformation if requested
    est1_final <- if (EXP_TRANSFORM) exp(result1$est) else result1$est
    ci1_l_final <- if (EXP_TRANSFORM) exp(result1$conf[1]) else result1$conf[1]
    ci1_u_final <- if (EXP_TRANSFORM) exp(result1$conf[2]) else result1$conf[2]

    est2_final <- if (EXP_TRANSFORM) exp(result2$est) else result2$est
    ci2_l_final <- if (EXP_TRANSFORM) exp(result2$conf[1]) else result2$conf[1]
    ci2_u_final <- if (EXP_TRANSFORM) exp(result2$conf[2]) else result2$conf[2]

    # est_compare is ratio if EXP_TRANSFORM, difference otherwise
    est_compare_final <- if (EXP_TRANSFORM) exp(est_compare_original_scale) else est_compare_original_scale
    compare_ci_l_final <- if (EXP_TRANSFORM) exp(conf_compare_original_scale[1]) else conf_compare_original_scale[1]
    compare_ci_u_final <- if (EXP_TRANSFORM) exp(conf_compare_original_scale[2]) else conf_compare_original_scale[2]

    # --- Construct Results Data Frame (Two Arms) ---
    rslt <- data.frame(
      est1 = unname(est1_final),
      sd1 = unname(result1$sd), # SD is always on original scale
      ci1_l = unname(ci1_l_final), ci1_u = unname(ci1_u_final),
      est2 = unname(est2_final),
      sd2 = unname(result2$sd), # SD is always on original scale
      ci2_l = unname(ci2_l_final), ci2_u = unname(ci2_u_final),
      est_compare = unname(est_compare_final),
      sd_compare = unname(sd_compare_original_scale), # SD is always of the difference on original scale
      compare_ci_l = unname(compare_ci_l_final), compare_ci_u = unname(compare_ci_u_final)
    )
  } else { # SINGLE_ARM
    est1_final <- if (EXP_TRANSFORM) exp(result1$est) else result1$est
    ci1_l_final <- if (EXP_TRANSFORM) exp(result1$conf[1]) else result1$conf[1]
    ci1_u_final <- if (EXP_TRANSFORM) exp(result1$conf[2]) else result1$conf[2]

    rslt <- data.frame(
      est1 = unname(est1_final),
      sd1 = unname(result1$sd), # SD is always on original scale
      ci1_l = unname(ci1_l_final), ci1_u = unname(ci1_u_final),
      est2 = NA_real_, # Use NA_real_ for numeric NA
      sd2 = NA_real_,
      ci2_l = NA_real_,
      ci2_u = NA_real_,
      est_compare = NA_real_,
      sd_compare = NA_real_,
      compare_ci_l = NA_real_,
      compare_ci_u = NA_real_
    )
  }

  return(rslt)
}
