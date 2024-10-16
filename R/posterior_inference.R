#' Perform Posterior Inference
#'
#' This function calculates posterior estimates, standard deviations, and credible intervals for individual arms and the treatment effect, based on posterior distributions from the RBesT package.
#' If a second posterior distribution is not provided, it calculates these statistics for a single arm.
#' The function also supports exponentiating the results, which can be useful for endpoints like log-transformed data.
#'
#' @param post1 A RBesT mixture distribution object. The first posterior distribution (e.g., for the treatment arm).
#' @param post2 A RBesT mixture distribution object. The second posterior distribution (e.g., for the control arm). Default is `NULL`, in which case only posterior estimates for `post1` are calculated.
#' @param quantiles A vector. The upper and lower quantiles for the credible interval (e.g., `c(0.025, 0.975)` for a 95\% credible interval).
#' @param EXP_TRANSFORM A logical. If TRUE, the function exponentiates the results, which is useful for log-transformed data.
#'
#' @return A data frame with one row containing:
#'   - \code{est1}: Posterior mean for the first posterior distribution (e.g., treatment arm).
#'   - \code{sd1}: Posterior standard deviation for the first posterior distribution.
#'   - \code{ci1_l}: Lower bound of the credible interval for the first posterior distribution.
#'   - \code{ci1_u}: Upper bound of the credible interval for the first posterior distribution.
#'   - \code{est2}: Posterior mean for the second posterior distribution (if provided).
#'   - \code{sd2}: Posterior standard deviation for the second posterior distribution (if provided).
#'   - \code{ci2_l}: Lower bound of the credible interval for the second posterior distribution (if provided).
#'   - \code{ci2_u}: Upper bound of the credible interval for the second posterior distribution (if provided).
#'   - \code{est_compare}: Estimated treatment effect, either as a difference (for continuous data) or a ratio (for log-transformed data if `EXP_TRANSFORM = TRUE`).
#'   - \code{sd_compare}: Standard deviation of the estimated treatment effect.
#'   - \code{compare_ci_l}: Lower bound of the credible interval for the treatment effect.
#'   - \code{compare_ci_u}: Upper bound of the credible interval for the treatment effect.
#' @export
#'
#' @examples
#' # Example with one posterior distribution (single-arm inference):
#' post1 <- convert_RBesT_mix(post = data.frame(w = 0.5, mu1 = 0, sigma1 = 1,
#'                                              mu2 = 0.5, sigma2 = 1),
#'                            endpoint = "continuous")
#' post_inference_result <- posterior_inference(post1, quantiles = c(0.025, 0.975),
#'                                              EXP_TRANSFORM = FALSE)
#'
#' # Example with two posterior distributions (treatment vs. control):
#' post2 <- convert_RBesT_mix(post = data.frame(w = 0.9, mu1 = 0, sigma1 = 1,
#'                                              mu2 = 0.5, sigma2 = 1),
#'                            endpoint = "continuous")
#' post_inference_result <- posterior_inference(post1, post2,
#'                                              quantiles = c(0.025, 0.975),
#'                                              EXP_TRANSFORM = TRUE)
#'
#' @import RBesT
posterior_inference <- function(post1, post2 = NULL, quantiles, EXP_TRANSFORM = FALSE) {

  calculate_post <- function(post, EXP_TRANSFORM = FALSE) {
    est <- summary(post)["mean"]
    sd <- summary(post)["sd"]
    conf <- qmix(mix = post, quantiles)
    if (EXP_TRANSFORM) {
      est <- exp(est)
      conf <- exp(conf)
    }
    return(list(est = est, conf = conf, sd = sd))
  }

  SINGLE_ARM <- is.null(post2)

  result1 <- calculate_post(post1, EXP_TRANSFORM)

  if (!SINGLE_ARM){
    result2 <- calculate_post(post2, EXP_TRANSFORM)

    est_compare <- if (EXP_TRANSFORM) {
      result1$est / result2$est
    } else {
      result1$est - result2$est
    }

    conf_compare <- tryCatch({
      if (EXP_TRANSFORM) {
        exp(qmixdiff(post1, post2, quantiles))
      } else {
        qmixdiff(post1, post2, quantiles)
      }
    }, error = function(e) {
      NA
    })

    sd_compare <- sqrt((result1$sd)^2 +(result2$sd)^2)

  }

  rslt <- data.frame(
    est1 = result1$est,
    sd1 = result1$sd,
    ci1_l = result1$conf[1], ci1_u = result1$conf[2],
    est2 = if (!SINGLE_ARM) result2$est else NA,
    sd2 = if (!SINGLE_ARM) result2$sd else NA,
    ci2_l = if (!SINGLE_ARM) result2$conf[1] else NA,
    ci2_u = if (!SINGLE_ARM) result2$conf[2] else NA,
    est_compare = if (!SINGLE_ARM) est_compare else NA,
    sd_compare = if (!SINGLE_ARM) sd_compare else NA,
    compare_ci_l = if (!SINGLE_ARM) conf_compare[1] else NA,
    compare_ci_u = if (!SINGLE_ARM) conf_compare[2] else NA
  )

  return(rslt)
}
