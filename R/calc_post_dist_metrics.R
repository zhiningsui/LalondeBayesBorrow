#' Calculate Posterior Distribution Metrics
#'
#' @description
#' This function calculates various evaluation metrics for the posterior distribution, including the average bias, standard error, empirical standard deviation, and nominal coverage probability.
#' For the "g-score" endpoint, it calculates the log-median ratio, and for the "OR" endpoint, it calculates the rate difference.
#'
#' @param endpoint A string. Specifies the type of endpoint ("g-score" or "OR").
#' @param true_value A scalar. The true median ratio for the "g-score" endpoint or the true rate difference for the "OR" endpoint.
#' @param post_est_ci A data frame. Contains posterior inference obtained from `bayesian_lalonde_decision()`, including posterior estimates, standard errors, and 95\% credible intervals for each simulation repetition.
#'
#' @return A data frame containing the evaluation metrics for the posterior distribution, including:
#'   - `bias_avg`: The average bias of the estimate.
#'   - `sd_avg`: The average standard error of the estimate.
#'   - `sd_empirical`: The empirical standard deviation of the estimate.
#'   - `cp`: The nominal coverage probability (95\%).
#'   - Additional metrics specific to the "g-score" endpoint. More details can be found in `eval_gscore_approx_dist()`.
#' @export
#'
#' @examples
#' # Example with a continuous endpoint:
#' post1 <- convert_RBesT_mix(post = data.frame(w = 0.5, mu1 = 0, sigma1 = 1, mu2 = 0.5, sigma2 = 1),
#'                            endpoint = "continuous")
#' post2 <- convert_RBesT_mix(post = data.frame(w = 0.9, mu1 = 0, sigma1 = 1, mu2 = 0.5, sigma2 = 1),
#'                            endpoint = "continuous")
#' post_inference_result <- posterior_inference(post1, post2, quantiles = c(0.025, 0.975),
#'                          EXP_TRANSFORM = TRUE)
#'
#' calc_post_dist_metrics("g-score", true_value = 1.2, post_est_ci = post_inference_result)
#'
#' @import stats
calc_post_dist_metrics <- function(endpoint, true_value, post_est_ci) {
  if (endpoint == "g-score") {
    # Parameter of interest: log_delta = log(theta_t/theta_c)
    log_delta <- log(true_value)

    # Point estimate for each run: log_delta_hat = log(theta_t_hat) - log(theta_c_hat)
    log_delta_hat <- log(post_est_ci$est1) - log(post_est_ci$est2)

    # Average Bias of log_delta_hat
    bias_avg <- mean(log_delta_hat - log_delta)

    # Average Bias of estimated median ratio
    bias_avg_median_ratio1 <- exp(mean(log_delta_hat)) - true_value
    median_ratio_hat <- post_est_ci$est_compare
    bias_avg_median_ratio2 <- mean(median_ratio_hat - true_value)

    # Average SD (i.e., SE) of log_delta_hat
    sds <- post_est_ci$sd_compare
    sd_avg <- mean(sds)

    # Empirical SD of log_delta_hat
    sd_empirical <- sd(log_delta_hat)

    # Coverage probability
    ci <- data.frame(cil = log(post_est_ci$compare_ci_l),
                     ciu = log(post_est_ci$compare_ci_u))
    ci$coverage <- ifelse(ci$cil < log_delta & ci$ciu > log_delta, 1, 0)
    cp <- mean(ci$coverage)

    metrics <- data.frame(log_delta = log_delta,
                          bias_avg = bias_avg,
                          bias_avg_median_ratio1 = bias_avg_median_ratio1,
                          bias_avg_median_ratio2 = bias_avg_median_ratio2,
                          sd_avg = sd_avg,
                          sd_empirical = sd_empirical,
                          cp = cp)

  } else if (endpoint == "OR") {
    # Parameter of interest: delta = rate_t - rate_c
    delta <- true_value

    # Point estimate for each run: delta_hat = rate_t_hat - rate_c_hat
    delta_hat <- post_est_ci$est_compare

    # Average Bias of estimated rate difference
    bias_avg_rate_diff <- mean(delta_hat - delta)

    # Average SD (i.e., SE) of delta_hat
    sds <- post_est_ci$sd_compare
    sd_avg <- mean(sds)

    # Empirical SD of delta_hat
    sd_empirical <- sd(delta_hat)

    # Coverage probability
    ci <- data.frame(cil = post_est_ci$compare_ci_l,
                     ciu = post_est_ci$compare_ci_u)
    ci$coverage <- ifelse(ci$cil < delta & ci$ciu > delta, 1, 0)
    cp <- mean(ci$coverage)

    metrics <- data.frame(delta = delta,
                          bias_avg_rate_diff = bias_avg_rate_diff,
                          sd_avg = sd_avg,
                          sd_empirical = sd_empirical,
                          cp = cp)
  }

  return(metrics)
}


