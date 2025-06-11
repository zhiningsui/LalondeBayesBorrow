#' Calculate Posterior Distribution Metrics
#'
#' @description
#' This function calculates various evaluation metrics for the posterior distribution, including the average bias, standard error, empirical standard deviation, and nominal coverage probability.
#' For the "g-score" endpoint, it calculates the log-median ratio, and for the "OR" endpoint, it calculates the rate difference.
#'
#' @details
#' For the `"g-score"` endpoint, it calculates the log-median ratio, and for the `"OR"` endpoint, it calculates the rate difference.
#' - Average bias
#' - Standard error
#' - Empirical standard deviation
#' - Nominal coverage probability
#'
#' @param endpoint A string. Specifies the type of endpoint (`"g-score"` or `"OR"` or `"continuous"`).
#' @param true_value A scalar. The true median ratio for the `"g-score"` endpoint, the true rate difference for the `"OR"` endpoint, or the true mean difference for the `"continuous"` endpoint.
#' @param post_est_ci A data frame. Contains posterior inference obtained from `bayesian_lalonde_decision()`, including posterior estimates, standard errors, and 95% credible intervals for each simulation repetition.
#' @param remove.na Logical. Whether to remove NA values before calculation. Defaults to `FALSE`.
#'
#' @return A data frame containing the evaluation metrics for the posterior distribution, including:
#'   * `bias_avg`: The average bias of the estimate (specific column names like `bias_avg_median_ratio1`, `bias_avg_median_ratio2`, `bias_avg_rate_diff`, `bias_avg_mean_diff` depend on the endpoint).
#'   * `sd_avg`: The average standard error of the estimate (`sd_compare`).
#'   * `sd_empirical`: The empirical standard deviation of the estimate (`est_compare`).
#'   * `cp`: The nominal coverage probability (95%).
#'   * Additional metrics specific to the `"g-score"` endpoint (see `eval_gscore_approx_dist()`).
#' @export
#' @seealso `bayesian_lalonde_decision`, `eval_gscore_approx_dist`
#' @importFrom stats na.omit sd
calc_post_dist_metrics <- function(endpoint, true_value, post_est_ci, remove.na = FALSE) {
  if (endpoint == "g-score") {

    delta <- true_value
    # Parameter of interest: log_delta = log(theta_t/theta_c)
    log_delta <- log(delta)

    # Point estimate for each run: log_delta_hat = log(theta_t_hat) - log(theta_c_hat)
    log_delta_hat <- log(post_est_ci$est1) - log(post_est_ci$est2)

    # Average Bias of log_delta_hat
    bias_avg <- mean(log_delta_hat - log_delta)

    # Average Bias of estimated median ratio
    bias_avg_median_ratio1 <- exp(mean(log_delta_hat)) - delta
    median_ratio_hat <- post_est_ci$est_compare
    bias_avg_median_ratio2 <- mean(median_ratio_hat - delta)

    # Average SD (i.e., SE) of log_delta_hat
    sds <- post_est_ci$sd_compare
    sd_avg <- mean(sds)

    # Empirical SD of log_delta_hat
    sd_empirical <- stats::sd(log_delta_hat)

    # Coverage probability
    if(remove.na) {
      ci <- data.frame(cil = log(post_est_ci$compare_ci_l),
                       ciu = log(post_est_ci$compare_ci_u)) %>%
        na.omit()
    }else{
      ci <- data.frame(cil = log(post_est_ci$compare_ci_l),
                       ciu = log(post_est_ci$compare_ci_u))
    }

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
    sd_empirical <- stats::sd(delta_hat)

    # Coverage probability
    if(remove.na) {
      ci <- data.frame(cil = post_est_ci$compare_ci_l,
                       ciu = post_est_ci$compare_ci_u) %>%
        na.omit()
    }else{
      ci <- data.frame(cil = post_est_ci$compare_ci_l,
                       ciu = post_est_ci$compare_ci_u)
    }

    ci$coverage <- ifelse(ci$cil < delta & ci$ciu > delta, 1, 0)
    cp <- mean(ci$coverage)

    metrics <- data.frame(delta = delta,
                          bias_avg_rate_diff = bias_avg_rate_diff,
                          sd_avg = sd_avg,
                          sd_empirical = sd_empirical,
                          cp = cp)
  } else if (endpoint == "continuous") {

    delta <- true_value
    delta_hat <- post_est_ci$est_compare
    bias_avg_mean_diff <- mean(delta_hat - delta)

    sds <- post_est_ci$sd_compare
    sd_avg <- mean(sds)
    sd_empirical <- stats::sd(delta_hat)

    if(remove.na) {
      ci <- data.frame(cil = post_est_ci$compare_ci_l,
                       ciu = post_est_ci$compare_ci_u) %>% na.omit()
    } else {
      ci <- data.frame(cil = post_est_ci$compare_ci_l,
                       ciu = post_est_ci$compare_ci_u)
    }

    ci$coverage <- ifelse(ci$cil < delta & ci$ciu > delta, 1, 0)
    cp <- mean(ci$coverage)

    metrics <- data.frame(
      delta = delta,
      bias_avg_mean_diff = bias_avg_mean_diff,
      sd_avg = sd_avg,
      sd_empirical = sd_empirical,
      cp = cp
    )
  }

  return(metrics)
}


