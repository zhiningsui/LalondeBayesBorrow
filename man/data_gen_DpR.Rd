% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_gen.R
\name{data_gen_DpR}
\alias{data_gen_DpR}
\title{Generate Constrained Depth of Response (DpR) data for multiple arms across multiple simulations}
\usage{
data_gen_DpR(n_arms, nsim, n_list, p_list, arm_names, mu_list, sigma_list)
}
\arguments{
\item{n_arms}{Number of arms (can be inferred from arm_names, but kept for consistency with template)}

\item{nsim}{Number of simulations}

\item{n_list}{Named list of sample sizes per arm (names should match arm_names)}

\item{p_list}{Named list of proportion of -1 per arm (names should match arm_names)}

\item{arm_names}{Character vector of arm names}

\item{mu_list}{Named list of means per arm (theta_x) (names should match arm_names)}

\item{sigma_list}{Named list of standard deviations per arm (sigma_x) (names should match arm_names)}
}
\value{
A list containing:
\itemize{
\item \code{data}: A data.table with columns \code{nsim}, \code{id}, \code{arm}, \code{value}, containing all simulated DpR values.
\item \code{summary}: A data.table with columns \code{nsim} and arm-specific summary statistics (n, mu_hat, sd, s).
\item \code{true_value}: A data frame with true mean differences (currently only 'treatment' vs 'control').
}
}
\description{
Generate Constrained Depth of Response (DpR) data for multiple arms across multiple simulations
}
\seealso{
\code{bayesian_lalonde_decision}
}
