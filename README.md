# LalondeBayesBorrow: Bayesian Go/No-Go Decision Making for Clinical Trials
`LalondeBayesBorrow` is an R package that provides a comprehensive framework for designing and analyzing clinical trials using the Bayesian Lalonde decision-making framework. It is particularly well-suited for early-stage trials where leveraging historical control data can improve efficiency and decision-making. The package supports both binary and continuous endpoint types.

`LalondeBayesBorrow` can be used for both simulation-based trial design and for the analysis of real-world clinical trial data.

## Key Features
- Bayesian Analysis: Perform Bayesian analysis for one or two-arm trials.
- Historical Data Borrowing: Incorporate historical data using the Self-Adaptive Mixture (SAM) prior and power priors, with a gated control.
- Go/No-Go Decisions: Apply the Lalonde decision framework to make clear "Go", "No-Go", or "Consider" decisions based on posterior probabilities and credible intervals.
- Flexible Data Generation: Simulate trial data for different endpoint types to evaluate operating characteristics (OCs) under different scenarios.
- Real-World Data Analysis: Easily analyze a single dataset from a completed or ongoing trial.
- Operating Characteristics: Calculate and visualize operating characteristics like the False Go Rate (FGR) and False Stop Rate (FSR).

## Installation  
You can install the development version of `LalondeBayesBorrow` from GitHub with:
```R
# install.packages("devtools")
devtools::install_github("zhiningsui/LalondeBayesBorrow")
```

## Getting Started: A Quick Example

Here's a simple example of how to use the package to run a simulation and visualize the results.

```R
library(LalondeBayesBorrow)
library(dplyr)
library(ggplot2)

# 1. Set up your simulation parameters
param_grid <- expand.grid(
  trt_n = c(45, 60),
  ctrl_p = seq(0.15, 0.45, 0.15),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n, trt_p = 0.2 + ctrl_p)

data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "binary")

# Define historical data and other parameters
data_gen_params_h <- list(control_h = list(name = "control_h", n = 180, p = 0.3))
prior_params <- list(control.delta_SAM = 0.1, control.ess_h = 45)
nsim <- 1000 # Use a larger number for more stable results

# 2. Run the simulation
bayes_results <- LalondeBayesBorrow::run_simulation(
  nsim = nsim,
  data_gen_params_list = data_gen_params_list,
  data_gen_params_h = data_gen_params_h,
  prior_params = prior_params,
  endpoint = "binary",
  lrv = 0.1,
  tv = 0.2
)

# 3. Visualize the results
# Plot the operating characteristics (Go/Consider/No-Go probabilities)
plot_oc(bayes_results$oc_all, x_var = "control.p", facet_vars = c("trt_n", "true_value.compare_true"))

# Plot the risk profile (Type I and Type II errors)
plot_risk(bayes_results$oc_all, x_var = "control.p.diff", facet_vars = "trt_n")
```

