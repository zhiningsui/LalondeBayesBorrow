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
Here's a simple example of how to use the package to analyze a single dataset with a continuous endpoint.
```R
library(LalondeBayesBorrow)

# 1. Prepare your data
# For a single analysis, you just need a data frame with summary statistics.
# Let's assume you have a treatment and a control arm.
current_data <- data.frame(
  treatment.n = 50,
  treatment.mu_hat = 1.8,
  treatment.s = 0.5,
  control.n = 50,
  control.mu_hat = 1.5,
  control.s = 0.4
)

# You can also include historical data to borrow from.
historical_data <- data.frame(
  control_h.n = 100,
  control_h.mu_hat = 1.6,
  control_h.s = 0.4
)

# 2. Set up your prior parameters
# These parameters define the non-informative and informative priors.
# We'll use a SAM prior, so we set `control.w = NULL`.
prior_params <- list(
  treatment.theta0 = 0,   # Non-informative prior mean for treatment
  treatment.s0 = 100,     # Non-informative prior SD for treatment
  control.theta0 = 1.6,   # Non-informative prior mean for control
  control.s0 = 100,       # Non-informative prior SD for control
  control.w = NULL,       # Use SAM prior for weight
  control.delta_SAM = 0.2,  # Clinically significant difference for SAM
  control.ess_h = 50      # Effective sample size from historical data
)

# 3. Run the analysis
# The main function is `lalonde_analysis`.
results <- lalonde_analysis(
  endpoint = "continuous",
  current_data = current_data,
  historical_data = historical_data,
  prior_params = prior_params,
  lrv = 0,         # Lower reference value
  tv = 0.2,        # Target value
  fgr = 0.2,       # False Go Rate
  fsr = 0.1        # False Stop Rate
)

# 4. View the results
# The results list contains the posterior inference and the go/no-go decision.
print(results$post_inference)
print(results$decision)
```
