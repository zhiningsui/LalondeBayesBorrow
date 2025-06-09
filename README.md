# LalondeBayesBorrow: Bayesian Go/No-Go Decision Making for Clinical Trials
`LalondeBayesBorrow` is an R package that provides a comprehensive framework for designing and analyzing clinical trials using the Bayesian Lalonde go/no-go decision-making framework. It is particularly well-suited for early-stage trials where leveraging historical control data can improve efficiency and decision-making.

The package allows you to simulate trial data, perform Bayesian analysis with historical borrowing, and evaluate the operating characteristics of your trial design. It is designed to be flexible, allowing for various endpoint types and borrowing strategies.

## Key Features
- Bayesian Analysis: Perform Bayesian analysis for one- or two-arm trials for continuous and binary endpoints.
- Historical Data Borrowing: Incorporate historical data using the Self-Adaptive Mixture (SAM) prior and power priors, with a gated control.
- Go/No-Go Decisions: Apply the Lalonde decision framework to make clear "Go", "No-Go", or "Consider" decisions based on posterior probabilities and credible intervals.
- Flexible Data Generation: Simulate trial data for different endpoint types to evaluate operating characteristics (OCs) under different scenarios.
- Real-World Data Analysis: Easily analyze a single dataset from a completed or ongoing trial.
- Publication-Ready Visualizations: Generate insightful plots of operating characteristics, risk profiles (FGR/FSR), and posterior distributions with a suite of user-friendly plotting functions.

## Installation  
You can install the development version of `LalondeBayesBorrow` from GitHub with:
```R
# install.packages("devtools")
devtools::install_github("zhiningsui/LalondeBayesBorrow")
```

## Usage: A Complete Workflow Example

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

# 2. Define Historical Data and Borrowing Strategies
data_gen_params_h <- list(control_h = list(name = "control_h", n = 180, p = 0.3))

prior_params_list <- list(
  Borrowing = list(control.delta_SAM = 0.1, control.ess_h = 45, control.w = NULL),
  No_Borrowing = list(control.w = 0)
)

# 3. Run the Simulation
# This is a small simulation for example purposes.
# For stable results, a larger nsim (e.g., 1000+) is recommended.
bayes_results_raw <- LalondeBayesBorrow::run_simulation(
  nsim = 100,
  data_gen_params_list = data_gen_params_list,
  data_gen_params_h = data_gen_params_h,
  prior_params_list = prior_params_list,
  endpoint = "binary",
  lrv = 0.1,
  tv = 0.2
)

# 4. Post-Process the Results
# This step consolidates the raw output and calculates operating characteristics.
processed_results <- process_sim_results(bayes_results_raw)

# 5. Visualize the Results
# Create a plot of the operating characteristics (Go/No-Go/Consider probabilities).
oc_plot <- plot_oc(
  oc_data = processed_results$oc_all,
  x_var = "borrow",
  facet_formula = . ~ true_value.compare_true,
  labs(
    title = "Operating Characteristics by Borrowing Strategy",
    x = "Borrowing Strategy",
    y = "Probability of Decision"
  )
)

print(oc_plot)

# Create a risk plot (False Go Rate and False Stop Rate)
risk_plot <- plot_risk(
  oc_all = processed_results$oc_all,
  lrv = 0.1,
  tv = 0.2,
  x_var = "borrow",
  labs(
      title = "Type I and Type II Error Rates",
      x = "Borrowing Strategy"
  )
)
print(risk_plot)

# 6. Create a Summary Table of the Posterior Mean Difference (PMD)
pmd_summary_df <- create_pmd_summary(
  post_inference_all = processed_results$post_inference_all,
  borrow_to_compare = c("Yes", "No")
)
print(pmd_summary_df)
```

