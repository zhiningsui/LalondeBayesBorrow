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
# First, ensure all necessary libraries are loaded
library(LalondeBayesBorrow)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 1. Set up your simulation parameters
# Define the grid of parameters for the current trial arms
param_grid <- expand.grid(
  trt_n = 60, # Using a single sample size for a simpler example
  ctrl_p = seq(0.15, 0.45, 0.15),
  stringsAsFactors = FALSE
) %>%
  mutate(ctrl_n = trt_n, trt_p = 0.2 + ctrl_p)

# Create a list of data generation parameter sets from the grid
data_gen_params_list <- lapply(apply(param_grid, 1, as.list),
                               create_data_gen_params, endpoint = "binary")

# 2. Define Historical Data and Borrowing Strategies
# Use create_data_gen_params for historical data as well
data_gen_params_h <- create_data_gen_params(
  list(ctrl_h_n = 180, ctrl_h_p = 0.3),
  endpoint = "binary"
)

# Define prior parameters for each strategy (Borrowing vs. No Borrowing)
prior_params_list <- list(
  Borrowing = list(
    control.w = NULL, # Use NULL for SAM Prior
    control.delta_SAM = 0.1,
    control.delta_gate = 0.1,
    control.ess_h = 45
  ),
  No_Borrowing = list(control.w = 0)
)

# 3. Run the Simulation
# The `run_simulation` function from the README does not exist.
# The correct method is to loop through scenarios and call `bayesian_lalonde_decision`.

nsim <- 500 # For stable results, a larger nsim (e.g., 1000+) is recommended.
lrv <- 0.1
tv <- 0.2
bayes_results_raw <- list() # This list will store the results

# Loop over each data generation scenario
for (i in seq_along(data_gen_params_list)) {
  data_gen_params <- data_gen_params_list[[i]]

  # Generate data for the current scenario
  sim_data <- data_gen_binary(
    arm_names = names(data_gen_params),
    nsim = nsim,
    n_list = lapply(data_gen_params, `[[`, "n"),
    prob_list = lapply(data_gen_params, `[[`, "p"),
    seed = 123 # for reproducibility
  )

  # Prepare the summary data frame for the Bayesian function
  summary_h_df <- data.frame(
    control_h.n = data_gen_params_h$control_h$n,
    control_h.count = data_gen_params_h$control_h$n * data_gen_params_h$control_h$p
  )
  data_summary_df <- cbind(sim_data$summary, summary_h_df)

  # Loop over each prior parameter set (e.g., Borrowing vs. No Borrowing)
  for (k in seq_along(prior_params_list)) {
    prior_params <- prior_params_list[[k]]

    # Combine all settings for tracking purposes
    settings_df <- as.data.frame(
      c(list(true_value = sim_data$true_value), data_gen_params, data_gen_params_h, prior_params)
    )

    # Run the core Bayesian decision function
    post <- bayesian_lalonde_decision(
      endpoint = "binary",
      data_summary = data_summary_df,
      settings = settings_df,
      prior_params = prior_params,
      arm_names = c("treatment", "control", "control_h"),
      lrv = lrv,
      tv = tv
    )

    # Append the results of this run to the master list
    bayes_results_raw[[length(bayes_results_raw) + 1]] <- post
  }
}

# 4. Post-Process the Results
# This step consolidates the raw output and calculates operating characteristics.
processed_results <- process_sim_results(bayes_results_raw)

# 5. Visualize the Results

# a) Plot of Operating Characteristics (Go/No-Go/Consider probabilities)
oc_plot <- plot_oc(
  oc_data = processed_results$oc_all,
  x_var = "borrow",
  facet_formula = . ~ true_value.compare_true,
  labs(
    title = "Operating Characteristics by Borrowing Strategy",
    x = "Borrowing Strategy",
    y = "Probability of Decision",
    subtitle = "Faceted by True Difference in Objective Response Rate"
  )
)
print(oc_plot)


# b) Plot of FGR and FSR
# The function `plot_risk` does not exist.
# First, create the risk data frame, then plot it using ggplot2.
risk_df <- create_risk_df(
  oc_all = processed_results$oc_all,
  lrv = lrv,
  tv = tv
)

# Reshape data for plotting
risk_plot_df <- risk_df %>%
  pivot_longer(cols = c(FGR, FSR), names_to = "risk_type", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(borrow = factor(borrow, levels = c("No", "Yes")))

# Create the plot
risk_plot <- ggplot(risk_plot_df, aes(x = borrow, y = value, color = risk_type, group = risk_type)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~true_value.compare_true) + # Facet by the true effect
  geom_hline(aes(yintercept = fgr, linetype = "FGR Threshold"), color = "red") +
  geom_hline(aes(yintercept = fsr, linetype = "FSR Threshold"), color = "blue") +
  scale_y_continuous(labels = percent) +
  labs(
    title = "Type I (FGR) and Type II (FSR) Error Rates",
    subtitle = "Faceted by True Difference in Objective Response Rate",
    x = "Borrowing Strategy",
    y = "Rate",
    color = "Risk Type",
    linetype = "Thresholds"
  ) +
  theme_bw()

print(risk_plot)


# 6. Create a Summary Table of the Posterior Mean Difference (PMD)
pmd_summary_df <- create_pmd_summary(
  post_inference_all = processed_results$post_inference_all,
  borrow_to_compare = c("Yes", "No")
)

print(pmd_summary_df)
```

