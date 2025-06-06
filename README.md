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
```R
devtools::install_github("zhiningsui/LalondeBayesBorrow")
