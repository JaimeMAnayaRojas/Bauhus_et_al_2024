
# Statistical Analysis of Sleep-Like Behavior in Fish using Hidden Markov Models

This repository provides the statistical methods used in the paper "Effects of Parasites on Sleep-Like Behavior in Fish," where we model sleep-like behavior in fish using locomotor activity data and Hidden Markov Models (HMMs). The analysis was performed using the `hmmTMB` R package.

## Overview

In this study, we used locomotor activity per minute as time series data for each fish to fit a Hidden Markov Model (HMM), allowing us to identify periods of sleep-like behavior. HMMs are particularly useful for organisms unsuitable for electroencephalographic analyses, as they help model hidden patterns such as sleep based on observable behaviors like locomotion.

### Key Features

- **Handling Missing Data**: We assigned missing observations of locomotor activity to `NA` to preserve the time-series structure. A small fraction of the data (4.33%) was missing, and zero values (0.0057%) were adjusted slightly above zero to avoid additional parameters in the model.
  
- **State-Dependent Modeling**: We modeled the state-dependent distributions of locomotor activity as gamma distributions, with the assumption that all individuals followed the same state-dependent processes. We used a three-state HMM to represent:
  - **State 1**: Sleep-like behavior (lowest locomotion per minute)
  - **State 2**: Moderate activity (intermediate locomotion)
  - **State 3**: High activity (highest locomotion)

- **Diel Activity Patterns**: State-switching probabilities were modeled as a function of time of day, using trigonometric functions with a 24-hour wavelength as covariates. Random intercepts per fish were included to account for individual heterogeneity.

- **Model Fit and Validation**: The model fit was assessed by simulating data from the fitted HMM and calculating pseudo-residuals. Though some lack of fit was observed in the tails of the distribution, the overall fit was satisfactory.

### Tools Used
- [hmmTMB](https://cran.r-project.org/web/packages/hmmTMB/index.html): R package for Hidden Markov Models using TMB (Template Model Builder).

## Repository Structure

```
📁 data/          # Contains locomotor activity data for fish
📁 scripts/       # R scripts used for model fitting and analysis
📁 results/       # Outputs and figures generated from the analysis
📄 README.md      # Overview of the project (this file)
```

## Installation

To replicate the analysis, ensure you have R installed along with the `hmmTMB` package.

```r
install.packages("hmmTMB")
```

## Usage

1. **Prepare Data**: The data should be formatted as time series of locomotor activity per minute, with missing values denoted as `NA`.

2. **Fit the HMM**: Use the provided R scripts to fit a Hidden Markov Model to the data. The main script for fitting the model is `R/1_HMM_Modeling.R`.

3. **Analyze Results**: After fitting the model, use the scripts in the `results/` folder to decode the state sequences and generate figures for diel patterns, state distributions, and model diagnostics.

## Citation

If you use these methods in your research, please cite our paper:

> [Insert Paper Citation Here]

## Contributing

If you find any issues or have suggestions for improving the methods, feel free to open an issue or submit a pull request.

