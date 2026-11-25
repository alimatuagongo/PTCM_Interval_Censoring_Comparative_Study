# A Comparative Study of Promotion Time Cure Models for Interval-Censored Survival Data

This repository contains the R code for a comparative study of statistical and machine learning approaches for incidence modeling in the Promotion Time Cure Model (PTCM) framework for interval-censored survival data.

The estimation framework incorporates an Expectation-Maximization (EM) algorithm with a Weibull proportional hazards model for the latency component.

## Methods

Seven incidence modeling approaches are implemented and compared:

- Logistic Regression (LOGIT)
- Natural Spline Regression (SPLINE)
- Support Vector Machine (SVM)
- Decision Tree (DT)
- Random Forest (RF)
- Extreme Gradient Boosting (XGBoost)
- Neural Network (NN)

## Repository Contents

### ML_IC_PCM_Simulation.R

This script implements the simulation study, including:

- Generation of interval-censored survival data in the presence of a cured population.
- Three incidence settings with varying levels of complexity.
- EM-based estimation of incidence and latency components.
- Hyperparameter tuning for machine learning methods.
- Training and testing data evaluation.
- Estimation of uncured probabilities, population survival functions, susceptible survival functions, and promotion time survival functions.
- Monte Carlo performance evaluation.

The performance measures include:

- Bias
- Absolute bias
- Mean squared error (MSE)
- Area under the receiver operating characteristic curve (AUC)

The simulation is configured for 300 Monte Carlo replications. Sample sizes, simulation settings, and incidence modeling methods can be modified within the script.

### ML_IC_PCM_NDATA.R

This script implements the real-data application using the `hdsd` dataset available in the `thregI` R package.

The analysis includes:

- Data preparation and preprocessing.
- Implementation of all seven incidence modeling approaches.
- Weibull proportional hazards latency modeling.
- Stratified 10-fold cross-validation.
- Predictive performance comparison using AUC.
- Latency parameter estimation.
- Optional bootstrap-based statistical inference.
- Export of model comparison results.

## Required R Packages

The following R packages are used in the analysis:

- survival
- splines
- e1071
- rpart
- randomForest
- xgboost
- neuralnet
- MASS
- caTools
- dplyr
- PRROC
- pROC
- xtable
- readxl
- writexl
- thregI

## Running the Code

1. Download or clone this repository.
2. Install the required R packages.
3. Open the relevant R script in R or RStudio.
4. Specify the desired simulation settings or modeling approaches.
5. Run the script to obtain the corresponding results.

The simulation study may require substantial computational time, depending on the selected methods, sample sizes, and number of Monte Carlo replications.

## Author

**Alimatu Agongo**

Department of Mathematics  
University of Texas at Arlington
