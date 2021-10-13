# Two-Phase-Prediction
R code of the paper "Two-phase stratified sampling and analysis for predicting binary outcomes"

## R Package-TwoPhaseAccuracy
### TwoPhaseAccuracy-functions.R
* Key Notations
  *  Outcome: outcome information (1: case; 0: control)
  *  X: standard predictors whose information is available to all subjects
  *  Z: new biomarkers whose information is only available to a subset of subjects
  *  Phase_ID: information of whether a subject is Phase I or Phase II (1: phase I; 2: phase II)
  *  Stratum: stratification variables defined based on X
* Methods:
  *  Benchmark: fit logistic regression model to all subjects: logitP(Y=1|X) ~ X
  *  MLE: fit logistic regression model to Two-phase dataset with maximum likelihood: logitP(Y=1|X,Z) ~ X+Z
* Predictive Accuracy Measures:
  *  TPR(q), FPR(p), AUC
Main functions of estimating TPR, FPR & AUC with bootstrapping method for standard error calculations under two-phase study design 
* Functions:
  * evalTwoPhase(): estimate predictive accuracy measures
  * seTwoPhase(): estimate bootstrap standard errors
  * evalTwoPhase.fun(): helper function used in seTwoPhase()
  * summaryTwoPhase(): print results
* Input Arguments:
  * Outcome: vector of length n 
  * X: matrix of size n x t (t standard predictors included in the model)
  * Z: matrix of size n x m (m biomarkers included in the model); some subjects have NA values
  * Stratum: vector of length n
  * Phase_ID: vector of length n
  * namesX: names of X; vector of length t
  * namesZ: names of Z; vector of length m
  * method: "Benchmark", "ML"
  * q: threshold value for TPR(q)
  * p: threshold value for FPR(p)
  * numBoot: number of bootstrap samples
