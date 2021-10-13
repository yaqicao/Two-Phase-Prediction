# Two-Phase-Prediction
R code of the paper "Two-phase stratified sampling and analysis for predicting binary outcomes"

## Folder: R Package-TwoPhaseAccuracy
We have developed an R package "TwoPhaseAccuracy". In this package, function "evalTwoPhase" is for calculating AUC, TPR, and FPR, with two-phase estimators implemented via function tps() in  R package "osDesign" (Haneuse and others, 2011). Function "seTwoPhase" is for estimating the corresponding standard errors using bootstrap resampling. The input arguments of the two functions include threshold values for estimating TPR and FPR and a dataframe that contains case-control status, stratum membership, Phase I predictors, Phase II predictors, variable names, and a variable indicating whether a subject was selected into Phase II. Function "seTwoPhase" includes an additional input argument for the number of bootstrap samples. A third function, "summaryTwoPhase", summarizes the results from "evalTwoPhase" and "seTwoPhase" by outputing the estimates of AUC, TPR, and FPR together with their standard error estimates.

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

### TwoPhaseAccuracy-example.R
The file example.R includes one example of how to estimate TPR, FPR & AUC as well as calculating standard errors with bootstrapping under Two-phase study designs.

## Folder: Simulation-TwoPhasePrediction

### FUN_DATA_sampling.R
* Functions:
  * Dat_gen(): Generate full data
  * Dat_gencc(): Generate Phase II data with Case-control sampling design
  * Data_genBalanced(): Generate Phase II data with balanced design
  * Dat_genEb(): Generate Phase II data with R-balanced design(Proposed)
  * Dat_formatG()/Dat_format(): Format data with/without stratum G
 
### FUN_Estimation.R and FUN_Evaluation.R
* Functions:
  * Dat_fullFit(), Dat_aucfull(): Develop the logistic regression model and calculate the predictive performance measures using the full data as benchmark
  * Dat_mleFitG(), Dat_mleFit(): For the data with/without G, obtain the two-phase ML estimation of odds ratio parameters and the proposed estimator of F (C.D.F. of (X,Z))
  * Dat_auc_mle(): Estimation and asymptotic standard error of AUC.

## Main-Design.R

## Main-Post-stratification.R

* Case-control design + R-balanced (Balanced) post-stratification analysis method
* Balanced design + R-balanced post-stratification analysis method
