# A Public Framework for Testing Different Imputation Methods for Clinical Datasets


## Contents

1. [Project aim](#project-aim)
2. [Overview of folder structure](#folder-structure)
3. [Code files](#code-files)
4. [How to run the analysis](#how-to-run-the-analysis)


## Project aim

This repository provides a framework to evaluate imputation methods for 
different percentages of missing values. The analysis includes two aspects: 1) 
error analysis in terms of root mean squared errror (RMSE) and absolute bias 
and 2) predictive model analysis in terms of area under the curve (AUC) of a 
generalized linear model (GLM). Missing values are removed in two different 
fashions: 1) missing-at-random (MAR) and 2) missing-completely-at-random (MCAR). 
The following imputation methods are compared:

1. mean imputation
2. median imputation
3. multiple imputation by chained equations (MICE)
4. hot-deck ("sampling")
5. expectation maximization (EM)
6. listwise deletion (for predictive model analysis only)


## Folder structure

```
.
+-- code
+-- data
+-- imputations
|   +-- different-samples
|   +-- same-samples
+-- results
|   +-- error-analysis
|   |   +-- different-samples
|   |   +-- same-samples
|   +-- predictive-model
|   |   +-- across-parameters
|   |   +-- per-parameter
+-- figures
|   +-- error-analysis
|   |   +-- per-parameter
|   +-- predictive-model
|   |   +-- with-listwise
|   |   +-- without-listwise
```

All executable code can found in the code folder and is further described in the
[code files section](#code-files). The data folder contains the data that should
be analyzed. 

In the imputations folder all calculated imputations for the 
different imputation methods are stored. The imputations are computed in two 
different fashions: 1) different samples and 2) same samples. Using different 
samples, different missing values are imputed in each iteration for every 
imputation method. It is used to calculate the root mean squared error (RMSE). 
When applying same samples, the same missing values are imputed in each 
iteration for every imputation method. This is crucial in order to determine the 
absolute bias of the imputation methods.

The results folder contains both the results of the error analysis (RMSE and 
absolute bias) and the results of the predictive model analysis (area under the 
curve, AUC). The predictive model folder is split into across parameters and per 
parameter. Across parameters means that values are removed across parameters for
each percentage. In the per parameter fashion values are removed for one 
parameter at a time for each percentage.

In the figures folders all plots from the different analyses are saved.


## Code files

* __General code files__:
    * *config.R*: contains all configurations of both analyses.
    * *helper_functions.R*: contains all helper functions for both analyses.
* __Error analysis__:
    * *run_imputation_per_parameter.R*: calculates imputations for every method 
                                        and parameter and stores it.
    * *evaluate_error.R*: computes the RMSE and absolute bias and stores it.
    * *plot_error_results.R*: plots the error analysis averaged over numerical
                              and categorical data and stores it.
    * *plot_error_results_per_param.R*: plots the error analysis per parameter
                                        and stores it.
    * *run_error_statistical_analysis.R*: calculates the p values comparing the 
                                          RMSE of the different imputation 
                                          methods.
* __Predictive model analysis__:
    * *run_pred_model_across_parameters.R*: computes the AUCs for a GLM model 
                                            using imputed data. Missing values 
                                            were introduced across parameters.
                                            Resulting AUCs are then stored.
    * *run_pred_model_per_parameter.R*: calculates the AUCs for a GLM model 
                                        using imputed data. Missing values were
                                        introduced for one parameter at a time.
                                        Resulting AUCs are then stored.
    * *plot_perf_across_parameters.R*: plots the AUCs and corresponding p 
                                       values and saves the plot.
    * *plot_perf_per_parameter.R*: plots the AUCs resulted from 
                                   *run_pred_model_per_parameter.R*.
    * *run_pred_model_statistical_analysis.R*: runs a more detailed statistical
                                               analysis comparing the different
                                               imputation methods between each 
                                               other and to the full model 
                                               (model using full dataset).
                                               

## How to run the analysis

__General__:

1. Put your data into the data folder and change the file path in *config.R*.
2. Change parameters (if necessary) in the *config.R* file.
3. Adapt the *preprocess* and *definedatatypes* functions in the 
   *helper_functions.R* file according to your dataset. Decide which of your 
   parameters are numerical and which are categorical and adapt this in the
   *config.R* file.

__Error analysis__:

4. Run *run_imputation_per_parameter.R* and *evaluate_error.R* for your 
   numerical and categorical data to evaluate the error.
5. Run *plot_error_results.R* and *run_error_statistical_analysis.R* for 
   plotting and the respective p values.

__Predictive model analysis__:

4. Decide on a parameter that should be predicted and change it in the 
   respective predictive model files.
5. Run *run_pred_model_across_parameter.R* and *run_pred_model_per_parameter.R*
   for the AUC evaluation.
6. Run *plot_perf_across_parameters.R* and *plot_perf_per_parameter.R* to get 
   the plots.
7. Run *run_pred_model_statistical_analysis.R* for the statistical analysis.
