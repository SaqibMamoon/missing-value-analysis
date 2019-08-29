# CONFIGURATION FILE

################################################################################
################################# Data Paths ###################################
################################################################################

path <- "/home/tabea/Documents/Code/outcome-prediction/missing-data/"
data_path <- paste(path, "data/MissionTplus_missingdataanalysis.csv", sep="")
path_to_plots <- paste("/home/tabea/Documents/Code/outcome-prediction/",
                       "missing-data/figures/", sep="")


################################################################################
####################### General parameter configurations #######################
################################################################################

# numeric "num" or categorical "cat" parameters
param_mode <- "cat"
# missing value mechanisms: missing at random (MAR) and 
#                           missing completely at random (MCAR)
mv_modes <- c("MAR", "MCAR")
# decide on percentages of missing values
percent_mv_vec <- seq(1, 60, 1)
# imputation methods
imp_methods <- c("Mean", "Median", "MICE", "Hot-deck", "EM")
# imputation methods including listwise deletion
imp_methods_lw <- c("Mean", "Median", "MICE", "Hot-deck", "EM", 
                    "Listwise deletion")

# number of shuffles for analysis
nr_iterations <- 100
# number iterations for mice and expectation maximization
nr_iter_mice <- 10
nr_iter_em <- 10

# plot per clinical parameter
# if FALSE, the error is averaged for all num/cat parameters
per_param <- FALSE

# all numerical parameters
column_name_num <- c("HOURS.TO.MRT1", "age", 
                     "mRS_before_stroke_discuss_inclusion", "NIHSS_scan1",
                     "NIHSS_scan3", "mRS_discharge", 
                     "TOAST_discharge_discuss_inclusion")
# all categorical parameters
column_name_cat <- c("gender", "treatment", "scan1_occlusion",
                     "risk_factor_hyperlipidemia", "risk_factor_Diabetes",
                     "risk_factor_hypertonia")


################################################################################
##################### Parameters for the Predictive Model ######################
################################################################################

# all parameters that are used for the prediction
column_names_pred <- c("HOURS.TO.MRT1", "age",
                       "mRS_before_stroke_discuss_inclusion", "NIHSS_scan1",
                       "gender", "treatment", "scan1_occlusion", 
                       "risk_factor_hyperlipidemia", "risk_factor_Diabetes",
                       "risk_factor_hypertonia")

# all parameters that are used for the "per parameter predictive model" analysis
column_names_imp <- c("HOURS.TO.MRT1", "mRS_before_stroke_discuss_inclusion",
                      "NIHSS_scan1")
# train/test ratio
percent_train <- 0.8
# balance training set
balance_train <- TRUE
shuffle_train <- FALSE
# max percentage of missing values for listwise deletion
lwise_threshold <- 10
# max percentage of missing values for MAR alone
MAR_max <- 9
# max plotting percentage
plot_max_percentage <- 60


################################################################################
########################### Parameters for the Plots ###########################
################################################################################

plot_every_5 <- TRUE
plot_up_to_25 <- FALSE
plot_std <- TRUE
plot_line <- TRUE
without_listwise <- FALSE