# PREDICTIVE MODEL: script to run GLM per parameter

# import libraries
library(dplyr)
library(glmnet)
library(abind)
library(ROCR)
library(pROC)

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

# preprocessing of data
df <- preprocess(data_path)

# exclude all non-complete cases
df_cc <- na.omit(df)

# get labels
df_cc$mRS_discharge[df_cc$mRS_discharge<=2] <- 0 # good outcome
df_cc$mRS_discharge[df_cc$mRS_discharge>2] <- 1  # bad outcome
y <- df_cc$mRS_discharge

# get rid of discharge parameters
df_cc <- subset(df_cc, select=column_names_pred)

# normalize numerical data
df_norm <- normalizedata(df_cc,intersect(column_names_pred, column_name_num))
# reset row indices
rownames(df_norm) <- NULL

result_path <- paste(path, "results/predictive-model/", sep="")

for (mv_mode in mv_modes) {
  # for results storage
  blank_data <- rep(NaN, length(column_names_imp) * (length(imp_methods) + 1) *
                    nr_iterations * length(percent_mv_vec))
  AUC_all <- array(blank_data, c(length(percent_mv_vec),
                                 length(column_names_imp),
                                 length(imp_methods) + 1, nr_iterations))
  acc_all <- array(blank_data, c(length(percent_mv_vec), 
                                 length(column_names_imp), 
                                 length(imp_methods) + 1, nr_iterations))
  
  # loop over every percentage
  for (k in array(1:length(percent_mv_vec))) {
    
    blank_data <- rep(NaN, length(column_names_imp) * 
                      (length(imp_methods) + 1) * nr_iterations)
    AUC_per_percent <- array(blank_data, c(length(column_names_imp),
                                           length(imp_methods) + 1,
                                           nr_iterations))
    acc_per_percent <- array(blank_data, c(length(column_names_imp), 
                                           length(imp_methods) + 1,
                                           nr_iterations))
    
    percent_mv <- percent_mv_vec[k]
    
    print(paste("Percentage: ", percent_mv, "/", length(percent_mv_vec)))
    
    # load percentage file for numerical data
    filename_output <- get_output_path(percent_mv, "num", mv_mode,
                                       nr_iterations)
    load(paste(path, "imputations/different-samples/", filename_output, sep=""))
    # store everything in extra variables and get rid of discharge parameters
    idx_imp_num <- idx_imp_all[1:4,,]
    label_num <- label_all[1:4,,]
    imputed_values_num <- imputed_values[1:4,,,]
    
    # load percentage file for categorical data
    filename_output <- get_output_path(percent_mv, "cat", mv_mode, 
                                       nr_iterations)
    load(paste(path, "imputations/different-samples/", filename_output, sep=""))
    # store everything in extra variables
    idx_imp_cat <- idx_imp_all
    label_cat <- label_all
    imputed_values_cat <- imputed_values
    
    # concatenate num and cat data
    idx_imp_all <- abind(idx_imp_num, idx_imp_cat, along=1)
    label_all <- abind(label_num, label_cat, along=1)
    imputed_values <- abind(imputed_values_num, imputed_values_cat, along=1)
  
    # store original dataframe in df_real
    df_real <- df_norm
    # get df_imp as a dataframe for the imputations
    df_imp <- df_norm
    
    param_nr <- 0
    
    # parameter selection: based on GLM on all data, 
    # highest coefficients in the 1000Plus dataset: NIHSS_scan1, HOURS.TO.MRT1, 
    # mRS_before_stroke_discuss_inclusion
    for (p in 1:length(column_names_pred)) {
      
      # if param is part of the list
      if (column_names_pred[p] %in% column_names_imp) {
        param_nr <- param_nr + 1
        
        # for every imputation i
        for (i in 1:nr_iterations) {
          # for every imputation method m
          for (m in 1:(length(imp_methods) + 1)) {
            # replace
            if (m==(length(imp_methods)+1)) {
              replace_with_NA <- rep(NA, dim(imputed_values)[3])
              df_imp[[column_names_pred[p]]] <- 
                replace(df_imp[[column_names_pred[p]]], idx_imp_all[p,,i],
                        replace_with_NA)
            } else {
              df_imp[[column_names_pred[p]]] <- 
                replace(df_imp[[column_names_pred[p]]], idx_imp_all[p,,i],
                        imputed_values[p,m,,i])
            }    
            
            # make sure that both classes are in the training and testing set
            y_train <- 0
            y_test <- 0
            counter <- 0
            while ((length(unique(y_train))==1) || 
                   (length(unique(y_test))==1)) {
              # split randomly into train and test
              train_idx <- sample.int(dim(df_imp)[1], floor(dim(df_imp)[1] *
                                                      percent_train))
              # other 20% are testing data
              test_idx <- 1:dim(df_imp)[1]
              test_idx <- test_idx[-train_idx]
              # creating the subsets for data and labels
              df_train <- df_imp[train_idx,]  # imputations for train
              df_test <- df_real[test_idx,]   # real data set for test
              y_train <- y[train_idx]
              y_test <- y[test_idx]
              # for listwise del delete rows here
              if (m==(length(imp_methods)+1)) {
                df_tmp <- cbind(df_train, y_train)
                df_tmp <- na.omit(df_tmp)
                y_train <- df_tmp$y_train
                df_train <- subset(df_tmp, select=-c(y_train))
              }
              # reset indices
              rownames(df_train) <- NULL
              rownames(df_test) <- NULL
              counter <- counter + 1
              if (counter==10000) {
                print("Might be caught in infinite loop")
              }
            }
            
            # balance training labels 
            if (balance_train) {
              df_train_bad <- df_train[y_train==1,] # take all bad labels
              df_train_good <- df_train[y_train==0,]
              rownames(df_train_bad) <- NULL
              rownames(df_train_good) <- NULL
              good_sample_idx <- sample.int(dim(df_train_good)[1],
                                            dim(df_train_bad)[1])
              df_train_good <- df_train_good[good_sample_idx,]
              rownames(df_train_good) <- NULL
              df_train <- rbind(df_train_bad, df_train_good)
              y_train <- cbind(t(rep(1,dim(df_train_bad)[1])),
                               t(rep(0,dim(df_train_bad)[1])))
              y_train <- drop(y_train)
              rownames(df_train) <- NULL
            }
            
            # shuffle train set
            if (shuffle_train) {
              shuffle_idx <- sample.int(dim(df_train)[1], dim(df_train)[1])
              df_train <- df_train[shuffle_idx,]
              y_train <- y_train[shuffle_idx]
              rownames(df_train) <- NULL
            }
            
            # fit model on train and predict on test
            fit_glm <- glm(y_train~., family=binomial(link="logit"),
                           data=df_train)
            p_test <- predict(fit_glm, df_test, type="response")
            pr_test <- prediction(p_test, y_test)
            
            # evaluate on test
            auc_test <- performance(pr_test, measure="auc")
            AUC_per_percent[param_nr,m,i] <- auc_test@y.values[[1]]
            acc_test <- performance(pr_test, measure="acc")
            acc_per_percent[param_nr,m,i] <- max(acc_test@y.values[[1]])
          
            # "reset" data frames
            df_real <- df_norm
            df_imp <- df_norm
          }
        }
      }
    }
    # save AUC per percent
    save(AUC_per_percent, acc_per_percent,
         file=paste(result_path, "per-parameter/", "GLM_", mv_mode, "_",
                    percent_mv, ".RData", sep=""))
    
    AUC_all[k,,,] <- AUC_per_percent
    acc_all[k,,,] <- acc_per_percent
  }
  # save AUC for all
  save(AUC_all,acc_all,file=paste(result_path, "per-parameter/", "GLM_",
                                  mv_mode, "_all", ".RData", sep=""))
}