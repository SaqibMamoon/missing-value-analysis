# PREDICTIVE MODEL: script to run GLM across all parameters

# import libraries
library(lattice)
library(mice)
library(Amelia)
library(dplyr)

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

# preprocessing of data, dataset-specific
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
df_norm <- normalizedata(df_cc, intersect(column_names_pred, column_name_num))

# reset row indices
rownames(df_norm) <- NULL
df_real <- df_norm

# storage for AUC
blank_data <- rep(NaN, length(percent_mv_vec) * length(imp_methods_lw) *
                  nr_iterations)
AUC_per_percent <- array(blank_data, c(length(percent_mv_vec),
                                       length(imp_methods_lw), nr_iterations))

# get performance for model without imputation
AUC_noimp <- rep(NaN, nr_iterations)
for (i in array(1:nr_iterations)) {
  train_idx <- sample.int(dim(df_norm)[1], floor(dim(df_norm)[1] * 
                                           percent_train))
  # other 20% are testing data
  test_idx <- 1:dim(df_norm)[1]
  test_idx <- test_idx[-train_idx]
  # creating the subsets for data and labels
  df_train <- df_norm[train_idx,]  
  df_test <- df_norm[test_idx,]  
  y_train <- y[train_idx]
  y_test <- y[test_idx]
  AUC_noimp[i] <- glm_model(df_train, df_test, y_train, y_test)
}

for (mv_mode in mv_modes) {

  print(mv_mode)
  
  start_time <- Sys.time()
  
  for (k in array(1:length(percent_mv_vec))) {
    
    percent_mv <- percent_mv_vec[k]
    
    print(paste("Percentage: ", percent_mv, "/", max(percent_mv_vec)))
    
    for (i in array(1:nr_iterations)) {
      
      #if ((i%%10)==0) {
      #  print(paste(i,"/",nr_iterations,sep=""))
      #}
      
      # training/test split
      train_idx <- sample.int(dim(df_norm)[1], floor(dim(df_norm)[1] *
                                               percent_train))
      test_idx <- 1:dim(df_norm)[1]
      test_idx <- test_idx[-train_idx] # test and train disjoint
      df_train <- df_norm[train_idx,]  
      df_test <- df_norm[test_idx,]    
      y_train <- y[train_idx]
      y_test <- y[test_idx]
      # reset indices
      rownames(df_train) <- NULL
      rownames(df_test) <- NULL
      
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
        y_train <- cbind(t(rep(1, dim(df_train_bad)[1])), 
                         t(rep(0, dim(df_train_bad)[1])))
        y_train <- drop(y_train)
        rownames(df_train) <- NULL
      }
      
      # determine number of missing values that should be amputed
      nr_mv <- floor(nrow(df_train) * ncol(df_train) * (percent_mv / 100))
      # save original train data
      df_train_orig <- df_train
      
      # artificially remove data points
      if (mv_mode=="MCAR") {
        counter <- 0
        while (!sum(is.na(df_train))==nr_mv) {
          counter <- counter + 1
          df_train <- df_train_orig
          # introduce random missing values to fill in
          df_train[!is.na(df_train)][sample(seq(df_train[!is.na(df_train)]),
                                     nr_mv)] <- NA
          # define data types
          df_train <- definedatatypes_pred(df_train)
          if ((counter%%10000)==0) {
            print("Might be caught in infinite loop.")
          }
        }
      } else if (mv_mode=="MAR") {
        counter <- 0
        # unfortunately the ampute package is not deterministic in the sense of 
        # how many values it amputes -> while loop until the desired number of 
        # values are removed
        while (!sum(is.na(df_train))==nr_mv) {
          counter <- counter + 1
          df_train <- df_train_orig
          
          freq_vec <- rep(1 / dim(df_train)[2], dim(df_train)[2])
          if (percent_mv>MAR_max) {
            tmp <- ampute(df_train, prop=MAR_max/100, bycases=FALSE,
                          freq=freq_vec, mech=mv_mode)
          } else {
            tmp <- ampute(df_train, prop=percent_mv/100, bycases=FALSE,
                          freq=freq_vec, mech=mv_mode)
          }
          df_train <- tmp$amp
          
          # fill rest with MCAR if > max MAR percent
          if (percent_mv>MAR_max) {
            nr_mv_MAR_max <- floor(nrow(df_train) * ncol(df_train) *
                                   (MAR_max / 100))
            df_train[!is.na(df_train)][sample(seq(df_train[!is.na(df_train)]),
                                       nr_mv-nr_mv_MAR_max)] <- NA
          }
          # define data types
          df_train <- definedatatypes_pred(df_train)
          
          if ((counter%%10000)==0) {
            print("Might be caught in infinite loop.")
          }
        }
      }
      
      ### IMPUTATION
      
      # 1. mean imputation
      df_train_mean <- df_train
      for(c in 1:ncol(df_train_mean)){
        df_train_mean[is.na(df_train_mean[,c]),
                      c] <- round(mean(df_train_mean[,c], na.rm=TRUE))
      }
      df_train_mean <- definedatatypes_pred(df_train_mean)
      AUC_per_percent[k,1,i] <- glm_model(df_train_mean, df_test, y_train,
                                          y_test)
      
      # 2. median imputation
      df_train_median <- df_train
      for(c in 1:ncol(df_train_median)){
        df_train_median[is.na(df_train_median[,c]),
                        c] <- round(median(df_train_median[,c], na.rm=TRUE))
      }
      df_train_median <- definedatatypes_pred(df_train_median)
      AUC_per_percent[k,2,i] <- glm_model(df_train_median, df_test, y_train,
                                          y_test)
      
      # 3. MICE
      imp <- mice(df_train, m=nr_iter_mice, print=FALSE, threshold=10)
      # impute with last iteration
      df_train_mice <- complete(imp, action=nr_iter_mice)
      df_train_mice <- definedatatypes_pred(df_train_mice)
      AUC_per_percent[k,3,i] <- glm_model(df_train_mice, df_test, y_train,
                                          y_test)

      # 4. hot-deck
      imp <- mice(df_train, method="sample", m=1, print=FALSE, threshold=10)
      df_train_hd <- complete(imp, action=1)
      df_train_hd <- definedatatypes_pred(df_train_hd)
      AUC_per_percent[k,4,i] <- glm_model(df_train_hd, df_test, y_train, y_test)
      
      # 5. EM
      imp <- amelia(df_train, m=nr_iter_em, boot.type="none", p2s=0)
      df_train_em <- imp$imputations[[nr_iter_em]]
      df_train_em <- definedatatypes_pred(df_train_em)
      AUC_per_percent[k,5,i] <- glm_model(df_train_em, df_test, y_train, y_test)
      
      # 6. listwise deletion
      if (percent_mv<=lwise_threshold) {
        df_tmp <- cbind(df_train, y_train)
        df_tmp <- na.omit(df_tmp)
        y_train_lw <- df_tmp$y_train
        df_train_lw <- subset(df_tmp, select=-c(y_train))
        rownames(df_train_lw) <- NULL
        AUC_per_percent[k,6,i] <- glm_model(df_train_lw, df_test, y_train_lw,
                                            y_test)
      }
    }
      
    # store values
    filename <- paste(mv_mode, "_AUC.RData", sep="")
    save(AUC_noimp, AUC_per_percent, 
         file=paste(path, "results/predictive-model/across-parameters/", 
                    filename, sep=""))
  }
  
  end_time <- Sys.time()
  diff <- end_time - start_time
  print(paste("Calculation time:", diff))
}