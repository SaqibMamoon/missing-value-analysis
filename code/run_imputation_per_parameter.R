# IMPUTE: run this script to get imputations for each parameter X times

# imports
library(lattice)
library(mice)
library(Amelia)
library(dplyr)

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

if (param_mode=="num") {
  column_name_vec <- column_name_num
} else {
  column_name_vec <- column_name_cat
}


################################################################################
#################  PART I: Impute same missing values X times  ################# 
################################################################################

# preprocessing of data, specific to data set
df <- preprocess(data_path)
# exclude all non-complete cases
df_cc <- na.omit(df)

# normalize numerical data
if (param_mode=="num") {
  df_norm <- normalizedata(df_cc, column_name_vec)
} else {
  df_norm <- df_cc
}

for (mv_mode in mv_modes) {
  for (k in array(1:length(percent_mv_vec))) {
    
    percent_mv <- percent_mv_vec[k]
    
    print(paste("Percentage: ", percent_mv, "/", length(percent_mv_vec)))
    
    nr_mv <- floor(nrow(df_norm) * (percent_mv / 100))
    
    # define variables for storing the imputed values, indices and labels
    blank_data <- rep(NaN, length(column_name_vec) * length(imp_methods) *
                           nr_iterations * nr_mv)
    imputed_values <- array(blank_data, c(length(column_name_vec),
                                          length(imp_methods),nr_mv,
                                          nr_iterations))
    blank_data <- rep(NaN, length(column_name_vec) * nr_mv)
    label_all <- array(blank_data, c(length(column_name_vec), nr_mv)) 
    idx_imp_all <- array(blank_data, c(length(column_name_vec), nr_mv)) 
    
    filename <- get_output_path(percent_mv, param_mode, mv_mode, nr_iterations)
    
    if (!file.exists(paste(path, "imputations/same-samples/", filename, 
                           sep=""))) {
      
      start_time <- Sys.time()
      
      # loop over all parameters
      for (j in array(1:length(column_name_vec))) {
        
        column_name <- column_name_vec[j]
        
        print(paste("Parameter: ", column_name, " ", j, "/",
                    length(column_name_vec)))
        
        # get datatype of parameter
        datatype <- sapply(df_norm[[column_name]], class)[1]
        
        # artificially remove data points
        if (mv_mode=="MCAR") {
          # introduce random missing values to fill in
          idx_imp_all[j,] <- sample(1:nrow(df_norm), nr_mv)
          # store true labels
          label_all[j,] <- df_norm[[column_name]][idx_imp_all[j,]]
          df_norm[[column_name]][idx_imp_all[j,]] <- NA
          #na_value_idx <- is.na(df_norm[[column_name]])
        } else if (mv_mode=="MAR") {
          counter <- 0
          # unfortunately the ampute package is not deterministic in the sense 
          # that it does not always impute the same amount of values -> while 
          # loop until the desired number of values are removed
          while (!sum(is.na(df_norm))==nr_mv) {
            counter <- counter + 1
            freq_vec <- rep(0, dim(df)[2])
            freq_vec[grep(column_name, colnames(df_norm))] <- 1
            tmp <- ampute(df_norm, prop=percent_mv/100, bycases=TRUE,
                          freq=freq_vec)
            df_norm <- tmp$amp
            # define data types
            df_norm <- definedatatypes(df_norm)
            # check if nr_mv is correct
            if (!sum(is.na(df_norm[[column_name]]))==nr_mv) {
              df_cc <- na.omit(df)
              # normalize numerical data
              if (param_mode=="num") {
                df_norm <- normalizedata(df_cc, column_name_vec)
              } else {
                df_norm = df_cc
              }
            } else {
              # store true labels
              idx_imp_all[j,] <- which(is.na(df_norm[[column_name]]))
              label_all[j,] <- tmp$data[[column_name]][idx_imp_all[j,]]
            }
            if (counter==10000) {
              print("Might be caught in infinite loop.")
            }
          }
        }
  
        # IMPUTATION
        
        # 1. & 2. mean and median imputation (same for every iteration)
        if (datatype=="numeric") {
          # calculate the imputed values
          imputed_values[j,1,,] <- rep(vector(mode="numeric", length=nr_mv) +
                                       mean(df_norm[[column_name]], na.rm=TRUE),
                                       each=nr_iterations)
          imputed_values[j,2,,] <- rep(vector(mode="numeric", length=nr_mv) +
                                       median(df_norm[[column_name]], 
                                              na.rm=TRUE),
                                       each=nr_iterations)
        } else if (datatype=="logical") {
          imputed_values[j,1,,] <- rep(vector(mode="numeric", length=nr_mv) +
                                       round(mean(df_norm[[column_name]],
                                                  na.rm=TRUE)),
                                       each=nr_iterations)
          imputed_values[j,2,,] <- rep(vector(mode="numeric", length=nr_mv) + 
                                       median(df_norm[[column_name]], 
                                              na.rm=TRUE),
                                       each=nr_iterations)
        } else if (datatype=="factor") {
          df_norm_tmp <- as.numeric(df_norm[[column_name]])
          imputed_values[j,1,,] <- rep(vector(mode="numeric", length=nr_mv) +
                                       round(mean(df_norm_tmp, na.rm=TRUE)),
                                       each=nr_iterations)
          imputed_values[j,2,,] <- rep(vector(mode="numeric", length=nr_mv) +
                                       median(df_norm_tmp, na.rm=TRUE),
                                       each=nr_iterations)
        } 
        
        for (i in array(1:nr_iterations)) {
  
          # 3. mice
          imp <- mice(df_norm, m=nr_iter_mice, print=FALSE)
          imputed_values[j,3,,i] <- imp$imp[[column_name]][,nr_iter_mice]
          
          # 4. sample
          imp_sample <- mice(df_norm, m=nr_iter_mice, method="sample",
                             print=FALSE)
          imputed_values[j,4,,i] <- imp_sample$imp[[column_name]][,nr_iter_mice]
          
          # 5. EM
          imp <- amelia(df_norm, m=nr_iter_em, boot.type="none", p2s=0)
          df_tmp <- imp$imputations[[nr_iter_mice]]
          df_tmp <- definedatatypes(df_tmp)
          imputed_values[j,5,,i] <- df_tmp[[column_name]][idx_imp_all[j,]]  
        }
        
        df_cc <- na.omit(df)
        
        # normalize numerical data
        if (param_mode=="num") {
          df_norm <- normalizedata(df_cc, column_name_vec)
        } else {
          df_norm <- df_cc
        }
      }
      
      # store values
      save(imputed_values, label_all, idx_imp_all,
           file=paste(path, "imputations/same-samples/", filename, sep=""))
      
      end_time <- Sys.time()
      diff <- end_time - start_time
      print(paste("Calculation time:", diff))
    }
  }
}


################################################################################
##############  PART II: Impute different missing values X times  ############## 
################################################################################

# preprocessing of data
df <- preprocess(data_path)

# exclude all non-complete cases
df_cc <- na.omit(df)

# normalize numerical data
if (param_mode=="num") {
  df_norm <- normalizedata(df_cc, column_name_vec)
} else {
  df_norm <- df_cc
}

for (mv_mode in mv_modes) {
  for (k in array(1:length(percent_mv_vec))) {
    
    percent_mv <- percent_mv_vec[k]
    
    print(paste("Percentage: ", percent_mv, "/", length(percent_mv_vec)))
    
    nr_mv <- floor(nrow(df_norm) * (percent_mv / 100))
    
    # define variables for storing the imputed values, indices and labels
    blank_data <- rep(NaN, length(column_name_vec) * length(imp_methods) *
                      nr_iterations * nr_mv)
    imputed_values <- array(blank_data, c(length(column_name_vec),
                                          length(imp_methods), nr_mv,
                                          nr_iterations))
    blank_data <- rep(NaN, length(column_name_vec) * nr_mv * nr_iterations)
    label_all <- array(blank_data, c(length(column_name_vec), nr_mv,
                                     nr_iterations)) 
    idx_imp_all <- array(blank_data, c(length(column_name_vec), nr_mv,
                                       nr_iterations)) 
    
    filename <- get_output_path(percent_mv, param_mode, mv_mode, nr_iterations)
    
    if (!file.exists(paste(path, "imputations/different-samples/", filename,
                           sep=""))) {
      
      start_time <- Sys.time()
      
      # loop over all parameters
      for (j in array(1:length(column_name_vec))) {
        
        column_name <- column_name_vec[j]
        
        print(paste("Parameter: ", column_name, " ", j, "/",
                    length(column_name_vec)))
        
        # get datatype of parameter
        datatype <- sapply(df_norm[[column_name]], class)[1]
        
        for (i in array(1:nr_iterations)) {
        
          # artificially remove data points
          if (mv_mode=="MCAR") {
            # introduce random missing values to fill in
            idx_imp_all[j,,i] <- sample(1:nrow(df_norm), nr_mv)
            # store true labels
            label_all[j,,i] <- df_norm[[column_name]][idx_imp_all[j,,i]]
            df_norm[[column_name]][idx_imp_all[j,,i]] <- NA
          } else if (mv_mode=="MAR") {
            counter <- 0
            while (!sum(is.na(df_norm))==nr_mv) {
              counter <- counter + 1
              freq_vec <- rep(0, dim(df)[2])
              freq_vec[grep(column_name, colnames(df_norm))] <- 1
              tmp <- ampute(df_norm, prop=percent_mv/100, bycases=TRUE,
                            freq=freq_vec)
              df_norm <- tmp$amp
              # define data types
              df_norm <- definedatatypes(df_norm)
              # check if nr_mv is correct
              if (!sum(is.na(df_norm[[column_name]]))==nr_mv) {
                df_cc <- na.omit(df)
                # normalize numerical data
                if (param_mode=="num") {
                  df_norm <- normalizedata(df_cc, column_name_vec)
                } else {
                  df_norm <- df_cc
                }
              } else {
                # store true labels
                idx_imp_all[j,,i] <- which(is.na(df_norm[[column_name]]))
                label_all[j,,i] <- tmp$data[[column_name]][idx_imp_all[j,,i]]
              }
              if (counter==100) {
                print("Might be caught in infinite loop.")
              }
            }
          }
        
          # IMPUTATION
          
          # 1. & 2. mean and median imputation (same for every iteration)
          if (datatype=="numeric") {
            # calculate the imputed values
            imputed_values[j,1,,i] <- mean(df_norm[[column_name]], na.rm=TRUE)
            imputed_values[j,2,,i] <- median(df_norm[[column_name]], na.rm=TRUE)
            #print(imputed_values[j,1,,])
          } else if (datatype=="logical") {
            imputed_values[j,1,,] <- round(mean(df_norm[[column_name]],
                                                na.rm=TRUE))
            imputed_values[j,2,,] <- median(df_norm[[column_name]], na.rm=TRUE)
          } else if (datatype=="factor") {
            df_norm_tmp <- as.numeric(df_norm[[column_name]])
            imputed_values[j,1,,] <- round(mean(df_norm_tmp, na.rm=TRUE))
            imputed_values[j,2,,] <- median(df_norm_tmp, na.rm=TRUE)
          } 
          
          # 3. mice
          imp <- mice(df_norm, m=nr_iter_mice, print=FALSE)
          imputed_values[j,3,,i] <- imp$imp[[column_name]][,nr_iter_mice]
          
          # 4. sample
          imp_sample <- mice(df_norm, m=nr_iter_mice, method="sample",
                             print=FALSE)
          imputed_values[j,4,,i] <- imp_sample$imp[[column_name]][,nr_iter_mice]
  
          # 5. EM
          imp <- amelia(df_norm, m=nr_iter_mice, boot.type="none", p2s=0)
          df_tmp <- imp$imputations[[nr_iter_mice]]
          df_tmp <- definedatatypes(df_tmp)
          imputed_values[j,5,,i] <- df_tmp[[column_name]][idx_imp_all[j,,i]]  
          
          df_cc <- na.omit(df)
          
          # normalize numerical data
          if (param_mode=="num") {
            df_norm <- normalizedata(df_cc, column_name_vec)
          } else {
            df_norm = df_cc
          }
        }
      }
      
      # store values
      save(imputed_values, label_all, idx_imp_all, 
           file=paste(path, "imputations/different-samples/", filename, sep=""))
      
      end_time <- Sys.time()
      diff <- end_time - start_time
      print(paste("Calculation time:", diff))
    }
  }
}