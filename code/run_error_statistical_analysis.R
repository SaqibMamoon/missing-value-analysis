# STATISTICAL ANALYSIS: for error measurement RMSE

library(ggplot2)
library(dplyr)

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

param_modes <- c("num", "cat")

blank_data <- rep(NaN, length(mv_modes) * length(param_modes) * 
                  length(percent_mv_vec) * length(imp_methods) *
                  length(imp_methods))
p_vals_rmse <- array(blank_data, c(length(mv_modes), length(param_modes),
                                   length(percent_mv_vec), length(imp_methods),
                                   length(imp_methods)))

for (mv_m in array(1:length(mv_modes))) {
  for (p in array(1:length(param_modes))) {
    
    if (param_modes[p]=="num") {
      column_name_vec <- column_name_num
    } else {
      column_name_vec <- column_name_cat
    }
    
    for (k in array(1:length(percent_mv_vec))) {
      
      percent_mv <- percent_mv_vec[k]
      
      # load results
      filename_output <- get_output_path(percent_mv, param_modes[p], 
                                         mv_modes[mv_m], nr_iterations)
      load(paste(path, "imputations/different-samples/", filename_output,
                 sep=""))
      
      print(paste("Percentage: ", percent_mv, "/", length(percent_mv_vec)))
      
      nr_mv <- dim(idx_imp_all)[2]
      
      # storage for errors
      blank_data <- rep(NaN, length(column_name_vec) * length(imp_methods) *
                        nr_iterations * nr_mv)
      if (param_modes[p]=="num") {
        rmse_all <- array(blank_data, c(length(column_name_num),
                                        length(imp_methods), nr_mv,
                                        nr_iterations))
      } else {
        rmse_all <- array(blank_data, c(length(column_name_cat), 
                                        length(imp_methods), nr_mv,
                                        nr_iterations))
      }
      # for each imputation method
      for (m in array(1:length(imp_methods))) {
        label_rep <- rep(label_all, nr_iterations)
        if (param_modes[p]=="num") {
          label_rep <- array(label_rep, c(length(column_name_num), nr_mv,
                                          nr_iterations))
        } else {
          label_rep <- array(label_rep, c(length(column_name_cat), nr_mv,
                                          nr_iterations))
        }
        rmse_all[,m,,] <- sqrt((label_rep - imputed_values[,m,,])**2) 
      }
      # mean
      rmse_mean <- apply(rmse_all, c(2,4), mean)
      for (m in array(1:length(imp_methods))) {
        for (n in array(1:length(imp_methods))) {
          p_vals_rmse[mv_m,p,k,m,n] <- wilcox.test(rmse_mean[m,], rmse_mean[n,],
                                                   alternative="less",
                                                   paired=TRUE)$p.value
        }
      }
    }
  }
}