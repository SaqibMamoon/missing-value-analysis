# EVALUATE: script to evaluate imputation methods using RMSE and bias

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

if (param_mode=="num") {
  column_name_vec <- column_name_num
} else {
  column_name_vec <- column_name_cat
}

# storage for means and standard deviations
if (per_param) {
  blank_data <- rep(NaN, length(percent_mv_vec) * length(imp_methods) *
                    length(column_name_vec))
  rmse_avg <- array(blank_data, c(length(percent_mv_vec), 
                                  length(column_name_vec), length(imp_methods)))
  rmse_acc <- array(blank_data, c(length(percent_mv_vec),
                                  length(column_name_vec), length(imp_methods)))
  rmse_sd <- array(blank_data, c(length(percent_mv_vec),
                                 length(column_name_vec), length(imp_methods)))
  bias_avg <- array(blank_data, c(length(percent_mv_vec), 
                                  length(column_name_vec), length(imp_methods)))
  bias_acc <- array(blank_data, c(length(percent_mv_vec),
                                  length(column_name_vec), length(imp_methods)))
  bias_sd <- array(blank_data, c(length(percent_mv_vec), 
                                 length(column_name_vec), length(imp_methods)))
} else {
  blank_data <- rep(NaN, length(percent_mv_vec) * length(imp_methods))
  rmse_avg <- array(blank_data, c(length(percent_mv_vec), length(imp_methods)))
  rmse_acc <- array(blank_data, c(length(percent_mv_vec), length(imp_methods)))
  rmse_sd <- array(blank_data, c(length(percent_mv_vec), length(imp_methods)))
  bias_avg <- array(blank_data, c(length(percent_mv_vec), length(imp_methods)))
  bias_acc <- array(blank_data, c(length(percent_mv_vec), length(imp_methods)))
  bias_sd <- array(blank_data, c(length(percent_mv_vec), length(imp_methods)))
}
  
for (mv_mode in mv_modes) {
  if (per_param) {
    filename <- paste("results_", toString(mv_mode), "_", toString(param_mode),
                      "_percentage_", toString(max(percent_mv_vec)),
                      "_nriter_", toString(nr_iterations), "_per_param.RData",
                      sep="")
  } else {
    filename <- paste("results_", toString(mv_mode), "_", toString(param_mode),
                      "_percentage_", toString(max(percent_mv_vec)), "_nriter_",
                      toString(nr_iterations), ".RData", sep="")
  }
  
  # ABSOLUTE BIAS
  
  if (!file.exists(paste(path, "results/error-analysis/same-samples/", filename,
                         sep=""))) {
    
    for (k in array(1:length(percent_mv_vec))) {
      
      percent_mv <- percent_mv_vec[k]
      
      # load results
      filename_output <- get_output_path(percent_mv, param_mode, mv_mode,
                                         nr_iterations)
      load(paste(path, "imputations/same-samples/", filename_output, sep=""))
      
      print(paste("Percentage: ", percent_mv, "/", length(percent_mv_vec)))
      
      nr_mv <- dim(idx_imp_all)[2]
      
      blank_data <- rep(NaN, length(column_name_vec) * 
                        length(imp_methods) * nr_mv)
      bias_all <- array(blank_data, c(length(column_name_vec),
                                      length(imp_methods), nr_mv))
      # calculate average estimates per parameter per methods and per value 
      # missing
      avg_estimates <- apply(imputed_values, c(1,2,3), mean)
      avg_tmp <- apply(avg_estimates, c(2,3), mean)
      # for each imputation method
      for (m in array(1:length(imp_methods))) {
        bias_all[,m,] <- abs(avg_estimates[,m,] - label_all)
      }
      # store in mean and sd in variables
      if (per_param) {
        bias_avg[k,,] <- apply(bias_all, c(1,2), mean)
        bias_sd[k,,] <- apply(bias_all, c(1,2), sd)  
      } else {
        bias_avg[k,] <- apply(bias_all, c(2), mean)
        # the variance is the variance of the estimated samples
        bias_sd[k,] <- apply(avg_tmp, c(1), sd)
        bias_acc[k,] <- apply(bias_all, c(2), sum)
      }  
    }
  
    # store everything
    save(bias_avg, bias_sd, bias_acc, 
         file=paste(path, "results/error-analysis/same-samples/", filename, 
                    sep=""))
  }
}

# RMSE

for (mv_mode in mv_modes) {
  if (per_param) {
    filename <- paste("results_", toString(mv_mode), "_", toString(param_mode),
                      "_percentage_", toString(max(percent_mv_vec)), 
                      "_nriter_", toString(nr_iterations), "_per_param.RData",
                      sep="")
  } else {
    filename <- paste("results_", toString(mv_mode), "_", toString(param_mode),
                      "_percentage_", toString(max(percent_mv_vec)), "_nriter_",
                      toString(nr_iterations), ".RData", sep="")
  }
  
  
  if (!file.exists(paste(path, "results/error-analysis/different-samples/",
                         filename, sep=""))) {
    
    for (k in array(1:length(percent_mv_vec))) {
      
      percent_mv <- percent_mv_vec[k]
      
      # load results
      filename_output <- get_output_path(percent_mv, param_mode, mv_mode,
                                         nr_iterations)
      load(paste(path, "imputations/different-samples/", filename_output,
                 sep=""))
      
      print(paste("Percentage: ", percent_mv, "/", length(percent_mv_vec)))
      
      nr_mv <- dim(idx_imp_all)[2]
      
      # 1. RMSE
      # storage for errors
      blank_data <- rep(NaN, length(column_name_vec) * length(imp_methods) *
                        nr_iterations * nr_mv)
      rmse_all <- array(blank_data, c(length(column_name_vec),
                                      length(imp_methods), nr_mv,
                                      nr_iterations))
      # for each imputation method
      for (m in array(1:length(imp_methods))) {
        label_rep <- rep(label_all, nr_iterations)
        label_rep <- array(label_rep, c(length(column_name_vec), nr_mv, 
                                        nr_iterations))
        rmse_all[,m,,] <- sqrt((label_rep - imputed_values[,m,,])**2) 
      }
      # store in mean and sd in variables
      if (per_param) {
        rmse_avg[k,,] <- apply(rmse_all, c(1,2), mean)
        rmse_sd[k,,] <- apply(rmse_all, c(1,2), sd)      
      } else {
        rmse_avg[k,] <- apply(rmse_all, c(2), mean)
        rmse_tmp <- apply(rmse_all, c(2,4), mean)
        rmse_acc[k,] <- apply(rmse_all, c(2), sum)
        rmse_sd[k,] <- apply(rmse_tmp, c(1), sd)
      }
    }
    
    # store everything
    save(rmse_avg, rmse_sd, rmse_acc, 
         file=paste(path, "results/error-analysis/different-samples/", filename,
                    sep=""))
  }
}