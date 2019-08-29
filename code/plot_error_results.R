# PLOT: script to plot error and absolute bias results

library(ggplot2)
library(dplyr)

# get functions from helper script and variable from config script
source("helper_functions.R")
source("config.R")

param_modes <- c("num", "cat") 

for (param_mode in param_modes) {
  
  rmse_avg_df <- data.frame()
  bias_avg_df <- data.frame()
  rmse_acc_df <- data.frame()
  bias_acc_df <- data.frame()
  
  for (i in array(1:length(mv_modes))) {
    
    # load one of the four result files
    filename <- (paste("results_", toString(mv_modes[i]), "_", param_mode,
                       "_percentage_", toString(max(percent_mv_vec)),
                       "_nriter_", toString(nr_iterations), ".RData", sep=""))
    load(paste(path, "results/error-analysis/same-samples/", filename,
               sep=""))
    load(paste(path, "results/error-analysis/different-samples/", filename,
               sep=""))
    
    # create dataframe for each error measurement
    rmse_avg_tmp <- create_df(rmse_avg,rmse_sd, imp_methods)
    bias_avg_tmp <- create_df(bias_avg,bias_sd, imp_methods)
    rmse_acc_tmp <- create_df(rmse_acc,rmse_sd, imp_methods)
    bias_acc_tmp <- create_df(bias_acc,bias_sd, imp_methods)
    
    # add missing value mechanism MCAR or MAR
    rmse_avg_tmp["missingvaluemechanism"] <- mv_modes[i]
    bias_avg_tmp["missingvaluemechanism"] <- mv_modes[i]
    rmse_acc_tmp["missingvaluemechanism"] <- mv_modes[i]
    bias_acc_tmp["missingvaluemechanism"] <- mv_modes[i]
    
    # merge datasets
    rmse_avg_df <- rbind(rmse_avg_df, rmse_avg_tmp)
    bias_avg_df <- rbind(bias_avg_df, bias_avg_tmp)
    rmse_acc_df <- rbind(rmse_acc_df, rmse_acc_tmp)
    bias_acc_df <- rbind(bias_acc_df, bias_acc_tmp)
  }
  
  rmse_avg_df["errorcalc"] <- "Average"
  rmse_acc_df["errorcalc"] <- "Accumulated average"
  bias_avg_df["errorcalc"] <- "Average"
  bias_acc_df["errorcalc"] <- "Accumulated average"
  rmse <- rbind(rmse_avg_df, rmse_acc_df)
  bias <- rbind(bias_avg_df, bias_acc_df)
  
  rmse <- rmse %>% filter(percentage<=plot_max_percentage)
  bias <- bias %>% filter(percentage<=plot_max_percentage)
  rmse <- rmse %>% filter(Method!="Median")
  bias <- bias %>% filter(Method!="Median")
  
  if (plot_every_5) {
    rmse <- rmse %>% filter(((percentage %% 5)==0) | (percentage==1))
    bias <- bias %>% filter(((percentage %% 5)==0) | (percentage==1))
  }
  
  if (plot_up_to_25) {
    rmse <- rmse %>% filter(((percentage<26)))
    bias <- bias %>% filter(((percentage<26)))
  }
  
  if (param_mode=="num") {
    rmse_plot_num <- plot_error_bias(rmse, "Root mean squared error")
    bias_plot_num <- plot_error_bias(bias, "Absolute bias")
  } else {
    rmse_plot_cat <- plot_error_bias(rmse, "% of misclassified samples")
    bias_plot_cat <- plot_error_bias(bias, "Absolute bias")
  }

}

path_to_plot <- paste(path_to_plots, "error-analysis/", sep="")
ggsave(paste(path_to_plot, "rmse_plot_num.png", sep=""), plot=rmse_plot_num,
       device="png", width=10, height=6)
ggsave(paste(path_to_plot, "rmse_plot_cat.png", sep=""), plot=rmse_plot_cat,
       device="png", width=10, height=6)

ggsave(paste(path_to_plot, "bias_plot_num.png", sep=""), plot=bias_plot_num,
       device="png", width=10, height=6)
ggsave(paste(path_to_plot, "bias_plot_cat.png", sep=""), plot=bias_plot_cat,
       device="png", width=10, height=6)