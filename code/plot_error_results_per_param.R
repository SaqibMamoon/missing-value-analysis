# PLOT: plot error per parameter

library(ggplot2)
library(dplyr)

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

if (param_mode=="num") {
  column_name_vec <- column_name_num
} else {
  column_name_vec <- column_name_cat
}

# per param
for (p in array(1:length(column_name_vec))) {
  
  rmse <- data.frame()
  bias <- data.frame()

  for (i in array(1:length(mv_modes))) {
    
    # load one of the four result files
    filename <- (paste("results_", toString(mv_modes[i]), "_",
                       toString(param_mode), "_percentage_", 
                       toString(max(percent_mv_vec)), "_nriter_",
                       toString(nr_iterations), "_per_param.RData", sep=""))

    load(paste(path, "results/error-analysis/different-samples/",
               filename, sep=""))
    load(paste(path, "results/error-analysis/same-samples/", filename,
               sep=""))
    
    # create dataframe for each error measurement
    rmse_tmp <- create_df(rmse_avg[,p,], rmse_sd[,p,])
    bias_tmp <- create_df(bias_avg[,p,], bias_sd[,p,])
  
    # add missing value mechanism MCAR or MAR
    rmse_tmp["missingvaluemechanism"] <- mv_modes[i]
    bias_tmp["missingvaluemechanism"] <- mv_modes[i]

    # merge datasets
    rmse <- rbind(rmse, rmse_tmp)
    bias <- rbind(bias, bias_tmp)
  }
  
  if (plot_every_5) {
    rmse <- rmse %>% filter((percentage %% 5)==0)
    bias <- bias %>% filter((percentage %% 5)==0)
  }
  
  # plot for each error measurement
  plot_path <- paste(path, "figures/error-analysis/per-parameter/",
                     column_name_vec[p], "_rmse_num.jpg", sep="")
  rmse_plot_num <- plot_error(rmse,param_modes[j], "rmse", "RMSE")
  ggsave(plot_path, rmse_plot_num, width=8, height=5, dpi=300, units="in")
  plot_path <- paste(path, "figures/error-analysis/per-parameter/",
                     column_name_vec[p], "_bias_num.jpg", sep="")
  bias_plot_num <- plot_error(bias,param_modes[j], "bias", "mean absolute bias")
  ggsave(plot_path, bias_plot_num, width=8, height=5, dpi=300, units="in")
}




