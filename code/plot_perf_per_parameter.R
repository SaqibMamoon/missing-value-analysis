# PLOT: plot performance of GLM per parameter

library(ggplot2)
library(dplyr)

path = "/home/tabea/Documents/Code/outcome-prediction/missing-data/"

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

result_path <- paste(path, "results/predictive-model/", sep="")

mv_modes <- c("MCAR", "MAR")

for (param in 1:length(column_names_imp)) {
  
  AUC <- data.frame()
  acc <- data.frame()
  
  param_str <- column_names_imp[param]
  print(param_str)
  
  for (mv_m in mv_modes) {
    
    load(paste(result_path, "per-parameter/GLM_", mv_m, "_all", ".RData",
               sep=""))
    
    # AUC
    
    AUC_mean <- apply(AUC_all[,param,,], c(1,2), mean)
    AUC_sd <- apply(AUC_all[,param,,], c(1,2), sd)
    AUC_tmp <- create_df(AUC_mean, AUC_sd)
  
    AUC_tmp["missingvaluemechanism"] <- mv_m
    
    AUC <- rbind(AUC, AUC_tmp)
    
    # accuracy
    
    acc_mean <- apply(acc_all[,param,,], c(1,2), mean)
    acc_sd <- apply(acc_all[,param,,], c(1,2), sd)
    acc_tmp <- create_df(acc_mean, acc_sd)
    
    acc_tmp["missingvaluemechanism"] <- mv_m
    
    acc <- rbind(acc, acc_tmp)
  }
  
  if (plot_every_5) {
    AUC <- AUC %>% filter(((percentage %% 5)==0))
    acc <- acc %>% filter(((percentage %% 5)==0))
  }
  
  if (without_listwise) {
    AUC <- AUC %>% filter((method!="lwise del"))
    acc <- acc %>% filter((method!="lwise del"))
  }
  
  AUC_plot <- plot_error(AUC, param_str, param_str, "AUC", imp_methods_lw)
  acc_plot <- plot_error(acc, param_str, param_str, "Accuracy", imp_methods_lw)
  
  path_to_plot <- paste(path_to_plots, "predictive-model/", sep="")
  if (without_listwise) {
    ggsave(paste(path_to_plot, "without-listwise/AUC_plot_", param_str, ".png",
                 sep=""), plot=AUC_plot, device="png", width=10)
    ggsave(paste(path_to_plot, "without-listwise/acc_plot_", param_str, ".png",
                 sep=""), plot=acc_plot, device="png", width=10)
  } else {
    ggsave(paste(path_to_plot, "with-listwise/AUC_plot_", param_str, ".png",
                 sep=""), plot=AUC_plot, device="png", width=10)
    ggsave(paste(path_to_plot, "with-listwise/acc_plot_", param_str, ".png",
                 sep=""), plot=acc_plot, device="png", width=10)
  }
}