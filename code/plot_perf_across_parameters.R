# PLOT: plot performance of GLM across parameters

library(ggplot2)
library(dplyr)

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

result_path <- paste(path, "results/predictive-model/", sep="")

for (mv_m in mv_modes) {
  
  AUC <- data.frame()
  
  filename <- paste(mv_m, "_AUC.RData", sep="")
  load(paste(path, "results/predictive-model/across-parameters/", filename,
             sep=""))  
  
  AUC_mean <- apply(AUC_per_percent[,,], c(1,2), mean)
  AUC_sd <- apply(AUC_per_percent[,,], c(1,2), sd)
  AUC_tmp <- create_df(AUC_mean, AUC_sd, imp_methods_lw)
  
  AUC_tmp["missingvaluemechanism"] <- mv_m

  AUC <- rbind(AUC, AUC_tmp)
  
  # add no imputation AUC value to each method
  AUC_noimp_mean <- mean(AUC_noimp)
  AUC_noimp_sd <- sd(AUC_noimp)
  for (m in imp_methods_lw) {
    AUC[nrow(AUC)+1,] = list(value=AUC_noimp_mean, percentage=0, Method=m, 
                             sd=AUC_noimp_sd, missingvaluemechanism=mv_m)
  }
  
  AUC <- AUC %>% filter(Method!="Median")
  
  if (without_listwise) {
    AUC <- AUC %>% filter((Method!="Listwise deletion"))
    #acc <- acc %>% filter((method!="lwise del"))
    AUC_plot <- plot_error(AUC, param_str, param_str, "AUC", imp_methods_lw)
  } else {
    AUC <- AUC %>% filter(!((percentage>lwise_threshold) & 
                              (Method=="Listwise deletion")))

    AUC_zoomed <- AUC
    
    AUC["zoomed"] <- "Performance up to 70% missing values"
    
    AUC_zoomed <- AUC_zoomed %>% filter(percentage<=lwise_threshold)
    AUC_zoomed["zoomed"] <- "Performance zoomed in"
    
    AUC_all <- rbind(AUC, AUC_zoomed)
    
    AUC_all <- AUC_all %>% filter((percentage<=plot_max_percentage))
    #AUC <- AUC %>% filter(((percentage %% 5)==0) | 
    #         ((percentage<lwise_threshold)&((percentage %% 2)==0)))
    AUC_all <- AUC_all %>% filter(((percentage %% 5)==0) | 
                                    (percentage<lwise_threshold))
    
    AUC_all["pval"] <- FALSE
    
    # get p values for plot (in comparison of best method mean in our case)
    imp_methods_plot <- imp_methods_lw[-2]
    AUC_per_percent <- AUC_per_percent[,-2,]
    p_vals_mean <- data.frame(matrix(ncol=7, nrow=0))
    colnames(p_vals_mean) <- c("percentage", "Method", "missingvaluemechanism",
                               "value", "zoomed", "pval", "sd")
    for (m in array(2:(length(imp_methods_lw)-1))) {
      # fill percentage 0 with 1s
      p_vals_mean[nrow(p_vals_mean)+1,] = list(percentage=0, 
                                               Method=imp_methods_plot[m], 
                                               missingvaluemechanism=mv_m, 
                                               value=1, 
                                               zoomed=paste("Performance up to",
                                                      " 70% missing values",
                                                      sep=""),
                                               pval=paste("Mean imputation vs.",
                                                    " other methods", sep=""),
                                               sd=NA)
      for (k in array(1:dim(AUC_per_percent)[1])) {
        if (sum(is.na(AUC_per_percent[k,m,]))==0) {
          p_val <- wilcox.test(AUC_per_percent[k,1,], AUC_per_percent[k,m,],
                               alternative="greater", paired=TRUE)$p.value
          p_vals_mean[nrow(p_vals_mean)+1,] = list(percentage=k, 
                                                   Method=imp_methods_plot[m], 
                                                   missingvaluemechanism=mv_m, 
                                                   value=p_val, 
                                                   zoomed=paste("Performance ",
                                                          "up to 70% missing ",
                                                          "values", sep=""), 
                                                   pval=paste("Mean imputation",
                                                        " vs. other methods", 
                                                        sep=""), 
                                                   sd=NA)
        }
      }
    }
    p_vals_mean <- p_vals_mean %>% filter(!((percentage>lwise_threshold) & 
                                              (Method=="Listwise deletion")))
    p_vals_mean <- p_vals_mean %>% filter((percentage %% 5)==0)
    
    p_vals_lw <- data.frame(matrix(ncol=7, nrow=0))
    colnames(p_vals_lw) <- c("percentage", "Method", "missingvaluemechanism",
                             "value", "zoomed", "pval", "sd")
    for (m in array(1:(length(imp_methods_lw)-2))) {
      p_vals_lw[nrow(p_vals_lw)+1,] = list(percentage=0, 
                                           Method=imp_methods_plot[m], 
                                           missingvaluemechanism=mv_m, 
                                           value=1, 
                                           zoomed="Performance zoomed in", 
                                           pval=paste("Other methods vs. ",
                                                "listwise deletion",sep=""),
                                           sd=NA)
      for (k in array(1:lwise_threshold)) {
        if (sum(is.na(AUC_per_percent[k,m,]))==0) {
          p_val <- wilcox.test(AUC_per_percent[k,m,], AUC_per_percent[k,5,],
                               alternative="greater", paired=TRUE)$p.value
          p_vals_lw[nrow(p_vals_lw) + 1,] = list(percentage=k, 
                                                 Method=imp_methods_plot[m], 
                                                 missingvaluemechanism=mv_m, 
                                                 value=p_val, 
                                                 zoomed="Performance zoomed in", 
                                                 pval=paste("Other methods vs.",
                                                      " listwise deletion",
                                                      sep=""),
                                                 sd=NA)
        }
      }
    }
    
    p_vals_all <- rbind(p_vals_mean, p_vals_lw)
    p_plot <- plot_pvals(p_vals_all, "p value")
    AUC_plot <- plot_perf(AUC_all,"AUC")
    
    both_plots <- grid.arrange(AUC_plot, p_plot, nrow=2, heights=c(5,2),
                               widths=17)
  }
  
  path_to_plot <- paste(path_to_plots, "predictive-model/", sep="")
  if (without_listwise) {
    ggsave(paste(path_to_plot, "without-listwise/AUC_plot_across_params_", mv_m,
                 ".png", sep=""),
           plot=AUC_plot, device="png", width=10, height=10)
  } else {
    ggsave(paste(path_to_plot, "with-listwise/AUC_and_p_plot_across_params_",
                 mv_m, ".png", sep=""), plot=both_plots, device="png", width=15,
           height=10)
    ggsave(paste(path_to_plot, "with-listwise/AUC_plot_across_params_", mv_m,
                 ".png", sep=""), plot=AUC_plot, device="png", width=15,
           height=10)
    ggsave(paste(path_to_plot, "with-listwise/p_plot_across_params_", mv_m,
                 ".png", sep=""), plot=p_plot, device="png", width=15, height=3)
  }
}