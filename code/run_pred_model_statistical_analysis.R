# STATISTICAL ANALYSIS: for predictive model across parameters

library(ggplot2)
library(dplyr)

# get functions from helper script and variables from config script
source("helper_functions.R")
source("config.R")

result_path <- paste(path, "results/predictive-model/", sep="")

AUC <- data.frame()

# 1. compare every method for every percentage with the non-imputed data
for (mv_m in mv_modes) {
  filename <- paste(mv_m, "_AUC.RData", sep="")
  load(paste(path, "results/predictive-model/across-parameters/", filename,
             sep=""))  
  
  blank_data <- rep(NaN, dim(AUC_per_percent)[1] * length(imp_methods_lw))
  p_vals <- array(blank_data, c(dim(AUC_per_percent)[1],
                                length(imp_methods_lw)))
  for (m in array(1:length(imp_methods_lw))) {
    for (k in array(1:dim(AUC_per_percent)[1])) {
      if (sum(is.na(AUC_per_percent[k,m,]))==0) {
        p_vals[k,m] <- wilcox.test(AUC_noimp, AUC_per_percent[k,m,],
                                   alternative="greater", paired=TRUE)$p.value
      }
    }
    plot(array(1:dim(AUC_per_percent)[1]), p_vals[,m], main=imp_methods_lw[m],
         xlab="% missing values", ylab="p value", type="l")
    lines(array(1:dim(AUC_per_percent)[1]), rep(0.05,dim(AUC_per_percent)[1]), 
          col="red")
  }
  print(mv_m)
  significant <- p_vals<0.05
  for (m in array(1:length(imp_methods_lw))) {
    print(imp_methods_lw[m])
    print(paste("First significant percentage: ", min(which(significant[,m])),
                sep=""))
    print(paste("Last non-significant percentage: ", 
                max(which(!significant[,m])), sep=""))
  }
}

# 2. compare every method for each percentage with each other
for (mv_m in mv_modes) {
  filename <- paste(mv_m, "_AUC.RData", sep="")
  load(paste(path, "results/predictive-model/across-parameters/",
             filename, sep=""))  
  
  blank_data <- rep(NaN, dim(AUC_per_percent)[1] * length(imp_methods_lw) *
                    length(imp_methods_lw))
  p_vals <- array(blank_data, c(dim(AUC_per_percent)[1], length(imp_methods_lw),
                                length(imp_methods_lw)))
  for (k in array(1:dim(AUC_per_percent)[1])) {
    for (m in array(1:length(imp_methods_lw))) {
      for (n in array(1:length(imp_methods_lw))) {
        if ((sum(is.na(AUC_per_percent[k,m,]))==0) & 
            (sum(is.na(AUC_per_percent[k,n,]))==0)) {
          p_vals[k,m,n] <- wilcox.test(AUC_per_percent[k,m,],
                                       AUC_per_percent[k,n,],
                                       alternative="greater",
                                       paired=TRUE)$p.value
        }
      }
    }
  }
}