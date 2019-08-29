library(ROCR)
library(pROC)
library(glmnet)


# preprocessing function specific to the 1000Plus stroke dataset 
preprocess <- function(filepath) {
  df <- read.csv(filepath)
  
  # exclude columns with a lot of NA values + patient IDs
  df <- subset(df, select = -c(mssid, PatientID, symptomatic_bleed_during_stay,
                               NIHSS_post_treatment))
  
  # turn 9's and 8's into NA values
  df$scan1_occlusion[df$scan1_occlusion==9] <- NA
  df$scan1_occlusion[df$scan1_occlusion==8] <- NA
  df$mRS_before_stroke_discuss_inclusion[
    df$mRS_before_stroke_discuss_inclusion==9] <- NA
  df$risk_factor_hyperlipidemia[df$risk_factor_hyperlipidemia==9] <- NA
  df$risk_factor_Diabetes[df$risk_factor_Diabetes==9] <- NA
  df$risk_factor_hypertonia[df$risk_factor_hypertonia==9] <- NA
  df$TOAST_discharge_discuss_inclusion[
    df$TOAST_discharge_discuss_inclusion==9] <- NA
  
  # delete negative values for hours to mri + remove outlier
  df$HOURS.TO.MRT1[df$HOURS.TO.MRT1<0] <- NA
  df$HOURS.TO.MRT1[df$HOURS.TO.MRT1>5000] <- NA
  
  df<- transform(df,gender=as.numeric(gender))
  df$gender <- df$gender-1

  # define data types
  df <- definedatatypes(df)

  return(df)
}


# definition of datatypes per column, also specific to 1000Plus dataset
definedatatypes <- function(df) {
  df <- transform(df, HOURS.TO.MRT1=as.numeric(HOURS.TO.MRT1),
                  age=as.numeric(age),
                  gender=as.logical(gender),
                  treatment=as.logical(treatment),
                  scan1_occlusion=as.logical(scan1_occlusion),
                  mRS_before_stroke_discuss_inclusion= 
                    as.numeric(mRS_before_stroke_discuss_inclusion),
                  risk_factor_hyperlipidemia= 
                    as.logical(risk_factor_hyperlipidemia),
                  risk_factor_Diabetes=as.logical(risk_factor_Diabetes),
                  risk_factor_hypertonia=as.logical(risk_factor_hypertonia),
                  NIHSS_scan1=as.numeric(NIHSS_scan1),
                  NIHSS_scan3=as.numeric(NIHSS_scan3),
                  mRS_discharge=as.numeric(mRS_discharge),
                  TOAST_discharge_discuss_inclusion= 
                    as.numeric(TOAST_discharge_discuss_inclusion))
  return(df)
}


# definition of datatypes for parameters that are used in the predictive model
# adjusted for 1000Plus dataset
definedatatypes_pred <- function(df) {
  df <- transform(df, HOURS.TO.MRT1=as.numeric(HOURS.TO.MRT1),
                  age=as.numeric(age),
                  gender=as.logical(gender),
                  treatment=as.logical(treatment),
                  scan1_occlusion=as.logical(scan1_occlusion),
                  mRS_before_stroke_discuss_inclusion= 
                    as.numeric(mRS_before_stroke_discuss_inclusion),
                  risk_factor_hyperlipidemia=
                    as.logical(risk_factor_hyperlipidemia),
                  risk_factor_Diabetes=as.logical(risk_factor_Diabetes),
                  risk_factor_hypertonia=as.logical(risk_factor_hypertonia),
                  NIHSS_scan1=as.numeric(NIHSS_scan1))
  return(df)
}


# normalize all columns that are specified in vector column_name_vec 
normalizedata <- function(df, column_name_vec) {
  for (i in array(1:length(column_name_vec))) {
    df[[column_name_vec[i]]] <- (df[[column_name_vec[i]]] - 
      min(df[[column_name_vec[i]]])) / (max(df[[column_name_vec[i]]]) - 
      min(df[[column_name_vec[i]]]))
  }
  return(df)
}


# output path for imputations
get_output_path <- function(percent_mv, param_mode, mv_mode, nr_iterations) {
  return(paste(toString(mv_mode), "_", toString(param_mode), "_percentage_",
               toString(percent_mv), "_nriter_", toString(nr_iterations), 
               ".RData", sep=""))
}


# create coherent dataframe 
create_df <- function(data_avg, data_sd, imp_methods=imp_methods) {
  df <- data.frame()
  for (m in array(1:dim(data_avg)[2])) {
    df_tmp <- as.data.frame(data_avg[,m])
    colnames(df_tmp) <- "value"
    df_tmp["percentage"] <- seq(1,dim(data_avg)[1],1)
    df_tmp["Method"] <- imp_methods[m]
    df_tmp["sd"] <- data_sd[,m]
    df <- rbind(df,df_tmp)
  }
  return(df)
}


# general plotting function 
plot_error <- function(error_df, param, title_string, y_string, imp_methods) {
  if (plot_std) {
    if (plot_line) {
      error_plot <- ggplot(error_df, aes(x=percentage, y=value, color=Method)) + 
        geom_line(size=1) + geom_ribbon(aes(ymin=value-0.5*sd, 
                                            ymax=value+0.5*sd,fill=Method),
                                        alpha=0.1, color=NA) + 
        facet_wrap(~missingvaluemechanism) + xlab("% missing values") + 
        ylab(y_string) + 
        theme(axis.title.x=element_text(size=20),
              axis.title.y=element_text(size=20),
              axis.text.x=element_text(size=16),
              plot.title=element_text(size=22, face="bold", hjust=0.5),
              axis.text.y=element_text(size=16), 
              strip.text.x=element_text(size=16, face="bold"),
              legend.text=element_text(size=16), 
              legend.title=element_text(size=16, face="bold"))
    } else {
      error_plot <- ggplot(error_df, aes(x=percentage, y=value, color=Method)) + 
        geom_point(size=1) + geom_errorbar(aes(ymin=value-0.5*sd,
                                               ymax=value+0.5*sd,
                                               color=Method)) + 
        facet_wrap(~missingvaluemechanism) + xlab("% missing values") + 
        ylab(y_string) + 
        theme(axis.title.x=element_text(size=20),
              axis.title.y=element_text(size=20),
              axis.text.x=element_text(size=16),
              plot.title=element_text(size=22, face="bold", hjust=0.5),
              axis.text.y=element_text(size=16),
              strip.text.x=element_text(size=16, face="bold"),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16,face="bold"))
    } 
  } else {
    if (plot_line) {
      error_plot <- ggplot(error_df, aes(x=percentage, y=value, color=Method)) + 
        geom_line(size=1) + facet_wrap(~missingvaluemechanism) + 
        xlab("% missing values") + ylab(y_string) + 
        theme(axis.title.x=element_text(size=20), 
              axis.title.y=element_text(size=20),
              axis.text.x=element_text(size=16),
              plot.title=element_text(size=22, face="bold", hjust=0.5),
              axis.text.y=element_text(size=16),
              strip.text.x=element_text(size=16, face="bold"),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16, face="bold"))
    } else {
      error_plot <- ggplot(error_df, aes(x=percentage, y=value, color=Method)) + 
        geom_point(size=1) + facet_wrap(~missingvaluemechanism) + 
        xlab("% missing values") + ylab(y_string) + 
        theme(axis.title.x=element_text(size=20), 
              axis.title.y=element_text(size=20),
              axis.text.x=element_text(size=16),
              plot.title=element_text(size=22, face="bold", hjust=0.5),
              axis.text.y=element_text(size=16), 
              strip.text.x=element_text(size=16, face="bold"),
              legend.text=element_text(size=16), 
              legend.title=element_text(size=16, face="bold"))
    }
  }
  
  return(error_plot)
}


# plotting function for error measurements (error and bias)
plot_error_bias <- function(data_df, y_string) {
  error_plot <- ggplot(data_df, aes(x=percentage, y=value, color=Method)) + 
    geom_line(size=1) + geom_ribbon(aes(ymin=value-0.5*sd, ymax=value+0.5*sd,
                                        fill=Method), alpha=0.1, color=NA) + 
    facet_grid(errorcalc~missingvaluemechanism,scale="free") + 
    xlab("% missing values") + ylab(y_string) + 
    theme(axis.title.x=element_text(size=20), 
          axis.title.y=element_text(size=20),
          axis.text.x=element_text(size=16), 
          plot.title=element_text(size=22, face="bold", hjust=0.5),
          axis.text.y=element_text(size=16), 
          strip.text.x=element_text(size=16, face="bold"),
          strip.text.y=element_text(size=16, face="bold"),
          legend.text=element_text(size=16), 
          legend.title=element_text(size=16, face="bold")) 
  return(error_plot)
}


plot_perf <- function(data_df, y_string) {
  perf_plot <- ggplot(data_df, aes(x=percentage,y=value,color=Method)) + 
    geom_line(size=1) + geom_ribbon(aes(ymin=value-0.5*sd, ymax=value+0.5*sd,
                                        fill=Method), alpha=0.1, color=NA) +    
    facet_wrap(~zoomed,scale="free") + xlab("% missing values") + 
    ylab(y_string) + 
    theme(axis.title.x=element_text(size=20), 
          axis.title.y=element_text(size=20),
          axis.text.x=element_text(size=16), 
          plot.title=element_text(size=22, face="bold", hjust=0.5),
          axis.text.y=element_text(size=16), 
          strip.text.x=element_text(size=16, face="bold"),
          strip.text.y=element_text(size=16, face="bold"),
          legend.text=element_text(size=16), 
          legend.title=element_text(size=16, face="bold")) 
  return(perf_plot)
}


# plotting function for p values
plot_pvals <- function(data_df, y_string) {
  perf_plot <- ggplot(data_df, aes(x=percentage, y=value, color=Method)) + 
    geom_line(size=1) + geom_ribbon(aes(ymin=value-0.5*sd, ymax=value+0.5*sd,
                                        fill=Method), alpha=0.1, color=NA) +                   
    facet_wrap(~pval,scale="free") + xlab("% missing values") + ylab(y_string) + 
    geom_hline(yintercept=0.05, linetype="dashed", color="black", size=1) +
    theme(axis.title.x=element_text(size=20), 
          axis.title.y=element_text(size=20),
          axis.text.x=element_text(size=16), 
          plot.title=element_text(size=22, face="bold", hjust=0.5),
          axis.text.y=element_text(size=16), 
          strip.text.x=element_text(size=16, face="bold"),
          strip.text.y=element_text(size=16, face="bold"),
          legend.text=element_text(size=16), 
          legend.title=element_text(size=16, face="bold")) 
  return(perf_plot)
}


# GLM model function for given training and test set (including labels)
glm_model <- function(df_train, df_test, y_train, y_test) {
  fit_glm <- glm(y_train~., family=binomial(link="logit"), data=df_train)
  p_test <- predict(fit_glm, df_test, type="response")
  pr_test <- prediction(p_test, y_test)
  auc_test <- performance(pr_test, measure="auc")
  return(auc_test@y.values[[1]])
}
