# Imports
library(readr)
library(dplyr)
library(lubridate)
library(grf)
library(ggplot2)
library(tidyverse)
library(caret)
library(npsurv)
library(condSURV)
library(survival)
library(survminer)
library(ranger)
library(ggfortify)

theme_set(theme_bw())

################################## Preprocess ##################################

# Read data
data <- read_csv('../../racat_prep.csv', show_col_types = FALSE)

################################################################################

# Shuffle and reduce
data <- data[sample(nrow(data)),]
#data <- data[sample(nrow(data)*.1),]

# Encode cardinals
encode_ordinal <- function(x, order = unique(x)) {
  x <- as.numeric(factor(x, levels = order, exclude = NULL))
  x
}
data$nivah <- encode_ordinal(data$nivah)
data$viscositat <- encode_ordinal(data$viscositat)

# Select and format other variables
covariates <- data[c("edat", "any_qx", "bmi_val", "charlindex", "codisexe", "nivah", "Diabetes", "Obesitat", "Rheumatoid_Arthritis", "smoking_value", "Durada_intervencio_minuts", "viscositat", "Alcohol_Abuse")]
covariates <- covariates %>% mutate_all(as.numeric)
event_time <- data$T
treatment <- data$Antibiotic
eti <- 'Ea'
event_type <- data %>% pull(eti)
Y=event_time
W=treatment
D=event_type
continuous_covars<-c("edat", "charlindex", "bmi_val", "Durada_intervencio_minuts")
categorical_covars<-c("nivah", "smoking_value", "any_qx", "viscositat")
binary_covars<-c( "codisexe", "Diabetes", "Obesitat", "Rheumatoid_Arthritis", "Alcohol_Abuse")


########################## CSF #################################

horizon_start <- 5
horizon_end <- 120
horizon_vec <- seq(horizon_start, horizon_end, length.out=1+(horizon_end-horizon_start)*0.2)

for (covariate in colnames(covariates)){
  assign(paste0(covariate,"_mat"), list())
}

for (horizon in horizon_vec){
  
  # fit model
  failure.time <- seq(0, horizon, length.out = horizon)
  cs.forest <- causal_survival_forest(X=covariates, Y=round(event_time,0), W=treatment, 
                                      D=event_type, target="survival.probability", # "RMST" "survival.probability"
                                      failure.times=round(failure.time,0), horizon=horizon, censoring_model='forest',
                                      mtry=3) 
  
  # predict: simulate all t=1 - t=0 (OOB)
  predicted <- predict(cs.forest, NULL, estimate.variance = TRUE)
  
  
  for (covariate in colnames(covariates)){
    
    if (covariate=='edat' | covariate=='bmi_val' | covariate=='Durada_intervencio_minuts'){
      
      breaks <- seq(min(covariates[[covariate]]), max(covariates[[covariate]]), length.out = 21)
      binned_covar <- cut(covariates[[covariate]], breaks = breaks)
      pred_and_covar <- data.frame(cbind(predicted[[1]], predicted[[2]], binned_covar))
      colnames(pred_and_covar) <- c("prediction", "std", "binned_covar")
      group_means <- pred_and_covar %>%
        group_by(binned_covar) %>% 
        summarise(cate = mean(prediction), cate_std = mean(std))
            
      assign(paste0(covariate,"_mat"), cbind(get(paste0(covariate,"_mat")), paste(group_means[['cate']], group_means[['cate_std']], sep=';')))
      
      # Save index for each covar
      assign(paste0(covariate,"_ind"), levels(binned_covar[group_means[['binned_covar']]]))
            
    } else {
     
      pred_and_covar <- data.frame(cbind(predicted[[1]], predicted[[2]], covariates[[covariate]]))
      colnames(pred_and_covar) <- c("prediction", "std", "covar")
      group_means <- pred_and_covar %>% 
        group_by(covar) %>% 
        summarise(cate = mean(prediction), cate_std = mean(std))
      
      assign(paste0(covariate,"_mat"), cbind(get(paste0(covariate,"_mat")), paste(group_means[['cate']], group_means[['cate_std']],sep=';')))
      
      # Save index for each covar
      assign(paste0(covariate,"_ind"), group_means[['covar']])
    }
  }
}

# Set index and columns and save
for (covariate in colnames(covariates)){
  data_matrix <- cbind(get(paste0(covariate,"_ind")),get(paste0(covariate,"_mat")))
  colnames(data_matrix) <- c(0,horizon_vec)
  print(data_matrix)
  if (covariate=='edat' | covariate=='bmi_val' | covariate=='Durada_intervencio_minuts'){
    write.csv(data_matrix, sprintf("results_%s_wci/ate_matrix_%s_%s.csv", eti, covariate, eti), quote=1)
  }
  else{
    write.csv(data_matrix, sprintf("results_%s_wci/ate_matrix_%s_%s.csv", eti, covariate, eti))
  }
}

