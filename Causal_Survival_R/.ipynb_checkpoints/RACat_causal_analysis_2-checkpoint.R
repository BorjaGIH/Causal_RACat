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
#data <- read_csv('supervivencia_75_v_prep.csv', col_types = cols(data_VIH = col_datetime()))
data <- read_csv('../racat_prep.csv')

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
#covariates <- data[c("edat", "any_qx", "bmi_val", "charlindex", "codisexe", "nivah", "Diabetes", "Obesitat", "Rheumatoid_Arthritis", "smoking_value", "Durada_intervencio_minuts", "viscositat", "Alcohol_Abuse")]
covariates <- data[c( "any_qx",  "nivah",  "Durada_intervencio_minuts", "viscositat")]
covariates <- covariates %>% mutate_all(as.numeric)
event_time <- data$T
treatment <- data$Antibiotic
eti <- 'Ei'
event_type <- data %>% pull(eti)
Y=event_time
W=treatment
D=event_type

############################## Causal Survival Forest ##########################

horizon_start <- 5
horizon_end <- 150
horizon_vec <- seq(horizon_start, horizon_end, length.out=1+(horizon_end-horizon_start)*1)
ate <- c()
ate_std <- c()
index <- c()

for (i in 1:length(horizon_vec)){
  
  horizon <- horizon_vec[i]
  failure.time <- seq(0, horizon, length.out = horizon)
  # tau(X) = P[T(1) > horizon | X = x] - P[T(0) > horizon | X = x]
  cs.forest <- causal_survival_forest(X=covariates, Y=round(event_time,0), W=treatment, 
                                          D=event_type, target="survival.probability", # "RMST" "survival.probability"
                                          failure.times=round(failure.time,0), horizon=horizon, censoring_model='forest',
                                          mtry=3)

  # E[Y(1) - Y(0)] (target.sample = all)
  ate_i <- average_treatment_effect(cs.forest)
  ate<-c(ate, ate_i[[1]])
  ate_std<-c(ate_std, ate_i[[2]])
  index<-c(index, horizon)
  
  print(paste0("time: ", round(horizon, 2), ", ate: ", round(ate_i[[1]],4)))
  
  # Save the csf to file
  #saveRDS(object, file = "object.rds")
}

ate_vector<-cbind(index,ate,ate_std)
# Save matrix of ates and sds
write.csv(ate_vector, sprintf("ate_vector_%s.csv", eti))


#ggplot(ate_vector, aes(x = index, y = ate)) +
#  geom_errorbar(aes(ymin = ate - ate_std, ymax = ate + ate_std), width = 0.2) +
#  geom_line(colour = "red") +
#  geom_point(colour = "red", size = 2) +
#  labs(x = "Index", y = "Value")
#  ggsave(sprintf("ate_vect_150_%s.pdf", eti))

################################################################################