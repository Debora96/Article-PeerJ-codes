install.packages("survminer")
install.packages("skimr")

library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(skimr)
library(tidyverse)
library(ggplot2)


luad_clin <- GDCquery_clinic(project = "TCGA-LUAD", type = "Clinical")
lusc_clin <- GDCquery_clinic(project = "TCGA-LUSC", type = "Clinical")

dataLuad <- (luad_clin)

# Preprocessing LUAD dataset 

dataLuad <- dataLuad %>%
  mutate(across(where(is.character), ~na_if(., "not reported")))
dataLuad <- dataLuad %>%
  mutate(across(where(is.character), ~na_if(., "Not Reported")))
dataLuad <- dataLuad %>%
  mutate(across(where(is.character), ~na_if(., 'american indian or alaska native')))
dataLuad <- dataLuad %>%
  mutate(across(where(is.character), ~na_if(., "asian")))
  
glimpse(dataLuad)
summary(dataLuad)

dataLuad %>%
  select_if(is.character) %>%
  skim()

dataLuad %>%
  select_if(is.numeric) %>%
  skim()

dataLuad <-  subset(dataLuad, select= -c(last_known_disease_status, created_datetime, state, classification_of_tumor,tumor_grade, progression_or_recurrence, alcohol_history, treatments_pharmaceutical_treatment_type, treatments_radiation_treatment_type, disease) ) 
dataLuad <-  subset(dataLuad, select= -c(year_of_death, alcohol_intensity, days_to_last_known_disease_status, treatments_pharmaceutical_days_to_treatment_end, treatments_pharmaceutical_days_to_treatment_start, treatments_pharmaceutical_regimen_or_line_of_therapy, treatments_pharmaceutical_therapeutic_agents, treatments_pharmaceutical_treatment_effect, treatments_pharmaceutical_initial_disease_status, treatments_pharmaceutical_treatment_intent_type, treatments_pharmaceutical_treatment_anatomic_site, treatments_pharmaceutical_treatment_outcome, treatments_radiation_days_to_treatment_end, treatments_radiation_days_to_treatment_start,  treatments_radiation_regimen_or_line_of_therapy,  treatments_radiation_therapeutic_agents, treatments_radiation_treatment_effect, treatments_radiation_initial_disease_status, treatments_radiation_treatment_intent_type, treatments_radiation_treatment_anatomic_site, treatments_radiation_treatment_outcome) )
dataLuad <-  subset(dataLuad, select= -c(diagnosis_id, exposure_id, demographic_id, treatments_pharmaceutical_treatment_id, treatments_radiation_treatment_id, bcr_patient_barcode, days_to_recurrence  ) )
dataLuad <-  subset(dataLuad, select= -c(updated_datetime, ajcc_staging_system_edition, NA.  ) )
dataLuad <-  subset(dataLuad, select= -c(days_to_diagnosis, prim_diagnose, year_of_diagnosis, cigar_day, years_smoked, pack_years, days_to_birth, age_at_diagnosis, year_of_birth ) )

dataLuad <- dataLuad %>%
  rename(sync_malign = 'synchronous_malignancy',
         tumor_stg ='ajcc_pathologic_stage',
         origin = 'tissue_or_organ_of_origin',
         prior_cancer = 'prior_malignancy',
         t_patholog = 'ajcc_pathologic_t',
         n_patholog = 'ajcc_pathologic_n',
         m_patholog = 'ajcc_pathologic_m',
         phamaceutic_treat= 'treatments_pharmaceutical_treatment_or_therapy',
         radiation_treat= 'treatments_radiation_treatment_or_therapy')

dataLuad <- dataLuad %>%
  mutate(tumor_stg = fct_collapse(tumor_stg,
                                  T1 = c('Stage I', 'Stage IA', 'Stage IB'),
                                  T2 = c('Stage II', 'Stage IIA', 'Stage IIB'),
                                  others = c('Stage IIIA', 'Stage IIIB','Stage IV')))
                      
dataLuad <- dataLuad %>%
  mutate(origin = fct_collapse(origin,
                               Lobes = c('Lower lobe, lung', 'Middle lobe, lung', 'Upper lobe, lung')))

# Change vital status to binary values
dataLuad$vital_status <- factor(dataLuad$vital_status, 
                     levels = c("Dead", "Alive"), 
                     labels = c("1", "0"))

# Division of patients into age groups
hist(dataLuad$age_at_index)
hist(dataLuad$time, labels = TRUE,  xlim = c(0,7248), ylim = c(0,200))
dataLuad <- dataLuad %>% mutate(age_at_index = ifelse(age_at_index >=66, "old", "young"))
dataLuad$age_at_index <- factor(dataLuad$age_at_index)

# Creating the time variable
dataLuad$time <- dataLuad$days_to_death
dataLuad$time[is.na(dataLuad$days_to_death)] <- dataLuad$days_to_last_follow_up[is.na(dataLuad$days_to_death)]

#Remove NAs of the time variable
dataLuad <- dataLuad[complete.cases(dataLuad[ , c('time')]), ] 

#Changing object classes

dataLuad$vital_status <- as.numeric(dataLuad$vital_status)
dataLuad$gender <- as.factor(dataLuad$gender)

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = dataLuad$time, event = dataLuad$vital_status)
surv_object 

# Building the survival graphics

fit1 <- survfit(surv_object ~ tumor_stg, data = dataLuad)
summary(fit1)
ggsurvplot(fit1, data = dataLuad, pval = TRUE)

fit2 <- survfit(surv_object ~ race, data = dataLuad)
summary(fit2)

ggsurvplot(fit2, data = dataLuad, pval = TRUE)

fit3 <- survfit(surv_object ~ gender, data = dataLuad)
summary(fit3)

ggsurvplot(fit3, data = dataLuad, pval = TRUE)

fit4 <- survfit(surv_object ~ ethnicity, data = dataLuad)
summary(fit4)

ggsurvplot(fit4, data = dataLuad, pval = TRUE)

fit5 <- survfit(surv_object ~ age_at_index, data = dataLuad)
summary(fit5)

ggsurvplot(fit5, data = dataLuad, pval = TRUE)

# Cox proportional hazards models LUAD (Figure 3A)

fit.coxph <- coxph(surv_object ~ tumor_stg + race + gender+ ethnicity+ age_at_index, 
                   data = dataLuad)
ggforest(fit.coxph, data = dataLuad)

export(dataLuad, file = "dataLuad.csv")


# Same Analysis for LUSC

dataLusc <- (lusc_clin)

# Preprocessing LUAD dataset 

dataLusc <- dataLusc %>%
  mutate(across(where(is.character), ~na_if(., "not reported")))
dataLusc <- dataLusc %>%
  mutate(across(where(is.character), ~na_if(., "Not Reported")))
dataLusc <- dataLusc %>%
  mutate(across(where(is.character), ~na_if(., 'american indian or alaska native')))
dataLusc <- dataLusc %>%
  mutate(across(where(is.character), ~na_if(., "asian")))

glimpse(dataLusc)
summary(dataLusc)

dataLusc %>%
  select_if(is.character) %>%
  skim()

dataLusc %>%
  select_if(is.numeric) %>%
  skim()

dataLusc <-  subset(dataLusc, select= -c(last_known_disease_status, created_datetime, state, classification_of_tumor,tumor_grade, progression_or_recurrence, alcohol_history, treatments_pharmaceutical_treatment_type, treatments_radiation_treatment_type, disease) ) 
dataLusc <-  subset(dataLusc, select= -c(year_of_death, alcohol_intensity, days_to_last_known_disease_status, treatments_pharmaceutical_days_to_treatment_end, treatments_pharmaceutical_days_to_treatment_start, treatments_pharmaceutical_regimen_or_line_of_therapy, treatments_pharmaceutical_therapeutic_agents, treatments_pharmaceutical_treatment_effect, treatments_pharmaceutical_initial_disease_status, treatments_pharmaceutical_treatment_intent_type, treatments_pharmaceutical_treatment_anatomic_site, treatments_pharmaceutical_treatment_outcome, treatments_radiation_days_to_treatment_end, treatments_radiation_days_to_treatment_start,  treatments_radiation_regimen_or_line_of_therapy,  treatments_radiation_therapeutic_agents, treatments_radiation_treatment_effect, treatments_radiation_initial_disease_status, treatments_radiation_treatment_intent_type, treatments_radiation_treatment_anatomic_site, treatments_radiation_treatment_outcome) )
dataLusc <-  subset(dataLusc, select= -c(diagnosis_id, exposure_id, demographic_id, treatments_pharmaceutical_treatment_id, treatments_radiation_treatment_id, bcr_patient_barcode, days_to_recurrence  ) )
dataLusc <-  subset(dataLusc, select= -c(updated_datetime, ajcc_staging_system_edition, NA.  ) )
dataLusc <-  subset(dataLusc, select= -c(days_to_diagnosis, prim_diagnose, year_of_diagnosis, cigar_day, years_smoked, pack_years, days_to_birth, age_at_diagnosis, year_of_birth ) )

dataLusc <- dataLusc %>%
  rename(sync_malign = 'synchronous_malignancy',
         tumor_stg ='ajcc_pathologic_stage',
         origin = 'tissue_or_organ_of_origin',
         prior_cancer = 'prior_malignancy',
         t_patholog = 'ajcc_pathologic_t',
         n_patholog = 'ajcc_pathologic_n',
         m_patholog = 'ajcc_pathologic_m',
         phamaceutic_treat= 'treatments_pharmaceutical_treatment_or_therapy',
         radiation_treat= 'treatments_radiation_treatment_or_therapy')

dataLusc <- dataLusc %>%
  mutate(tumor_stg = fct_collapse(tumor_stg,
                                  T1 = c('Stage I', 'Stage IA', 'Stage IB'),
                                  T2 = c('Stage II', 'Stage IIA', 'Stage IIB'),
                                  others = c('Stage III','Stage IIIA', 'Stage IIIB','Stage IV')))

dataLusc <- dataLusc %>%
  mutate(origin = fct_collapse(origin,
                               Lobes = c('Lower lobe, lung', 'Middle lobe, lung', 'Upper lobe, lung')))


# Change vital status to binary values
dataLusc$vital_status <- factor(dataLusc$vital_status, 
                                levels = c("Dead", "Alive"), 
                                labels = c("1", "0"))

dataLusc$race <- factor(dataLusc$race, 
                                levels = c("black or african amreican", "White"), 
                                levels = c("Negros", "Brancos"))

# Division of patients into age groups
hist(dataLusc$age_at_index)
hist(dataLusc$time, labels = TRUE,  xlim = c(0,5287), ylim = c(0,200))

dataLusc <- dataLusc %>% mutate(age_at_index = ifelse(age_at_index >=68, "old", "young"))
dataLusc$age_at_index <- factor(dataLusc$age_at_index)

# Creating the time variable
dataLusc$time <- dataLusc$days_to_death
dataLusc$time[is.na(dataLusc$days_to_death)] <- dataLusc$days_to_last_follow_up[is.na(dataLusc$days_to_death)]


#Remove NAs of the time variable
dataLusc <- dataLusc[complete.cases(dataLusc[ , c('time')]), ] 

#Changing object classes

dataLusc$vital_status <- as.numeric(dataLusc$vital_status)
dataLusc$gender <- as.factor(dataLusc$gender)

# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = dataLusc$time, event = dataLusc$vital_status)
surv_object 

# Building the survival graphics

fit1 <- survfit(surv_object ~ tumor_stg, data = dataLusc)
summary(fit1)
ggsurvplot(fit1, data = dataLusc, pval = TRUE)


fit2 <- survfit(surv_object ~ race, data = dataLusc)
summary(fit2)

ggsurvplot(fit2, data = dataLusc, pval = TRUE)

fit3 <- survfit(surv_object ~ gender, data = dataLuscfilter)
summary(fit3)

ggsurvplot(fit3, data = dataLuscfilter, pval = TRUE)

fit4 <- survfit(surv_object ~ ethnicity, data = dataLuscfilter)
summary(fit4)

ggsurvplot(fit4, data = dataLuscfilter, pval = TRUE)

fit5 <- survfit(surv_object ~ age_at_index, data = dataLuscfilter)
summary(fit5)

ggsurvplot(fit5, data = dataLuscfilter, pval = TRUE)

# Cox proportional hazards models LUAD (Figure 3B)

fit.coxph <- coxph(surv_object ~ tumor_stg + race + gender+ ethnicity+ age_at_index, 
                   data = dataLusc)
ggforest(fit.coxph, data = dataLusc)

export(dataLusc, file = "dataLusc.csv")



#Survival LUAD vs LUSC

clinical <- plyr::rbind.fill(luad_clin,lusc_clin)

glimpse(clinical)
summary(clinical)
View(clinical)

clinical %>%
  select_if(is.character) %>%
  skim()

clinical %>%
  select_if(is.numeric) %>%
  skim()

# Change vital status to binary values
clinical$vital_status <- factor(clinical$vital_status, 
                                levels = c("Dead", "Alive"), 
                                labels = c("1", "0"))

# Creating Time
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

#Remove NAs 
clinical <- clinical[complete.cases(clinical[ , c('time')]), ] 

#Changing classes 

clinical$vital_status <- as.numeric(clinical$vital_status)


# Fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = clinical$time, event = clinical$vital_status)
surv_object 


fit1 <- survfit(surv_object ~ disease, data = clinical)
summary(fit1)
ggsurvplot(fit1, data = clinical, pval = TRUE)

# Survival Figure 2
ggsurvplot(
  fit1,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,7500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  # in legend of risk table
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
 
)
