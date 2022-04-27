library(tidyverse)
library(skimr)
library(finalfit)
library(rio)


data <- read_csv("LUAD_clin.csv")
data2 <- read_csv("LUSC_clin.csv")

# Remove Not reported data

#LUAD
data <- data %>%
  mutate(across(where(is.character), ~na_if(., "not reported")))
data <- data %>%
  mutate(across(where(is.character), ~na_if(., "Not Reported")))

#LUSC
data2 <- data2 %>%
  mutate(across(where(is.character), ~na_if(., "not reported")))
data2 <- data2 %>%
  mutate(across(where(is.character), ~na_if(., "Not Reported")))

# Select variables based on NA count (> 50% complete is a good choice!)
# ToDo: function to select variables with % completeness  

#LUAD
NA_fifty <- dim(data)[1]/2

NA_sum <- colSums(is.na(data))
NA_sum <- as.data.frame(NA_sum)
NA_sum <- tibble::rownames_to_column(NA_sum, "variables")
NA_sum <- NA_sum[NA_sum$NA_sum < NA_fifty ,]

data <- subset(data, select = NA_sum$variables)

#LUSC
NA_fifty <- dim(data2)[1]/2

NA_sum <- colSums(is.na(data2))
NA_sum <- as.data.frame(NA_sum)
NA_sum <- tibble::rownames_to_column(NA_sum, "variables")
NA_sum <- NA_sum[NA_sum$NA_sum < NA_fifty ,]

data2<- subset(data2, select = NA_sum$variables)


## Saving dataset ------------------------------------

write_csv(data, path = "LUAD_clinreported.csv")
write_csv(data2, path = "LUSC_clinreported.csv")


#FinalFit LUAD

# Table 1 - Patient demographics by variable of interest ----
explanatory = c("sync_malign","tumor_stg","origin", "age_diagnose" ,"prim_diagnose", "prior_cancer","prior_treat","t_patholog", "n_patholog","m_patholog","cigar_day","pack_years","ethnicity","race" ,"gender","age", "radiation_treat")
dependent = "vital_status" 
data %>%
  summary_factorlist(dependent, explanatory,
                     p=TRUE, add_dependent_label=TRUE) -> t1
knitr::kable(t1, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

export(t1, file = "tabelaLUAD_clean.xlsx")

#FinalFit LUSC

# Table 2 - Patient demographics by variable of interest ----
explanatory = c("sync_malign","tumor_stg","origin", "age_diagnose" ,"prim_diagnose", "prior_cancer","prior_treat","t_patholog", "n_patholog","m_patholog","cigar_day","pack_years","ethnicity","race" ,"gender","age", "radiation_treat")
dependent = "vital_status" 
data2 %>%
  summary_factorlist(dependent, explanatory,
                     p=TRUE, add_dependent_label=TRUE) -> t2
knitr::kable(t2, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

export(t2, file = "tabelaLUSC_clean.xlsx")


