###### TCGA Clinical Data Analysis #######
# cleaning and exploring dataset 

install.packages("tidyverse", dependencies = TRUE)
install.packages("skimr", dependencies = TRUE)
install.packages("finalfit", dependencies = TRUE)
BiocManager::install("GenomicDataCommons")
BiocManager::install("TCGAbiolinks")

library(tidyverse)
library(skimr)
library(finalfit)
library(GenomicDataCommons)
library(TCGAbiolinks)

vignette('skimr')

## 1. Data importing and visualizing ---------------------------
# Download clinical data at GDC (select primary site, project and choose 'clinical' at exploration page).
LUADclinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical", save.csv = TRUE)

# Using read_csv or read_delim, na = c("", "NA")
LUADclinicalraw <- read_delim("TCGA-LUAD_clinical.tsv", "\t", 
                            na = c("", "NA"),
                            escape_double = FALSE, 
                            trim_ws = TRUE)

class(LUADclinicalraw) 
dim(LUADclinicalraw) 
names(LUADclinicalraw) 
glimpse(LUADclinicalraw)
skim(LUADclinicalraw) 
#View(LUADclinicalraw)

## 2. Cleaning data ---------------------------
# Select variables based on NA count (> 50% complete is a good choice!)
# ToDo: function to select variables with % completeness  

NA_fifty <- dim(LUADclinicalraw)[1]/2

NA_sum <- colSums(is.na(LUADclinicalraw))
NA_sum <- as.data.frame(NA_sum)
NA_sum <- tibble::rownames_to_column(NA_sum, "variables")
NA_sum <- NA_sum[NA_sum$NA_sum < NA_fifty ,]

LUAD_clean <- subset(LUADclinicalraw, select = NA_sum$variables)


## Remove duplicate observations
LUAD_clean0 <- LUAD_clean %>%
  distinct_at('submitter_id', .keep_all = TRUE)

## Remove numeric variables with unique observations  
# ToDo: function to select and remove variables with unique observations? 
LUAD_clean0 %>%
  select_if(is.numeric) %>%
  skim()

LUAD_clean1 <-  subset(LUAD_clean0, select= -c(days_to_diagnosis) ) 

# Remove character variables with unique observations 
LUAD_clean1 %>%
  select_if(is.character) %>%
  skim()

LUAD_clean2 <-  subset(LUAD_clean1, select= -c(last_known_disease_status, state, classification_of_tumor,tumor_grade, progression_or_recurrence, alcohol_history, treatments_pharmaceutical_treatment_type,treatments_radiation_treatment_type, disease) )

# Remove character variables with similar information - check each one!

LUAD_clean3 <-  subset(LUAD_clean2, select= -c(ajcc_pathologic_stage, site_of_resection_or_biopsy, bcr_patient_barcode) )

# Remove other variables not directly related to patient - check each one!

LUAD_clean4 <-  subset(LUAD_clean3, select= -c(updated_datetime, ajcc_staging_system_edition, icd_10_code, diagnosis_id, exposure_id, demographic_id, treatments_pharmaceutical_treatment_id, treatments_radiation_treatment_id  ) )
View(LUAD_clean4)

## 3. Changing variables names ---------------------------
# Using snake_style 

LUAD_clean4 <- LUAD_clean4 %>%
  rename(sync_malign = 'synchronous_malignancy',
         tumor_stg ='tumor_stage',
         origin = 'tissue_or_organ_of_origin',
         last_follow_up ='days_to_last_follow_up',
         age_diagnose = 'age_at_diagnosis',
         prim_diagnose = 'primary_diagnosis',
         prior_cancer = 'prior_malignancy',
         year_diagnose = 'year_of_diagnosis',
         prior_treat = 'prior_treatment',
         t_patholog = 'ajcc_pathologic_t',
         morphology = 'morphology',
         n_patholog = 'ajcc_pathologic_n',
         m_patholog = 'ajcc_pathologic_m',
         cigar_day = 'cigarettes_per_day',
         pack_years = 'pack_years_smoked',
         age = 'age_at_index',
         year_birth = 'year_of_birth',
         phamaceutic_treat= 'treatments_pharmaceutical_treatment_or_therapy',
         radiation_treat= 'treatments_radiation_treatment_or_therapy')
view(LUAD_clean4)

## 4. Taming data ------------------------------------------
# Use pck lubridate for dates
LUAD_clean4 <- LUAD_clean4 %>%
  mutate_if(is.character, as.factor) %>%
  mutate(patient_id = as.character(patient_id),
         age = as.integer(age),
         year_diagnose = as.integer(year_diagnose))


## 5. Checking NA patterns -----------------------------
# Check distincts types of NAs: MCAR, MAR, MNAR

LUAD_clean4  %>%
missing_plot()
missing_glimpse(LUAD_clean4)

LUAD_clean5 <- LUAD_clean4[!is.na(LUAD_clean4$tumor_stg),]

LUAD_clean5  %>%
  missing_plot()

missing_glimpse(LUAD_clean5)

## 6. Checking numeric variables -----------------------------
# Check data distribution, unplausible values, outliers.
# Never delete an unusual value if the value is a possible one. 
# Deleting unusual values will bias your results and cause you to underestimate the variability in the observations.

# Describing numeric variables with summary().
# If the median and mean are similar, the distribution is likely roughly symmetrical. 
# Otherwise, it will be skewed to the right or to the left.
LUAD_clean5 %>%
  select_if(is.numeric) %>%
  summary()

# Histograms or density plots
ggplot(LUAD_clean5, aes(age)) +
  geom_histogram(bins = 20, alpha = 0.8, color = "red")

ggplot(LUAD_clean5, aes(age_diagnose)) +
  geom_histogram(bins = 20, alpha = 0.8, color = "red")

ggplot(LUAD_clean5, aes(year_diagnose)) +
  geom_density(color = "red")

ggplot(LUAD_clean5, aes(year_birth)) +
  geom_density(color = "red")

ggplot(LUAD_clean5, aes(days_to_birth)) +
  geom_density(color = "red")

ggplot(LUAD_clean5, aes(pack_years)) +
  geom_density(color = "red")

ggplot(LUAD_clean5, aes(last_follow_up)) +
  geom_density(color = "red")

# Boxplots 
# Inter quartil range (IQR) = Q3 â€” Q1
# whiskers = Â±1.58 IQR / âˆšn âˆ— IQR, where â€˜nâ€™ = samples
# Outliers = values below or above min and max whiskers values, respectively

ggplot(LUAD_clean5, aes(x ='', y=cigar_day)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(LUAD_clean5$cigar_day)

## 7. Checking categorical variables --------------------------
# Check frequency, lables and levels 
# Cancer staging: https://www.cancer.gov/about-cancer/diagnosis-staging/staging

LUAD_clean5 %>%
  select_if(is.factor) %>%
  summary() 

# To do: fct_drop to drop unused levels

# agregating levels
LUAD_clin <- LUAD_clean5 %>%
  mutate(tumor_stg = fct_collapse(tumor_stg,
                                  I = c('stage i', 'stage ia', 'stage ib'),
                                  II = c('stage ii', 'stage iia', 'stage iib'),
                                  III = c('stage iiia', 'stage iiib'),
                                  IV= c('stage iv')))

LUAD_clin <- LUAD_clin %>%
  mutate(t_patholog = fct_collapse(t_patholog, 
                                   T1 = c('T1', 'T1a', 'T1b'),
                                   T2 = c('T2', 'T2a', 'T2b')))

LUAD_clin <- LUAD_clin %>%
  mutate(m_patholog = fct_collapse(m_patholog, 
                                   M1 = c('M1', 'M1a', 'M1b')))

# changing level names
LUAD_clin <- LUAD_clin %>%
  mutate(ethnicity = fct_recode(ethnicity, 'hispanic/latino'='hispanic or latino', 'not hispanic/latino'='not hispanic or latino'),
         race = fct_recode(race, 'American.indian/Alaska.native'='american indian or alaska native', 'Black/African.american'='black or african american'))

LUAD_clin %>%
  select_if(is.factor) %>%
  summary()

# to visualize

LUAD_clin %>% 
  ggplot(aes(x = origin, fill = origin)) + 
  geom_bar() +
  theme_bw(15) +
  xlab("tissue of origin") + 
  ylab("frequency") + 
  theme(legend.position = 'none')

LUAD_clin %>% 
  ggplot(aes(x = origin, y = age, fill = origin)) + 
  geom_boxplot() +
  theme_bw(15) +
  xlab("tissue of origin") + 
  ylab("age") + 
  facet_wrap(~ gender) +
  theme(legend.position = 'none')


## 8. Checking again -----------------------

skim(LUAD_clin)
view(LUAD_clin)


## 9. Saving dataset ------------------------------------
write_csv(LUAD_clin, path = "your_path_here/LUAD_clin.csv")

#rm(LUAD_clean3, LUAD_clean2,LUAD_clean1, LUAD_clean0, LUAD_clean, LUADclinicalraw, NA_sum, NA_fifty)
