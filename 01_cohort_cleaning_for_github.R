# 0_Description ---------------------------------------------------------------------------

# Migrants' primary care utilisation before and during the COVID-19 pandemic in England: An interrupted time series
# Cohort cleaning script 
# Date started: 25/03/2020
# Author(s): Claire Zhang / Yamina Boukari / Neha Pathak / Parth Patel
# QC (date): Yamina Boukari (25/02/2022)

# 1_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc, data.table)

# 2_Set working directory ------------------------------------------------------------------

setwd("filepath")


############################ CLEAN FULL ELIGIBLE COHORT ##################### ----------------------------------

# 3_Import full cohort data ----------------------------------------------------------------------------

# Cohort files
# Jan 2021 GOLD build
all_patients <- read_csv("filepath") 

# Events files
migration <- read_csv("filepath.csv")
ethnicity <- read.table("filepath.csv", sep= '\t', header= TRUE)

# 4_Tidy full cohort data  ------------------------------------------------------------------------------
## 4a_Cohort ------------------------------------------------------------------

# Check for duplicates (whole row and patid only)
n_distinct(all_patients) == count(all_patients)
n_distinct(distinct(all_patients, patid)) == count(all_patients)

# Check new variables and classes
glimpse(all_patients)

tidy_all_patients <-  all_patients

# Change variables to the correct data type based on CPRD GOLD data specification & mapping
tidy_all_patients$prac_region <- factor(tidy_all_patients$prac_region, levels = c(1:13), labels = c("North East",
                                                                                                    "North West", "Yorkshire & The Humber", "East Midlands", "West Midlands", "East of England",
                                                                                                    "South West", "South Central", "London", "South East Coast", "Northern Ireland", "Scotland", "Wales"))
tidy_all_patients$gender <- factor(tidy_all_patients$gender, levels = c(0,1,2,3,4),
                                   labels = c("Data Not Entered", "Male", "Female", "Indeterminate", "Unknown"))
tidy_all_patients$toreason <- factor(tidy_all_patients$toreason, levels = c(0:34), labels = c("Data Not Entered",
                                                                                              "Death", "Removal to new TP/HB/CSA", "Internal Transfer", "Mental Hospital", "Embarkation",
                                                                                              "New TP/HB/CSA/Same GP", "Adopted Child", "Services", "Deduction at GP's Request",
                                                                                              "Registration Cancelled","Service Dependant",  "Deduction at Patient's Request", "Other reason",
                                                                                              "Enlistment", "Institution", "Transfer within Practice", "Linkage", "Untraced - Miscellaneous",
                                                                                              "Untraced - Immig", "Untraced - GP Resign", "Untraced - College", "Untraced - outwith Practice",
                                                                                              "Untraced - outwith HB", "Multiple Transfer", "Intra-consortium transfer", "Returned Undelivered",
                                                                                              "Internal Transfer - Address Change", "Internal Transfer within Partnership",
                                                                                              "Correspondence states gone away", "Practice advise outside their area",
                                                                                              "Practice advise patient no longer resident", "Practive advise removal via screening system",
                                                                                              "Practice advise removal via vaccination data", "Removal from Residential Institute"))
# Check factors correctly labelled
levels(tidy_all_patients$prac_region)
levels(tidy_all_patients$gender)
levels(tidy_all_patients$toreason)

# Find missing values
any(is.na(tidy_all_patients)) 
sum(is.na(tidy_all_patients)) 
summary(is.na(tidy_all_patients)) # No NAs for variables we need, therefore do not remove

# Look at summary of all variables for outliers and obvious errors
summary(tidy_all_patients)

# Look at summary of categorical variables for outliers and obvious errors where did not print in whole summary
summary(tidy_all_patients$prac_region)

# Select needed variables
all_patients_clean <- select(tidy_all_patients, c(patid, pracid, prac_region, gender, dob, data_start, data_end, frd, crd, prac_uts, prac_lcd, deathdate, tod))

# all_patients_clean ---

# Relevel prac_country and prac_region
all_patients_clean <- within(all_patients_clean, prac_region <- relevel (prac_region, ref="London"))

# Create yob variable for later use
all_patients_clean <- all_patients_clean %>%
  mutate(yob=year(dob))

n_distinct(all_patients_clean) 

# Check and save
all_patients_clean <- droplevels(all_patients_clean)
levels(all_patients_clean$prac_region)
levels(all_patients_clean$gender)

save(all_patients_clean, file = "filepath/all_patients_clean.Rdata") 

rm(tidy_all_patients, all_patients)

## 4b_Migration ---------------------------------------------------------------

# Remove whole row duplicates
n_distinct(migration) == count(migration)
distinct_migration <- migration %>% distinct()

# Remove duplicates based on patid + medcode + date as likely clinical coding errors
n_distinct(distinct_migration) == count(distinct(distinct_migration, patid, eventdate, medcode, category,  .keep_all = TRUE))
distinct_migration_2 <- distinct_migration %>% distinct(patid, eventdate, medcode, category,  .keep_all = TRUE)

# Select index event for duplicates based on patid and medcode
distinct_migration_3 <- distinct_migration_2 %>% arrange(eventdate) %>% distinct(patid, medcode, .keep_all = TRUE)

# Drop unneeded variables
tidy_migration <-  select(distinct_migration_3, -c(sysdate, constype, episode, enttype,
                                                   adid, data1, data2, data3, data4, data5, data6 , data7))

# Change variables to correct data type based on CPRD GOLD data specification & mapping
tidy_migration$category <- factor(tidy_migration$category, levels = c(1,2,3,4),
                                  labels = c("Non-UK origin", "Born outside of the UK", "First/main language not English",
                                             "Visa status indicating migration"))

# Check correct coercion of variable classes
glimpse(tidy_migration)

# Check factors correctly labelled
levels(tidy_migration$category)

# Find missing values
any(is.na(tidy_migration)) 
sum(is.na(tidy_migration)) 
summary(is.na(tidy_migration)) 

# Look at summary of all variables for outliers and obvious errors
summary(tidy_migration)

# Select needed variables
migration_clean <- select(tidy_migration, c(patid, category, eventdate))

# migration_clean ---

# Rename 'category' to avoid confusion with ethnicity 'category'
migration_clean <- rename(migration_clean, migcat = category)

# Create certainty of migration variable with only 1 record per patient - Display results as Definite (1) = category 2 (birth) & 4 (visa); Probable migrant (2) = category 3 (language); Possible migrants (3)= category 1 (origin)
n_distinct(migration_clean$patid) 
migration_clean <- migration_clean %>%
  mutate(migcertainty = as.integer(migcat))

migration_clean <- migration_clean %>%
  mutate(migcertainty=replace(migcertainty, migcertainty==2, 4)) %>%
  mutate(migcertainty=replace(migcertainty, migcertainty==4, 4))

migration_clean_eventdate <- migration_clean %>%
  select(c(patid, eventdate, migcertainty)) %>%
  group_by(patid, migcertainty) %>%
  summarise(eventdate=min(eventdate))

migration_clean <- migration_clean %>%
  group_by(patid) %>%
  summarise(migcertainty=max(migcertainty))
n_distinct(migration_clean) 

migration_clean <- left_join(migration_clean, migration_clean_eventdate, 
                             by = c("patid", "migcertainty"))

migration_clean$migcertainty <- factor(migration_clean$migcertainty, levels = c(4,3,1),
                                       labels = c("Definite", "Probable", "Possible"))

# Check and save
migration_clean <- droplevels(migration_clean)
levels(migration_clean$migcertainty)

save(migration_clean, file = "filepath/migration_clean.Rdata") 

rm(distinct_migration, distinct_migration_2,distinct_migration_3,
   migration,migration_clean_eventdate, tidy_migration)

## 4c_Ethnicity ---------------------------------------------------------------

# Remove duplicate rows
n_distinct(ethnicity) == count(ethnicity)
distinct_ethnicity <- ethnicity %>% distinct()

# Remove duplicates of combined patid + medcode + date (likely clinical coding errors)
n_distinct(distinct_ethnicity) == count(distinct(distinct_ethnicity, patid, eventdate, medcode, category,  .keep_all = TRUE))
distinct_ethnicity_2 <- distinct_ethnicity %>% distinct(patid, eventdate, medcode, category,  .keep_all = TRUE)

# Remove duplicate medcodes for unique patids by selecting first recorded medcode event
distinct_ethnicity_3 <- distinct_ethnicity_2 %>% arrange(eventdate) %>% distinct(patid, medcode, .keep_all = TRUE)

tidy_ethnicity <- distinct_ethnicity_3

# Check correct coercion of variable classes
glimpse(tidy_ethnicity)

# Find missing values
any(is.na(tidy_ethnicity)) 
sum(is.na(tidy_ethnicity))  
summary(is.na(tidy_ethnicity)) 

# Look at summary of all fields for outliers and obvious errors
summary(tidy_ethnicity)

# Select needed variables
ethnicity_clean <- select(tidy_ethnicity, c(patid, category, eventdate))

ethnicity_clean$category <- factor(ethnicity_clean$category, levels = c(1:19),
                                   labels = c("British", "Irish", "Gypsy or Irish Traveller", "Other White",
                                              "Mixed White and Black Caribbean", "Mixed White and Black African",
                                              "Mixed White and Asian", "Other Mixed",
                                              "Indian", "Pakistani", "Bangladeshi", "Chinese", "Other Asian",
                                              "African", "Caribbean", "Other Black", "Arab",
                                              "Any other ethnic group","Ethnic group not specified"))


# remove patients with duplicate records
n_distinct(distinct(ethnicity_clean, patid, category)) 
ethnicity_clean <- ethnicity_clean %>% distinct(patid, category, .keep_all = TRUE)

# Drop records with ethnicity recorded as ethnic group not specified or NA (i.e. no ethnicity data available for that patient) - note: later after joining to all_patient file, ethnicity N/A's retained as 'Unknown' category
sum(ethnicity_clean$category == "Ethnic group not specified") 
sum(is.na(ethnicity_clean$category))
ethnicity_NOS <- ethnicity_clean%>%
  filter(category == "Ethnic group not specified") 
ethnicity_clean <- ethnicity_clean%>%
  filter(category != "Ethnic group not specified") %>%
  filter(!is.na(category)) 

n_distinct(ethnicity_clean$patid)  

# Keep the most recent ethnic code for those with more than one ethnicity (for those with eventdate)
recent_ethnicity_clean <- ethnicity_clean %>%
  filter(!is.na(eventdate)) 
recent_ethnicity_clean$eventdate <- as.Date(recent_ethnicity_clean$eventdate) # convert eventdate from character to date format
recent_ethnicity_clean <- recent_ethnicity_clean %>%
  group_by(patid) %>%
  slice(which.max(eventdate))
n_distinct(recent_ethnicity_clean$patid) 
nrow(recent_ethnicity_clean) == n_distinct(recent_ethnicity_clean$patid) # T

# Select most frequently occuring ethnicity for those without eventdates
no_eventdate <- ethnicity_clean %>% 
  filter(is.na(eventdate)) 
n_distinct(no_eventdate$patid) 
no_eventdate$eventdate <- as.Date(no_eventdate$eventdate) # converts eventdate from a character variable to a date in order to later bind it with recent_ethnicity_clean

mode <- function(x) {
  ux <- unique(x)
  ux [which.max(tabulate(match(x,ux)))]
}

# no_eventdate <- no_eventdate %>% 
#   group_by(patid) %>%
#   mutate(category=mode(category)) # V2 28861 - selected the most frequently occuring category, V3 0 so not applicable
# 
# no_eventdate <- no_eventdate %>%
#   distinct(patid, category) # remove duplicates
# 
# nrow(no_eventdate) == n_distinct(no_eventdate$patid) # matches

ethnicity_clean <- bind_rows(recent_ethnicity_clean, no_eventdate) 
nrow(ethnicity_clean) == n_distinct(ethnicity_clean$patid)

# # Select eventdate ethnicity and remove no eventdate ethnicity code for patients who have both
# ethnicity_clean <- ethnicity_clean %>%
#   mutate(eventdate=replace(eventdate, is.na(eventdate), as.Date("1900-01-01")))
# 
# ethnicity_clean <- ethnicity_clean %>%
#   group_by(patid) %>%
#   slice(which.max(eventdate)) 
# 
# n_distinct(ethnicity_clean$patid) # 6360847 distinct patients
# nrow(ethnicity_clean) == n_distinct(ethnicity_clean$patid) # matches, and matches original distinct patient number before selections were made
# 
# ethnicity_clean <- ethnicity_clean %>%
#   mutate(eventdate=replace(eventdate, eventdate == "1900-01-01", NA))

# If needed for calculations like Rohini's paper, remove duplicates from ethnicity_NOS and add it back in at this stage

# Rename ethnicity category label
ethnicity_clean <- rename(ethnicity_clean, ethnicat = category)

# Create 6 group ethnicity variable 
ethnicity_clean <- ethnicity_clean %>%
  mutate(ethnicat6 = ethnicat)
levels(ethnicity_clean$ethnicat6) <-  c(levels(ethnicity_clean$ethnicat6),"White British", "White non-British", "Mixed/Multiple ethnic groups", 
                                        "Asian/Asian British", "Black/African/Caribbean/Black British", "Other ethnic group" )
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "British"] <- "White British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Irish"] <- "White non-British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Gypsy or Irish Traveller"] <- "White non-British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Other White"] <- "White non-British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Mixed White and Black Caribbean"] <- "Mixed/Multiple ethnic groups"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Mixed White and Black African"] <- "Mixed/Multiple ethnic groups"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Mixed White and Asian"] <- "Mixed/Multiple ethnic groups"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Other Mixed"] <- "Mixed/Multiple ethnic groups"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Indian"] <- "Asian/Asian British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Pakistani"] <- "Asian/Asian British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Bangladeshi"] <- "Asian/Asian British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Chinese"] <- "Asian/Asian British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Other Asian"] <- "Asian/Asian British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "African"] <- "Black/African/Caribbean/Black British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Caribbean"] <- "Black/African/Caribbean/Black British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Other Black"] <- "Black/African/Caribbean/Black British"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Arab"] <- "Other ethnic group"
ethnicity_clean$ethnicat6[ethnicity_clean$ethnicat6 == "Any other ethnic group"] <- "Other ethnic group"

# Relevel ethnicat6
ethnicity_clean <- within(ethnicity_clean, ethnicat6 <- relevel(ethnicat6, ref="White British") )

# Remove eventdate variable
ethnicity_clean <- subset(ethnicity_clean, select = -eventdate)

# Check and save
ethnicity_clean <- droplevels(ethnicity_clean) 
levels(ethnicity_clean$ethnicat)
levels(ethnicity_clean$ethnicat6)

save(ethnicity_clean, file = "filepath/ethnicity_clean.Rdata") # 6385555 obs. of 3 variables

rm(list=ls())

# 5_Joins (cohort, ethnicity, migration files) -----------------------------------------------------------------

load(file = "filepath/ethnicity_clean.Rdata")
load(file = "filepath/migration_clean.Rdata")
load(file = "filepath/all_patients_clean.Rdata")

# Remove migration_clean eventdate (only needed for secular trends check)
migration_clean <- dplyr::select(migration_clean, -eventdate)

# All patient joins ---

# Join all patient cohort file to ethnicity and migration
all_patients_alldata <- left_join(all_patients_clean, migration_clean, by = c("patid" = "patid"))
all_patients_alldata <- left_join(all_patients_alldata, ethnicity_clean, by = c("patid" = "patid"))

sum(is.na(all_patients_alldata$ethnicat))  
n_distinct(all_patients_alldata$patid) 
glimpse(all_patients_alldata)
sum(is.na(all_patients_alldata)) 
summary(is.na(all_patients_alldata)) 

# create migrant vs non-migrant variable for all_patients_alldata
sum(is.na(all_patients_alldata$migcertainty))  
all_patients_alldata <- all_patients_alldata %>%
  mutate(migrant_status= migcertainty)

levels(all_patients_alldata$migrant_status) <-  c(levels(all_patients_alldata$migrant_status),"Non-migrant", "Migrant")
all_patients_alldata$migrant_status[is.na(all_patients_alldata$migrant_status)] <- "Non-migrant"
all_patients_alldata$migrant_status[all_patients_alldata$migrant_status == "Definite"] <- "Migrant"
all_patients_alldata$migrant_status[all_patients_alldata$migrant_status == "Probable"] <- "Migrant"
all_patients_alldata$migrant_status[all_patients_alldata$migrant_status == "Possible"] <- "Migrant"

sum(all_patients_alldata$migrant_status == "Non-migrant") == sum(is.na(all_patients_alldata$migcertainty))

# Replace N/As (i.e. non-migrants) for migcertainty variable with "Non-migrant" certainty category
levels(all_patients_alldata$migcertainty)  <- c(levels(all_patients_alldata$migcertainty), "Non-migrant")
all_patients_alldata <- all_patients_alldata %>%
  mutate(migcertainty=replace(migcertainty, is.na(migcertainty), "Non-migrant"))

glimpse(all_patients_alldata)

# Relevel migcertainty and migrant_status
all_patients_alldata <- within(all_patients_alldata, migcertainty <- relevel (migcertainty, ref="Non-migrant"))
all_patients_alldata <- within(all_patients_alldata, migrant_status <- relevel (migrant_status, ref="Non-migrant"))

# Drop levels 
all_patients_alldata$migrant_status <- droplevels(all_patients_alldata$migrant_status)
levels(all_patients_alldata$migrant_status)
levels(all_patients_alldata$migcertainty)

# Recategorise NAs for ethnicity into separate category
levels(all_patients_alldata$ethnicat6)  <- c(levels(all_patients_alldata$ethnicat6), "Unknown")
levels(all_patients_alldata$ethnicat)  <- c(levels(all_patients_alldata$ethnicat), "Unknown")

all_patients_alldata <- all_patients_alldata %>%
  mutate(ethnicat6=replace(ethnicat6, is.na(ethnicat6), "Unknown")) %>%
  mutate(ethnicat=replace(ethnicat, is.na(ethnicat), "Unknown"))

sum(all_patients_alldata$ethnicat6== "Unknown") 

save(all_patients_alldata, file = "filepath/all_patients_alldata.Rdata")

rm(list = ls())
gc()

# 6_Check number of practices in total eligible population ----------------------

distinct_practices_all_patients <- n_distinct(all_patients_alldata$pracid)

############ CREATE 2015-2020 ENGLAND COHORT WITH RANDOM NMs & LINKAGES ############# -----------

# 7_Load randomly sampled group of non-migrants from England data linkages -------

# Load GOLD patid linkage file with randomly sampled non-migrant patids (derived in script 00)
# Ratio 1:4 migrants:non-migrants
gold_all_sets_patid_random <- read_tsv(file = "primary_care_cleaned_files/00_patid_lists_linkage/gold_all_sets_patid_random.txt")

# Load clean full eligible cohort file
load(file = "filepath/all_patients_alldata.Rdata")

# Join 
cohort_England <- gold_all_sets_patid_random %>%
  left_join(all_patients_alldata, by = c("patid" = "patid"))

glimpse(cohort_England)

# Check each row is a distinct patient
n_distinct(cohort_England$patid) == nrow(cohort_England)

# Check 1:4 ratio of migrants:non-migrants
4*(sum(cohort_England$migrant_status == "Migrant")) == sum(cohort_England$migrant_status == "Non-migrant")

# Check number of distinct practices
distinct_practices <- n_distinct(cohort_England$pracid)

# 8_IMD linkages  ------------------------------------------------------------
## 8a_Clean raw IMD data ---------------------------------------------------

# Raw inked data
patient_imd <- read_tsv("filepath.txt")
practice_imd <- read_tsv("filepath.txt")

# patient imd ---

# Check for any duplicate based on whole row
n_distinct(patient_imd) == count(patient_imd)

# Check for any duplicates based on combined patid and pracid
n_distinct(patient_imd) == count(distinct(patient_imd, patid,pracid,  .keep_all = TRUE))

# Change variables to correct data type based on CPRD GOLD data specification & mapping
tidy_patient_imd <-  patient_imd 
#   mutate_if(is.numeric, as.integer)

# Look at summary of all fields for outliers and obvious errors
summary(tidy_patient_imd)

tidy_patient_imd$imd2015_5 <- factor(tidy_patient_imd$imd2015_5, levels = c(1:5),
                                     labels = c("IMD 1", "IMD 2", "IMD 3", "IMD 4","IMD 5" ))

# Check correct coercion of variable classes
glimpse(tidy_patient_imd)

# Check factors correctly labelled
levels(tidy_patient_imd$imd2015_5)

# Find missing values
any(is.na(tidy_patient_imd))
sum(is.na(tidy_patient_imd))
summary(is.na(tidy_patient_imd)) # 1146 no info on imd

# Select needed variables
patient_imd_clean <- tidy_patient_imd

# Practice IMD ---

# Check for any duplicate based on whole row
n_distinct(practice_imd) == count(practice_imd)

# Check for any duplicates based on practice id
n_distinct(practice_imd) == count(distinct(practice_imd, pracid,  .keep_all = TRUE))

# drop unneeded variables: imd for NI, scotland, walers (beacuse not provided by CPRD)
tidy_practice_imd <- dplyr::select(practice_imd, -c(country, ni2017_imd_5,s2016_imd_5, w2014_imd_5))

# Look at summary of all fields for outliers and obvious errors
summary(tidy_practice_imd)

# Change variables to correct data type based on CPRD GOLD data specification & mapping
tidy_practice_imd$e2015_imd_5 <- factor(tidy_practice_imd$e2015_imd_5, levels = c(1:5),
                                        labels = c("IMD 1", "IMD 2", "IMD 3", "IMD 4","IMD 5" ))

# Check correct coercion of variable classes
glimpse(tidy_practice_imd)

# Check factors correctly labelled
levels(tidy_practice_imd$e2015_imd_5)

# Find missing values
any(is.na(tidy_practice_imd)) # F
sum(is.na(tidy_practice_imd))
summary(is.na(tidy_practice_imd))

# Select needed variables
practice_imd_clean <- tidy_practice_imd

# patient and practice IMD rename factors ---

# patient_imd_clean
patient_imd_clean <- rename(patient_imd_clean, patimd = imd2015_5)

# practice_imd_clean
practice_imd_clean <- rename(practice_imd_clean, pracimd = e2015_imd_5)

# Check and save
patient_imd_clean <- droplevels(patient_imd_clean)
practice_imd_clean <- droplevels(practice_imd_clean)
levels(patient_imd_clean$patimd)
levels(practice_imd_clean$pracimd)

save(patient_imd_clean, file = "filepath/patient_imd_clean.Rdata") 
save(practice_imd_clean, file = "filepath/practice_imd_clean.Rdata") 

## 8b_Join cleaned IMD data to full cohort file ----------------------------

load(file = "filepath/patient_imd_clean.Rdata")
load(file = "filepath/practice_imd_clean.Rdata")

# Join patient_imd_clean and practice_imd_clean ---
imd_clean <- inner_join(patient_imd_clean, practice_imd_clean, by = c("pracid" = "pracid"))
glimpse(imd_clean)
sum(is.na(imd_clean))
summary(is.na(imd_clean))

cohort_England <- left_join(cohort_England, imd_clean, by = c("patid" = "patid", "pracid" = "pracid"))

# Merge IMD (patient IMD if available, otherwise practice IMD)

cohort_England <- cohort_England %>%
  mutate(imd = ifelse(!is.na(patimd), patimd, pracimd))
cohort_England$imd <- factor(cohort_England$imd, levels = c(1:5),
                                   labels = c("IMD 1", "IMD 2", "IMD 3", "IMD 4","IMD 5" ))

sum(is.na(cohort_England$imd)) # No missing IMD

# # Add "unknown" category for IMD
# levels(cohort_England$imd)  <- c(levels(cohort_England$imd), "Unknown")
# 
# cohort_England <- cohort_England %>%
#   mutate(imd=replace(imd, is.na(imd), "Unknown"))

glimpse(cohort_England)

cohort_England <- cohort_England %>%
  select(-c(hes_e_21, death_e_21,lsoa_e_21, hes_e_19, hes_e_18))

# Save for importing into acute covid cohort cleaning script 
save(cohort_England, file = "filepath/cohort_England.Rdata")

# 9_Restrict to 2015-2020 ---------------------------------------------------

# Drop patients whose data_end is before 01/01/2015 and data_start after 26/12/2020
sum(cohort_England$data_end < as.Date("2015-01-01"))
sum(cohort_England$data_start > as.Date("2020-12-26"))  

cohort_England_2015_2020 <- cohort_England %>%
  filter(data_end >= as.Date("2015-01-01") & data_start <= as.Date("2020-12-26"))

# 10_Calculate patient follow-up time over entire period   ------

cohort_England_2015_2020$data_start <- as_date(cohort_England_2015_2020$data_start)
cohort_England_2015_2020$data_end <- as_date(cohort_England_2015_2020$data_end)
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(cohort_entry = as_date(ifelse(data_start > "2014-12-31", data_start,as_date("2015-01-01")))) %>%
  mutate(cohort_exit = as_date(ifelse(data_end > "2020-12-26", as_date("2020-12-26"),data_end))) %>%
  mutate(total_pdays = (cohort_exit - cohort_entry) + 1) %>% # Add 1 day for those who are only in the cohort for 1 day, who otherwise would have had follow-up time of 0 pdays
  mutate(total_pyears = total_pdays/365) # Close approximate 365 days per year since 2x leap years but 5 days missing from 2020
cohort_England_2015_2020$total_pyears <- as.numeric(cohort_England_2015_2020$total_pyears)
cohort_England_2015_2020$total_pdays <- as.numeric(cohort_England_2015_2020$total_pdays)

glimpse(cohort_England_2015_2020)

# 11_Cohort exclusions and save final cohort file for demographics analysis ----

## Total non-progressive n,% of patients with each exclusion characteristic

# Create denominator of all patients 2015-2020
total_2015_to_2020 <- n_distinct(cohort_England_2015_2020$patid)

# data_start > data_end errors
data_start_data_end_errors <- sum(cohort_England_2015_2020$data_start > cohort_England_2015_2020$data_end) 
data_start_data_end_errors_percent <- data_start_data_end_errors/total_2015_to_2020 

# indeterminate gender
indeterminate_gender <- sum(cohort_England_2015_2020$gender == 'Indeterminate') 
indeterminate_gender_percent <- indeterminate_gender/total_2015_to_2020 

# possible migrants
possible_migrant <- sum(cohort_England_2015_2020$migcertainty == 'Possible') 
possible_migrant_percent <- possible_migrant/total_2015_to_2020 

## Exclusions made and progressive n,% calculated (using progressive denominator)

# # Remove patients with data_start > data_end errors (as we can't accurately calculate their follow-up time)
# total_prog <- n_distinct(cohort_England_2015_2020$patid)
# data_start_data_end_errors_prog <- sum(cohort_England_2015_2020$data_start > cohort_England_2015_2020$data_end) # 0
# data_start_data_end_errors_percent_prog <- data_start_data_end_errors_prog/total_prog # 0%
# cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
#   filter(data_start <= data_end)
# summary(cohort_England_2015_2020$pyears)
# 
# # Remove patients with indeterminate gender
# total_prog2 <- n_distinct(cohort_England_2015_2020$patid)
# indeterminate_gender_prog <- sum(cohort_England_2015_2020$gender == 'Indeterminate') # 
# indeterminate_gender_percent_prog <- indeterminate_gender_prog/total_prog2 # %
# cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
#   filter(gender != "Indeterminate")
# cohort_England_2015_2020$gender <- droplevels(cohort_England_2015_2020$gender)


# Exclude possible migrants
total_prog3 <- n_distinct(cohort_England_2015_2020$patid)
possible_migrant_prog <- sum(cohort_England_2015_2020$migcertainty == 'Possible') 
possible_migrant_percent_prog <- possible_migrant_prog/total_prog3 
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  filter(migcertainty != "Possible")

cohort_England_2015_2020$migcertainty <- droplevels(cohort_England_2015_2020$migcertainty)
levels(cohort_England_2015_2020$migcertainty)
glimpse(cohort_England_2015_2020) 

# Save cohort file 
save(cohort_England_2015_2020, file = "filepath/cohort_England_2015_2020.Rdata")


####################### CREATE ANNUAL RECORDS ######################### --------------------


# 12_Create annual cohort records  ----

# Create variable with repeated sequence of years
year_variable <- as.data.frame(rep(seq(2015,2020,1), times = nrow(cohort_England_2015_2020)))
colnames(year_variable) <- 'studyyear'

# Replicate each patient 6 times and label by year 2015-2020
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  slice(rep(1:n(), each=6))

# Add the repeated year variable
cohort_England_2015_2020$studyyear <- year_variable$studyyear

# Filter to include patients only in the years that they were active in CPRD 
cohort_England_2015_2020_annual_records <- filter(cohort_England_2015_2020, studyyear >= year(cohort_entry) & studyyear <= year(cohort_exit))

# Calculate yearly follow-up time for each patient 
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  mutate(annual_cohort_entry = as_date(ifelse(year(cohort_entry) == studyyear, cohort_entry,as_date(paste0(studyyear,"-01-01"))))) %>%
  mutate(annual_cohort_exit = as_date(ifelse(year(cohort_exit) > studyyear,as_date(paste0(studyyear,"-12-31")) ,cohort_exit))) %>%
  mutate(annual_cohort_exit = replace(annual_cohort_exit, annual_cohort_exit > "2020-12-26", as_date("2020-12-26")))%>%
  mutate(pdays = (annual_cohort_exit - annual_cohort_entry) + 1) %>%
  mutate(pyears = pdays/365)  
cohort_England_2015_2020_annual_records$pyears <- as.numeric(cohort_England_2015_2020_annual_records$pyears)
cohort_England_2015_2020_annual_records$pdays <- as.numeric(cohort_England_2015_2020_annual_records$pdays)

glimpse(cohort_England_2015_2020_annual_records)
summary(cohort_England_2015_2020_annual_records$pdays)
summary(cohort_England_2015_2020_annual_records$pyears)

# 13_Create age variables for annual records and cohort file ----

# Create age variable (5 year increments)
cohort_England_2015_2020_annual_records$studyyear <- as.integer(cohort_England_2015_2020_annual_records$studyyear)
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  mutate(age = studyyear - yob)

## Exclusions made and progressive n,% calculated (using progressive denominator)

total_prog4 <- n_distinct(cohort_England_2015_2020_annual_records$patid)
age_prog <- cohort_England_2015_2020_annual_records %>% filter(age >=100)
age_prog <- n_distinct(age_prog$patid) 
age_prog_percent <- age_prog/total_prog4 
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>% # 91833931 - 91782923 = 51008 obs excluded 
  filter(age < 100)

# Re-load cohort_England_2015_2020
load(file = "filepath/cohort_England_2015_2020.Rdata")

# Filter cohort_England_2015_2020 by those in annual cohort file and save for cohort characteristics 
cohort_England_2015_2020_allages <- n_distinct(cohort_England_2015_2020$patid)
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  filter(patid %in% cohort_England_2015_2020_annual_records$patid)
cohort_England_2015_2020_under100 <- n_distinct(cohort_England_2015_2020$patid)
excluded_100_and_over <- cohort_England_2015_2020_allages - cohort_England_2015_2020_under100 # 225 distinct patients fully excluded
excluded_100_and_over_percent <- excluded_100_and_over/total_prog4 # 0.0374% fully excluded
n_distinct(cohort_England_2015_2020_annual_records$patid) == nrow(cohort_England_2015_2020) # T

save(cohort_England_2015_2020, file = "filepath/cohort_England_2015_2020.Rdata") # 601033 patients in final cohort
rm(cohort_England_2015_2020)

# Create studyyear_agecat variable
cohort_England_2015_2020_annual_records <-cohort_England_2015_2020_annual_records %>%
  mutate(studyyear_agecat = age)

cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, studyyear_agecat <6, 1)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 6, 10), 2)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 11, 15), 3)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 16, 19), 4)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 20, 24), 5)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 25, 29), 6)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 30, 34), 7)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 35, 39), 8)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 40, 44), 9)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 45, 49), 10)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 50, 54), 11)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 55, 59), 12)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 60, 64), 13)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 65, 69), 14)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 70, 74), 15)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 75, 79), 16)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 80, 84), 17)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 85, 89), 18)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 90, 94), 19)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 95, 99), 20)) 

cohort_England_2015_2020_annual_records$studyyear_agecat <- factor(cohort_England_2015_2020_annual_records$studyyear_agecat, 
                                                               levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                                                               labels = c("0-5 years", "6-10 years", "11-15 years",
                                                                          "16-19 years", "20-24 years", "25-29 years",
                                                                          "30-34 years", "35-39 years", "40-44 years",
                                                                          "45-49 years", "50-54 years", "55-59 years",
                                                                          "60-64 years", "65-69 years", "70-74 years",
                                                                          "75-79 years", "80-84 years", "85-89 years",
                                                                          "90-94 years", "95-99 years"))


# Create age_subcohort variable

cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  mutate(age_subcohort = age)

cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  mutate(age_subcohort=replace(age_subcohort, age_subcohort <16, 1)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 16, 24), 2)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 25, 34), 3)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 35, 49), 4)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 50, 64), 5)) %>%
  mutate(age_subcohort=replace(age_subcohort, age_subcohort >= 65, 6)) 

cohort_England_2015_2020_annual_records$age_subcohort <- factor(cohort_England_2015_2020_annual_records$age_subcohort,levels = c(1,2,3,4,5,6),
                                                            labels = c("0-15 years", "16-24 years", "25-34 years",
                                                                       "35-49 years","50-64 years", ">=65 years"))

glimpse(cohort_England_2015_2020_annual_records)
any(is.na(cohort_England_2015_2020_annual_records$studyyear_agecat)) 
any(is.na(cohort_England_2015_2020_annual_records$age_subcohort)) 

save(cohort_England_2015_2020_annual_records, file = "filepath/cohort_England_2015_2020_annual_records.Rdata")


############################ EXACT COHORT MATCHING ########################## ----------------------


# 14_Annual analysis ----
## 14a_Sensitivity analysis cohort 1: Matched on age at data_start, year at data_start and prac_region) ------------------------------------------------------

load(file = "filepath/cohort_England_2015_2020.Rdata")

# Derive age at data_start
cohort_England_2015_2020_test_ds <- cohort_England_2015_2020 %>%
  mutate(dob = lubridate::ymd(yob,truncated=2L))

calc_age <- function(birthDate, refDate = Sys.Date()) {
  require(lubridate)
  period <- as.period(interval(birthDate, refDate),unit = "year")
  period$year
}

cohort_England_2015_2020_test_ds <- cohort_England_2015_2020_test_ds %>%
  mutate(age_data_start= calc_age(dob, data_start))

# Derive year of data_start
cohort_England_2015_2020_test_ds <- cohort_England_2015_2020_test_ds %>%
  mutate(year_data_start= lubridate::year(ymd(data_start)))

# turn migrant status into binary integer
levels(cohort_England_2015_2020_test_ds$migrant_status) # level 1 = Non-migrant, level 2 = Migrant
cohort_England_2015_2020_test_ds <- cohort_England_2015_2020_test_ds %>%
  mutate(migrant_status_binary = as.numeric(migrant_status)) %>%
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 1, 0)) %>% # Non-migrant
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 2, 1)) # Migrant

# select variables needed
cohort_England_2015_2020_test_ds_2 <- select(cohort_England_2015_2020_test_ds, c(patid, migrant_status_binary,year_data_start, age_data_start, prac_region))

# Turn df into data table

setDT(cohort_England_2015_2020_test_ds_2)

# stratifies dataset and then selects random observations within strata
# data = dataset containing:
# - treatment/exposure variable 'mvar' (a string specifying variable name).
# - matching variable 'mvar' (a string specifying variable name). If you want to match on multiple variables, concatenate them first.
# other inputs are:
# - ratio of cases:controls (an integer > 0). You can also set ratio between 0 and 1, e.g. to 0.5 to get 2 cases/control, but the results aren't perfect
# - seed for fixing random selection of cases/controls (an integer; default NULL means no seed). Choice of seed is arbitrary. Seed will be reset globally at the end of the function (i.e. the function will overwrite a previous seed)
# returns data.table of matched observations, with additional variable 'id' for use in paired/grouped analyses
# the speed of the function mostly depends on the number of strata (i.e. levels of the matching variable)
# there is no problem converting the results to a standard data frame with data.frame(matched_data) if you prefer working in this way

smatch <- function (data, treat, mvar, ratio = 1, seed = NULL) {
  targ <- data[, .(case = sum(get(treat)), control = sum(!get(treat))), mvar]
  targ[, cst := floor(pmin(control / ratio, case))]
  targ[, cnt := cst * ratio]
  targ <- targ[cst > 0]
  setnames(targ, mvar, 'mvar')
  l2 <- cumsum(targ$cst)
  ids <- mapply(':', c(0, l2[-nrow(targ)]), l2-1)
  names(ids) <- targ$mvar
  case <- NULL
  control <- NULL
  set.seed(seed)
  on.exit(set.seed(NULL))
  for(i in targ$mvar) {
    case[[i]] <- data[get(treat) == T & get(mvar) == i][sample(.N, targ$cst[targ$mvar == i])]
    case[[i]][, id := ids[[i]]]
    control[[i]] <- data[get(treat) == F & get(mvar) == i][sample(.N, targ$cnt[targ$mvar == i])]
    control[[i]][, id := rep(ids[[i]], each = ratio)]
  }
  rbindlist(c(case, control))
}

# concatenate multiple variables
cohort_England_2015_2020_test_ds_2[, age_year_region := paste0(age_data_start, '-', year_data_start, '-', prac_region)]

# create matched dataset

# matched_data <- smatch(cohort_England_2015_2020_2, 'migrant_status_binary', 'age_year_region') # gives a different matched set each time
matched_data <- smatch(cohort_England_2015_2020_test_ds_2, 'migrant_status_binary', 'age_year_region', seed = 5) # gives the same matched set each time (the number 5 is arbitrary)

# check balance
dcast(matched_data, age_year_region ~ migrant_status_binary, value.var = 'age_data_start', fun.aggregate = length)

# Turn back into dataframe
matched_cohort_England_2015_2020_test_ds <- as.data.frame(matched_data)

# Rejoin matched cohort to rest of dataset
matched_cohort_England_2015_2020_test_ds <- select(matched_cohort_England_2015_2020_test_ds, -c(prac_region, year_data_start, age_data_start, migrant_status_binary))
matched_cohort_England_2015_2020_test_ds <- left_join(matched_cohort_England_2015_2020_test_ds, cohort_England_2015_2020_test_ds, by = c("patid"="patid"))

matched_cohort_England_2015_2020_test_ds <- select(matched_cohort_England_2015_2020_test_ds, -c(migrant_status_binary, age_year_region))
glimpse(matched_cohort_England_2015_2020_test_ds)
sum(matched_cohort_England_2015_2020_test_ds$migrant_status == "Non-migrant") == sum(matched_cohort_England_2015_2020_test_ds$migrant_status == "Migrant") 
save(matched_cohort_England_2015_2020_test_ds, file = "filepath/matched_cohort_England_2015_2020_test_ds.Rdata")

# Create annual records file
load(file = "filepath/cohort_England_2015_2020_annual_records.Rdata")
matched_cohort_England_2015_2020_annual_records_test_ds <- cohort_England_2015_2020_annual_records %>%
  filter(patid %in% matched_cohort_England_2015_2020_test_ds$patid)
n_distinct(matched_cohort_England_2015_2020_annual_records_test_ds$patid) == nrow(matched_cohort_England_2015_2020_test_ds) 

# Check no one 100 years old or over
sum(matched_cohort_England_2015_2020_annual_records_test_ds$age >= 100)

# Save annual cohort file 
save(matched_cohort_England_2015_2020_annual_records_test_ds, file = "filepath/matched_cohort_England_2015_2020_annual_records_test_ds.Rdata")

## 14b_Sensitivity analysis cohort 2: Matched on person-years and prac_region------------------------------------------------------

load(file = "filepath/cohort_England_2015_2020.Rdata")

# turn migrant status into binary integer
levels(cohort_England_2015_2020$migrant_status) # level 1 = Non-migrant, level 2 = Migrant
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(migrant_status_binary = as.numeric(migrant_status)) %>%
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 1, 0)) %>% # Non-migrant
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 2, 1)) # Migrant

# select variables needed
cohort_England_2015_2020_2 <- select(cohort_England_2015_2020, c(patid, migrant_status_binary,prac_region, total_pyears))

# Turn df into data table

setDT(cohort_England_2015_2020_2)

# concatenate multiple variables
cohort_England_2015_2020_2[, region_followup := paste0(prac_region, '-', total_pyears)]

# create matched dataset

matched_data <- smatch(cohort_England_2015_2020_2, 'migrant_status_binary', 'region_followup', seed = 5) # gives the same matched set each time (the number 5 is arbitrary)

# check balance
dcast(matched_data, region_followup ~ migrant_status_binary, value.var = 'total_pyears', fun.aggregate = length)

# Turn back into dataframe
matched_cohort_England_2015_2020_fu <- as.data.frame(matched_data)

# Rejoin matched cohort to rest of dataset
matched_cohort_England_2015_2020_fu <- select(matched_cohort_England_2015_2020_fu, -c(prac_region, total_pyears, migrant_status_binary))
matched_cohort_England_2015_2020_fu <- left_join(matched_cohort_England_2015_2020_fu, cohort_England_2015_2020, by = c("patid"="patid"))

matched_cohort_England_2015_2020_fu <- select(matched_cohort_England_2015_2020_fu, -c(migrant_status_binary, region_followup))
glimpse(matched_cohort_England_2015_2020_fu)
sum(matched_cohort_England_2015_2020_fu$migrant_status == "Non-migrant") == sum(matched_cohort_England_2015_2020_fu$migrant_status == "Migrant")
save(matched_cohort_England_2015_2020_fu, file = "filepath/matched_cohort_England_2015_2020_fu.Rdata")

# Create annual records file
load(file = "primary_care_cleaned_files/01_Cohort_Rdata/cohort_England_2015_2020_annual_records.Rdata")
matched_cohort_England_2015_2020_annual_records_fu <- cohort_England_2015_2020_annual_records %>%
  filter(patid %in% matched_cohort_England_2015_2020_fu$patid)
n_distinct(matched_cohort_England_2015_2020_annual_records_fu$patid) == nrow(matched_cohort_England_2015_2020_fu)

# Check no one 100 years old or over
sum(matched_cohort_England_2015_2020_annual_records_fu$age >= 100)

# Save annual cohort file 
save(matched_cohort_England_2015_2020_annual_records_fu, file = "filepath/matched_cohort_England_2015_2020_annual_records_fu.Rdata")

# 15_ITS analysis -------
## 15a_Main ITS analysis cohort matched on age at study start, prac_region, gender and imd ----

### Create cohort----

# Load cohort and annual records file
load(file = "filepath/cohort_England_2015_2020_annual_records.Rdata")
load(file = "filepath/cohort_England_2015_2020.Rdata")

# Derive age at study start
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(dob = lubridate::ymd(yob,truncated=2L))

cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(age_study_start= calc_age(dob, as.Date('2015-01-01')))

# turn migrant status into binary integer
levels(cohort_England_2015_2020$migrant_status) # level 1 = Non-migrant, level 2 = Migrant
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(migrant_status_binary = as.numeric(migrant_status)) %>%
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 1, 0)) %>% # Non-migrant
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 2, 1)) # Migrant

# select variables needed
cohort_England_2015_2020_2 <- select(cohort_England_2015_2020, c(patid, migrant_status_binary, age_study_start, prac_region, gender, imd))

# Turn df into data table

setDT(cohort_England_2015_2020_2)

# concatenate multiple variables
cohort_England_2015_2020_2[, age_region_gender_imd := paste0(age_study_start, '-', prac_region, '-', gender, '-', imd)]

# create matched dataset

matched_data <- smatch(cohort_England_2015_2020_2, 'migrant_status_binary', 'age_region_gender_imd', seed = 5) # gives the same matched set each time (the number 5 is arbitrary)

# check balance
dcast(matched_data, age_region_gender_imd ~ migrant_status_binary, value.var = 'age_study_start', fun.aggregate = length)

# Turn back into dataframe
matched_cohort_England_2015_2020 <- as.data.frame(matched_data)

# Rejoin matched cohort to rest of dataset
matched_cohort_England_2015_2020 <- select(matched_cohort_England_2015_2020, -c(prac_region, age_study_start ,migrant_status_binary))
matched_its_cohort <- left_join(matched_cohort_England_2015_2020, cohort_England_2015_2020, by = c("patid"="patid", 'gender'= 'gender',
                                                                                                  'imd' = 'imd'))

matched_its_cohort <- select(matched_its_cohort, -c(migrant_status_binary, age_region_gender_imd))
glimpse(matched_its_cohort) 
sum(matched_its_cohort$migrant_status == "Non-migrant") == sum(matched_its_cohort$migrant_status == "Migrant") 
save(matched_its_cohort, file = "filepath/matched_its_cohort.Rdata")

### Create weekly records  ---------------------------

# set up week to date mapping files
## create data_week_mapping_1 - rename eventdate into data_start, week into data_start_week
## create data_week_mapping_2 - rename eventdate into data_end, week into data_end_week
## data_week_mapping create third column (varname = first_date_week) which is = to the first date of that week

date_to_week_conversion_2015_2020 <- data.frame(date = seq(as.Date('2015-01-01'), as.Date('2020-12-26'), by='days'))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(week = cut.Date(date, breaks = '1 week', labels = FALSE, start.on.monday = FALSE)) %>%
  arrange(date)

n_distinct(date_to_week_conversion_2015_2020$date) # 2187 days, 313 weeks

save(date_to_week_conversion_2015_2020, file = "filepath/date_to_week_conversion_2015_2020.Rdata") 

data_week_mapping_1 <- date_to_week_conversion_2015_2020 %>%
  rename(data_start = date) %>%
  rename(data_start_week = week)
data_week_mapping_2 <- date_to_week_conversion_2015_2020 %>%
  rename(data_end = date) %>%
  rename(data_end_week = week)
data_week_mapping_3 <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  summarise(first_date_week = first(date))

load(file= "filepath/matched_its_cohort.Rdata")

matched_its_weekly_records <- matched_its_cohort %>%
  select(patid, pracid, data_start, data_end, 
         migrant_status, migcertainty, ethnicat, ethnicat6, prac_region, imd, gender, yob)

rm(matched_its_cohort)

n_distinct(matched_its_weekly_records$patid) 
n_distinct(matched_its_weekly_records$patid[matched_its_weekly_records$migrant_status=="Migrant"]) 

# Create week variable 
studyweek <- as.data.frame(rep(1:313,times = nrow(matched_its_weekly_records))) # 312 weeks and 2 days between 1 Jan 2015 and 26 Dec 2020
colnames(studyweek) <- 'studyweek'

# Replicate patients for the number of weeks in the analysis 
matched_its_weekly_records <- matched_its_weekly_records %>%
  slice(rep(1:n(), each=313))
matched_its_weekly_records$studyweek <- studyweek$studyweek

# Select only weeks where patients were active in CPRD
## left join matched_its_weekly_records to data_week_mapping_1 by data_start and same for 2 by data_end
## replace data_start_week NAs with -1

matched_its_weekly_records <- left_join(matched_its_weekly_records,
                                                    data_week_mapping_1, by=c("data_start" = "data_start")) # NAs represent individuals who had a data_start before the study period
matched_its_weekly_records <- left_join(matched_its_weekly_records, 
                                                    data_week_mapping_2, by=c("data_end" = "data_end"))
sum(is.na(matched_its_weekly_records$data_start_week)) 
sum(is.na(matched_its_weekly_records$data_end_week)) 
matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(data_start_week = replace(data_start_week, is.na(data_start_week), -1))

n_distinct(matched_its_weekly_records$patid) 

matched_its_weekly_records <- filter(matched_its_weekly_records, 
                                                 studyweek >= data_start_week &
                                                   studyweek <= data_end_week)

n_distinct(matched_its_weekly_records$patid) 

# Derive follow-up time for each patient/week
## deselect eventdate or whatever it's called
## join by studyweek = week 

matched_its_weekly_records <- left_join(matched_its_weekly_records, data_week_mapping_3, by=c("studyweek" = "week"))

matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(pdays = ifelse((studyweek != 1 & data_start_week < studyweek & data_end_week > studyweek), 7,
                        ifelse((studyweek == 1 & data_start_week < studyweek & data_end_week > studyweek), 3,
                               ifelse((data_start_week == studyweek & data_end_week == studyweek), data_end-data_start+1,
                                      ifelse((studyweek !=1 & data_start_week == studyweek & data_end_week > studyweek), 7-(data_start-first_date_week),
                                             ifelse((studyweek ==1 & data_start_week == studyweek & data_end_week > studyweek), 3-(data_start-first_date_week),
                                                    ifelse((data_start_week < studyweek & data_end_week == studyweek), data_end-first_date_week+1, NA)))))))


matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(pyears = pdays/365,
         pweeks = pdays/7)

glimpse(matched_its_weekly_records)
summary(matched_its_weekly_records$pdays)
summary(matched_its_weekly_records$pyears)
summary(matched_its_weekly_records$pweeks)

# Create studymonth variable
matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(studymonth = month(first_date_week))

# Create studyyear variable 
matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(studyyear = year(first_date_week))


matched_its_weekly_records <- select(matched_its_weekly_records,-c(data_start_week, data_end_week))

# Create age variable (5 year increments)
matched_its_weekly_records$studyyear <- as.integer(matched_its_weekly_records$studyyear)
matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(age = studyyear - yob)
matched_its_weekly_records <- matched_its_weekly_records %>% 
  filter(age < 100)

# Create studyyear_agecat variable
matched_its_weekly_records <-matched_its_weekly_records %>%
  mutate(studyyear_agecat = age)

matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, studyyear_agecat <6, 1)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 6, 10), 2)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 11, 15), 3)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 16, 19), 4)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 20, 24), 5)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 25, 29), 6)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 30, 34), 7)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 35, 39), 8)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 40, 44), 9)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 45, 49), 10)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 50, 54), 11)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 55, 59), 12)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 60, 64), 13)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 65, 69), 14)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 70, 74), 15)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 75, 79), 16)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 80, 84), 17)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 85, 89), 18)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 90, 94), 19)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 95, 99), 20)) 

matched_its_weekly_records$studyyear_agecat <- factor(matched_its_weekly_records$studyyear_agecat, 
                                                         levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                                                         labels = c("0-5 years", "6-10 years", "11-15 years",
                                                                    "16-19 years", "20-24 years", "25-29 years",
                                                                    "30-34 years", "35-39 years", "40-44 years",
                                                                    "45-49 years", "50-54 years", "55-59 years",
                                                                    "60-64 years", "65-69 years", "70-74 years",
                                                                    "75-79 years", "80-84 years", "85-89 years",
                                                                    "90-94 years", "95-99 years"))


# Create age_subcohort variable

matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(age_subcohort = age)

matched_its_weekly_records <- matched_its_weekly_records %>%
  mutate(age_subcohort=replace(age_subcohort, age_subcohort <16, 1)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 16, 24), 2)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 25, 34), 3)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 35, 49), 4)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 50, 64), 5)) %>%
  mutate(age_subcohort=replace(age_subcohort, age_subcohort >= 65, 6)) 

matched_its_weekly_records$age_subcohort <- factor(matched_its_weekly_records$age_subcohort,levels = c(1,2,3,4,5,6),
                                                      labels = c("0-15 years", "16-24 years", "25-34 years",
                                                                 "35-49 years","50-64 years", ">=65 years"))

glimpse(matched_its_weekly_records)
any(is.na(matched_its_weekly_records$studyyear_agecat)) 
any(is.na(matched_its_weekly_records$age_subcohort))
n_distinct(matched_its_weekly_records$patid)

save(matched_its_weekly_records, file = "filepath/01_Cohort_Rdata/matched_its_weekly_records.Rdata")

rm(list=ls())
gc()



## 15b_Sensitivity ITS analysis cohort matched on age at data_start, year at data_start, gender, prac_region and imd -----
### Create cohort  ------------------------------------------------------

# Load cohort and annual records file
load(file = "filepath/cohort_England_2015_2020_annual_records.Rdata")
load(file = "filepath/cohort_England_2015_2020.Rdata")

# Derive age at study start
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(dob = lubridate::ymd(yob,truncated=2L))

calc_age <- function(birthDate, refDate = Sys.Date()) {
  require(lubridate)
  period <- as.period(interval(birthDate, refDate),unit = "year")
  period$year
}

cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(age_data_start= calc_age(dob, as.Date(data_start)))

# Derive year of cohort entry
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(year_data_start = lubridate::year(ymd(data_start)))

# turn migrant status into binary integer
levels(cohort_England_2015_2020$migrant_status) # level 1 = Non-migrant, level 2 = Migrant
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(migrant_status_binary = as.numeric(migrant_status)) %>%
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 1, 0)) %>% # Non-migrant
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 2, 1)) # Migrant

# select variables needed
cohort_England_2015_2020_2 <- dplyr::select(cohort_England_2015_2020, c(patid, migrant_status_binary, age_data_start, year_data_start, prac_region, gender, imd))

# Turn df into data table

setDT(cohort_England_2015_2020_2)

# concatenate multiple variables
cohort_England_2015_2020_2[, age_year_region_gender_imd := paste0(age_data_start, '-', year_data_start, '-', prac_region, '-', gender, '-', imd)]

# create matched dataset

matched_data <- smatch(cohort_England_2015_2020_2, 'migrant_status_binary', 'age_year_region_gender_imd', seed = 5) # gives the same matched set each time (the number 5 is arbitrary)

# check balance
dcast(matched_data, age_year_region_gender_imd ~ migrant_status_binary, value.var = 'age_data_start', fun.aggregate = length)

# Turn back into dataframe
matched_cohort_England_2015_2020 <- as.data.frame(matched_data)

# Rejoin matched cohort to rest of dataset
matched_cohort_England_2015_2020 <- dplyr::select(matched_cohort_England_2015_2020, -c(prac_region, age_data_start, year_data_start ,migrant_status_binary))
matched_its_cohort_ds <- left_join(matched_cohort_England_2015_2020, cohort_England_2015_2020, by = c("patid"="patid", 'gender'= 'gender',
                                                                                                      'imd' = 'imd'))

matched_its_cohort_ds <- dplyr::select(matched_its_cohort_ds, -c(migrant_status_binary, age_year_region_gender_imd))
glimpse(matched_its_cohort_ds) 
sum(matched_its_cohort_ds$migrant_status == "Non-migrant") == sum(matched_its_cohort_ds$migrant_status == "Migrant") 
save(matched_its_cohort_ds, file = "filepath/matched_its_cohort_ds.Rdata")

### Create weekly records ------

# set up week to date mapping files
## create data_week_mapping_1 - rename eventdate into data_start, week into data_start_week
## create data_week_mapping_2 - rename eventdate into data_end, week into data_end_week
## data_week_mapping create third column (varname = first_date_week) which is = to the first date of that week

load(file = "filepath/date_to_week_conversion_2015_2020.Rdata") 

data_week_mapping_1 <- date_to_week_conversion_2015_2020 %>%
  rename(data_start = date) %>%
  rename(data_start_week = week)
data_week_mapping_2 <- date_to_week_conversion_2015_2020 %>%
  rename(data_end = date) %>%
  rename(data_end_week = week)
data_week_mapping_3 <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  summarise(first_date_week = first(date))


load(file= "filepath/matched_its_cohort_ds.Rdata")

matched_its_weekly_records_ds <- matched_its_cohort_ds %>%
  select(patid, pracid, data_start, data_end, 
         migrant_status, migcertainty, ethnicat6, prac_region,
         imd, gender, yob)

rm(matched_its_cohort_ds)

n_distinct(matched_its_weekly_records_ds$patid) 
n_distinct(matched_its_weekly_records_ds$patid[matched_its_weekly_records_ds$migrant_status=="Migrant"]) 

# Create week variable (or alternatively do years separately and then combine [as 2020 is a leap yr])
studyweek <- as.data.frame(rep(1:313,times = nrow(matched_its_weekly_records_ds))) # 312 weeks and 2 days between 1 Jan 2015 and 26 Dec 2020
colnames(studyweek) <- 'studyweek'

# Replicate patients for the number of months in the analysis 
matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  slice(rep(1:n(), each=313))
matched_its_weekly_records_ds$studyweek <- studyweek$studyweek
rm(studyweek)

# Select only months where patients were active in CPRD
## left join matched_its_weekly_records to data_week_mapping_1 by data_start and same for 2 by data_end
## replace data_start_week NAs with -1

matched_its_weekly_records_ds <- left_join(matched_its_weekly_records_ds,
                                           data_week_mapping_1, by=c("data_start" = "data_start"))
matched_its_weekly_records_ds <- left_join(matched_its_weekly_records_ds, 
                                           data_week_mapping_2, by=c("data_end" = "data_end"))
sum(is.na(matched_its_weekly_records_ds$data_start_week)) 
sum(is.na(matched_its_weekly_records_ds$data_end_week))
matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(data_start_week = replace(data_start_week, is.na(data_start_week), -1))

n_distinct(matched_its_weekly_records_ds$patid) 

matched_its_weekly_records_ds <- filter(matched_its_weekly_records_ds, 
                                        studyweek >= data_start_week &
                                          studyweek <= data_end_week)

n_distinct(matched_its_weekly_records_ds$patid) 

# Derive follow-up time for each patient/week
## deselect eventdate or whatever it's called
## join by studyweek = week 

matched_its_weekly_records_ds <- left_join(matched_its_weekly_records_ds, data_week_mapping_3, by=c("studyweek" = "week"))

matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(pdays = ifelse((studyweek != 1 & data_start_week < studyweek & data_end_week > studyweek), 7,
                        ifelse((studyweek == 1 & data_start_week < studyweek & data_end_week > studyweek), 3,
                               ifelse((data_start_week == studyweek & data_end_week == studyweek), data_end-data_start+1,
                                      ifelse((studyweek !=1 & data_start_week == studyweek & data_end_week > studyweek), 7-(data_start-first_date_week),
                                             ifelse((studyweek ==1 & data_start_week == studyweek & data_end_week > studyweek), 3-(data_start-first_date_week),
                                                    ifelse((data_start_week < studyweek & data_end_week == studyweek), data_end-first_date_week+1, NA)))))))


matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(pyears = pdays/365,
         pweeks = pdays/7)

# Create studymonth variable
matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(studymonth = month(first_date_week))

# Create studyyear variable 
matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(studyyear = year(first_date_week))

glimpse(matched_its_weekly_records_ds)
summary(matched_its_weekly_records_ds$pdays)
summary(matched_its_weekly_records_ds$pyears)
summary(matched_its_weekly_records_ds$pweeks)

matched_its_weekly_records_ds <- select(matched_its_weekly_records_ds,-c(data_start_week, data_end_week))

# Create age variable (5 year increments)
matched_its_weekly_records_ds$studyyear <- as.integer(matched_its_weekly_records_ds$studyyear)
matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(age = studyyear - yob)
matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>% 
  filter(age < 100)

# Create studyyear_agecat variable
matched_its_weekly_records_ds <-matched_its_weekly_records_ds %>%
  mutate(studyyear_agecat = age)

matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, studyyear_agecat <6, 1)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 6, 10), 2)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 11, 15), 3)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 16, 19), 4)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 20, 24), 5)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 25, 29), 6)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 30, 34), 7)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 35, 39), 8)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 40, 44), 9)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 45, 49), 10)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 50, 54), 11)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 55, 59), 12)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 60, 64), 13)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 65, 69), 14)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 70, 74), 15)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 75, 79), 16)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 80, 84), 17)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 85, 89), 18)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 90, 94), 19)) %>%
  mutate(studyyear_agecat=replace(studyyear_agecat, between (studyyear_agecat, 95, 99), 20)) 

matched_its_weekly_records_ds$studyyear_agecat <- factor(matched_its_weekly_records_ds$studyyear_agecat, 
                                                         levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                                                         labels = c("0-5 years", "6-10 years", "11-15 years",
                                                                    "16-19 years", "20-24 years", "25-29 years",
                                                                    "30-34 years", "35-39 years", "40-44 years",
                                                                    "45-49 years", "50-54 years", "55-59 years",
                                                                    "60-64 years", "65-69 years", "70-74 years",
                                                                    "75-79 years", "80-84 years", "85-89 years",
                                                                    "90-94 years", "95-99 years"))


# Create age_subcohort variable

matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(age_subcohort = age)

matched_its_weekly_records_ds <- matched_its_weekly_records_ds %>%
  mutate(age_subcohort=replace(age_subcohort, age_subcohort <16, 1)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 16, 24), 2)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 25, 34), 3)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 35, 49), 4)) %>%
  mutate(age_subcohort=replace(age_subcohort, between (age_subcohort, 50, 64), 5)) %>%
  mutate(age_subcohort=replace(age_subcohort, age_subcohort >= 65, 6)) 

matched_its_weekly_records_ds$age_subcohort <- factor(matched_its_weekly_records_ds$age_subcohort,levels = c(1,2,3,4,5,6),
                                                      labels = c("0-15 years", "16-24 years", "25-34 years",
                                                                 "35-49 years","50-64 years", ">=65 years"))

glimpse(matched_its_weekly_records_ds)
any(is.na(matched_its_weekly_records_ds$studyyear_agecat)) 
any(is.na(matched_its_weekly_records_ds$age_subcohort))
n_distinct(matched_its_weekly_records_ds$patid)

save(matched_its_weekly_records_ds, file = "filepath/matched_its_weekly_records_ds.Rdata")

rm(list=ls())
gc()

# 16_Cohorts not included in published paper----

## 16a_Annual analysis sensitivity analysis cohort 3: Matched on age at cohort entry, year at cohort entry and prac_region ------------------------------------------------------

load(file = "filepath/cohort_England_2015_2020.Rdata")

# Derive age at cohort entry
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(dob = lubridate::ymd(yob,truncated=2L))

calc_age <- function(birthDate, refDate = Sys.Date()) {
  require(lubridate)
  period <- as.period(interval(birthDate, refDate),unit = "year")
  period$year
}

cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(age_cohort_entry= calc_age(dob, cohort_entry))

# Derive year of cohort entry
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(year_cohort_entry = lubridate::year(ymd(cohort_entry)))

# turn migrant status into binary integer
levels(cohort_England_2015_2020$migrant_status) # level 1 = Non-migrant, level 2 = Migrant
cohort_England_2015_2020 <- cohort_England_2015_2020 %>%
  mutate(migrant_status_binary = as.numeric(migrant_status)) %>%
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 1, 0)) %>% # Non-migrant
  mutate(migrant_status_binary = replace(migrant_status_binary, migrant_status_binary == 2, 1)) # Migrant

# select variables needed
cohort_England_2015_2020_2 <- select(cohort_England_2015_2020, c(patid, migrant_status_binary,year_cohort_entry, age_cohort_entry, prac_region))

# function requires data.table

library(data.table)

# Turn df into data table

setDT(cohort_England_2015_2020_2)

# stratifies dataset and then selects random observations within strata
# data = dataset containing:
# - treatment/exposure variable 'mvar' (a string specifying variable name).
# - matching variable 'mvar' (a string specifying variable name). If you want to match on multiple variables, concatenate them first.
# other inputs are:
# - ratio of cases:controls (an integer > 0). You can also set ratio between 0 and 1, e.g. to 0.5 to get 2 cases/control, but the results aren't perfect
# - seed for fixing random selection of cases/controls (an integer; default NULL means no seed). Choice of seed is arbitrary. Seed will be reset globally at the end of the function (i.e. the function will overwrite a previous seed)
# returns data.table of matched observations, with additional variable 'id' for use in paired/grouped analyses
# the speed of the function mostly depends on the number of strata (i.e. levels of the matching variable)
# there is no problem converting the results to a standard data frame with data.frame(matched_data) if you prefer working in this way

smatch <- function (data, treat, mvar, ratio = 1, seed = NULL) {
  targ <- data[, .(case = sum(get(treat)), control = sum(!get(treat))), mvar]
  targ[, cst := floor(pmin(control / ratio, case))]
  targ[, cnt := cst * ratio]
  targ <- targ[cst > 0]
  setnames(targ, mvar, 'mvar')
  l2 <- cumsum(targ$cst)
  ids <- mapply(':', c(0, l2[-nrow(targ)]), l2-1)
  names(ids) <- targ$mvar
  case <- NULL
  control <- NULL
  set.seed(seed)
  on.exit(set.seed(NULL))
  for(i in targ$mvar) {
    case[[i]] <- data[get(treat) == T & get(mvar) == i][sample(.N, targ$cst[targ$mvar == i])]
    case[[i]][, id := ids[[i]]]
    control[[i]] <- data[get(treat) == F & get(mvar) == i][sample(.N, targ$cnt[targ$mvar == i])]
    control[[i]][, id := rep(ids[[i]], each = ratio)]
  }
  rbindlist(c(case, control))
}

# concatenate multiple variables
cohort_England_2015_2020_2[, age_year_region := paste0(age_cohort_entry, '-', year_cohort_entry, '-', prac_region)]

# create matched dataset

# matched_data <- smatch(cohort_England_2015_2020_2, 'migrant_status_binary', 'age_year_region') # gives a different matched set each time
matched_data <- smatch(cohort_England_2015_2020_2, 'migrant_status_binary', 'age_year_region', seed = 5) # gives the same matched set each time (the number 5 is arbitrary)

# check balance
dcast(matched_data, age_year_region ~ migrant_status_binary, value.var = 'age_cohort_entry', fun.aggregate = length)

# Turn back into dataframe
matched_cohort_England_2015_2020 <- as.data.frame(matched_data)

# Rejoin matched cohort to rest of dataset
matched_cohort_England_2015_2020 <- select(matched_cohort_England_2015_2020, -c(prac_region, year_cohort_entry, age_cohort_entry,migrant_status_binary))
matched_cohort_England_2015_2020 <- left_join(matched_cohort_England_2015_2020, cohort_England_2015_2020, by = c("patid"="patid"))

matched_cohort_England_2015_2020 <- select(matched_cohort_England_2015_2020, -c(migrant_status_binary, age_year_region))
glimpse(matched_cohort_England_2015_2020)
sum(matched_cohort_England_2015_2020$migrant_status == "Non-migrant") == sum(matched_cohort_England_2015_2020$migrant_status == "Migrant") # T
save(matched_cohort_England_2015_2020, file = "filepath/matched_cohort_England_2015_2020.Rdata")

# 264644 patients

# Create annual records file
load(file = "filepath/cohort_England_2015_2020_annual_records.Rdata")
matched_cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  filter(patid %in% matched_cohort_England_2015_2020$patid)
n_distinct(matched_cohort_England_2015_2020_annual_records$patid) == nrow(matched_cohort_England_2015_2020) # T

# Check no one 100 years old or over
sum(matched_cohort_England_2015_2020_annual_records$age >= 100)

# Save annual cohort file 
save(matched_cohort_England_2015_2020_annual_records, file = "filepath/matched_cohort_England_2015_2020_annual_records.Rdata")

