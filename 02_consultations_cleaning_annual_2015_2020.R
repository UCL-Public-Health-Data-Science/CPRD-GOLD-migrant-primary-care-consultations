# 0_Description ---------------------------------------------------------------------------

# Migrants' primary care utilisation before and during the COVID-19 pandemic in England: An interrupted time series
# Cleaning of CPRD GOLD consultations for annual consultation outcomes 2015-19
# Date started: 25/03/2020
# Author(s):  Claire Zhang / Yamina Boukari / Neha Pathak / Parth Patel 
# QC (date): Yamina BOukari (27/02/2022)

# 1_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc)

memory.limit(size=1000000000)

# 2_Set working directory ------------------------------------------------------------------

setwd("filepath")

######################## CLEAN RAW ANNUAL CONSULTATION EVENT FILES ############################## ---------------------------------------------------------
# 3_Clean each annual consultation event file ------------------------------------------------

## cons_2015  ------------------------------------------------------------

cons_2015 <- read_csv("filepath") 

# consultations 2015
glimpse(cons_2015)
names(cons_2015)
summary(cons_2015)
head(cons_2015)
tail(cons_2015)

# Check for duplicates observations (whole row)
n_distinct(cons_2015) == count(cons_2015) 

# Remove duplicates of combined patid, eventdate, constype, consid, staffid (likely clinical coding errors as consid links events at the same consultation)
n_distinct(cons_2015) == count(distinct(cons_2015, patid, eventdate, constype, consid, staffid,  .keep_all = TRUE))
distinct_cons_2015 <- cons_2015 %>% distinct(patid, eventdate, constype, consid, staffid,  .keep_all = TRUE) 
nrow(cons_2015) - nrow(distinct_cons_2015) 

# Drop unneeded variables
cons_2015 <-  select(distinct_cons_2015, -c(sysdate, staffid))

# Change variables to the correct data type based on CPRD GOLD data specification & mapping
cons_2015$constype <- factor(cons_2015$constype, levels = c(0:61),
                             labels = c("Data not entered", "Clinic", "Night visit, deputising service", "Follow-up/routine visit", "Night visit, local rota", "Mail from patient", "Night visit, practice",
                                        "Out of hours, practice", "Out of hours, non-practice", "Surgery consultation", "Telephone call from a patient", "Acute visit", "Discharge details",
                                        "Letter from outpatients", "Repeat issue", "Other", "Results recording", "Mail to patient", "Emergency consultation",
                                        "Administration", "Casualty attendance", "Telephone call to a patient", "Third party consultation", "Hospital admission", "Children's home visit",
                                        "Day case report", "GOS18 report", "Home visit", "Hotel visit", "NHS direct report", "Nursing home visit",
                                        "Residential home visit", "Twilight visit", "Triage", "Walk-in centre", "Co-op telephone advice", "Co-op surgery consultation",
                                        "Co-op home visit", "Minor injury service", "Medicine management", "Community clinic", "Community nursing note", "Community nursing report",
                                        "Data transferred from other system", "Health authority entry", "Health visitor note", "Health visitor report", "Hospital inpatient report", "Initial post discharge review",
                                        "Laboratory request", "Night visit", "Radiology request", "Radiology result", "Referral letter", "Social services report",
                                        "Telephone consultation", "Template entry", "GP to GP communication transaction", "Non-consultation medication data", "Non-consultation data", "ePharmacy message", "Extended Hours"))

# Check correct coercion of variable classes
glimpse(cons_2015)

# Check factors correctly labelled
levels(cons_2015$constype)

# Check for missing values
any(is.na(cons_2015)) 
sum(is.na(cons_2015)) 

# Select relevant variables 
cons_2015_clean <- dplyr::select(cons_2015, c(patid, constype, consid, eventdate))

rm(distinct_cons_2015, cons_2015)

# create 4 categories from 61 consultation type
cons_2015_clean <- cons_2015_clean %>%
  mutate(cons4type= constype)
levels(cons_2015_clean$cons4type) <-  c(levels(cons_2015_clean$cons4type),"Non-consultation", "Scheduled consultation", "Unscheduled/out-of-practice consultation", "Phone consultation")
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Clinic"] <- "Scheduled consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Follow-up/routine visit"] <- "Scheduled consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Surgery consultation"] <- "Scheduled consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Repeat issue"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Medicine management"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Community clinic"] <- "Scheduled consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Initial post discharge review"] <- "Scheduled consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Telephone call to a patient"] <- "Phone consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Telephone call from a patient"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Co-op telephone advice"] <- "Phone consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Telephone consultation"] <- "Phone consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Night visit, local rota"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Night visit, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Night visit, deputising service"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Out of hours, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Out of hours, non-practice"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Acute visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Emergency consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Casualty attendance"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Third party consultation"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Hospital admission"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Children's home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Hotel visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Residential home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Twilight visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Walk-in centre"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Co-op surgery consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Co-op home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Minor injury service"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Night visit"] <- "Unscheduled/out-of-practice consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Data not entered"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Mail from patient"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Discharge details"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Letter from outpatients"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Other"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Results recording"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Mail to patient"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Administration"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Day case report"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "GOS18 report"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "NHS direct report"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Triage"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Community nursing note"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Community nursing report"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Data transferred from other system"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Health authority entry"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Health visitor note"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Health visitor report"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Hospital inpatient report"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Laboratory request"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Radiology request"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Radiology result"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Referral letter"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Social services report"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Template entry"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "GP to GP communication transaction"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Non-consultation medication data"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Non-consultation data"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "ePharmacy message"] <- "Non-consultation"
cons_2015_clean$cons4type[cons_2015_clean$cons4type == "Extended Hours"] <- "Scheduled consultation"

# create combine 4 into 2 categories for consultation type
cons_2015_clean <- cons_2015_clean %>%
  mutate(cons2type= cons4type)
levels(cons_2015_clean$cons2type) <-  c(levels(cons_2015_clean$cons2type),"Direct", "Indirect")
cons_2015_clean$cons2type[cons_2015_clean$cons2type == "Non-consultation"] <- "Indirect"
cons_2015_clean$cons2type[cons_2015_clean$cons2type == "Scheduled consultation"] <- "Direct"
cons_2015_clean$cons2type[cons_2015_clean$cons2type == "Unscheduled/out-of-practice consultation"] <- "Direct"
cons_2015_clean$cons2type[cons_2015_clean$cons2type == "Phone consultation"] <- "Direct" 

# Drop unused factor levels
cons_2015_clean$cons4type <- droplevels(cons_2015_clean$cons4type)
cons_2015_clean$cons2type <- droplevels(cons_2015_clean$cons2type)

# Exclude indirect consultations
cons_2015_clean$cons2type <- as.integer(cons_2015_clean$cons2type)
sum(cons_2015_clean$cons2type==2) 
cons_2015_clean <- cons_2015_clean %>%
  filter(cons2type==1)

## Total non-progressive n,% of direct consultation events with each exclusion characteristic

# Exclude events belonging to patients not in the final included cohort
load(file = "primary_care_cleaned_files/01_Cohort_Rdata/cohort_England_2015_2020_annual_records.Rdata")
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  filter(studyyear == 2015)
before_exclusion_events <- nrow(cons_2015_clean)
cons_2015_clean <- cons_2015_clean %>%
  filter(patid %in% cohort_England_2015_2020_annual_records$patid)
total_events_2015 <- nrow(cons_2015_clean)  
before_exclusion_events - total_events_2015 
rm(cohort_England_2015_2020_annual_records)

# Select needed variables
load(file = "primary_care_cleaned_files/01_Cohort_Rdata/cohort_England_2015_2020.Rdata")
data_start_end <- cohort_England_2015_2020 %>%
  select(c(patid, pracid, data_start, data_end))
cons_2015_clean <- cons_2015_clean %>%
  left_join(data_start_end, by = c("patid" = "patid"))

# Events outside of data_start and data_end
outside_data_start_end <- sum(cons_2015_clean$eventdate < cons_2015_clean$data_start | cons_2015_clean$eventdate > cons_2015_clean$data_end)
outside_data_start_end_percent <- outside_data_start_end/total_events_2015 * 100 

# Save excluded events for later assessment of bias
cons_2015_excluded <- cons_2015_clean %>%
  filter(eventdate < data_start | eventdate > data_end)
save(cons_2015_excluded, file = "filepath/cons_2015_excluded.Rdata") 

# Exclude events that are outside data_start and data_end
cons_2015_clean <- cons_2015_clean %>%
  filter(eventdate >= data_start) %>%
  filter(eventdate <= data_end)

# Select needed variables
cons_2015_clean <- select(cons_2015_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate))

# create face-to-face and phone categories
cons_2015_clean <- cons_2015_clean %>%
  mutate(facetoface = cons4type) %>%
  mutate(phone = cons4type)
cons_2015_clean$facetoface <- as.integer(cons_2015_clean$facetoface)
cons_2015_clean$phone <- as.integer(cons_2015_clean$phone)
cons_2015_clean <- cons_2015_clean %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(1,4), 0)) %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(2,3), 1)) %>%
  mutate(phone=replace(phone, phone %in% c(1,2,3),0)) %>%
  mutate(phone=replace(phone, phone %in% c(4),1))

# Drop unused factor levels
cons_2015_clean <- droplevels(cons_2015_clean)

# Check factors correctly labelled
levels(cons_2015_clean$constype)

# Check correct coercion of variable classes
glimpse(cons_2015_clean) 

# Check for NAs
sum(is.na(cons_2015_clean))

# Export as .Rdata
save(cons_2015_clean, file = "filepath") 

rm(list=ls())
gc()

## cons_2016  ------------------------------------------------------------

cons_2016 <- read_csv("filepath") 

# consultations 2016
glimpse(cons_2016)
names(cons_2016)
summary(cons_2016)
head(cons_2016)
tail(cons_2016)

# Check for duplicates observations (whole row)
n_distinct(cons_2016) == count(cons_2016) 

# Remove duplicates of combined patid, eventdate, constype, consid, staffid (likely clinical coding errors as consid links events at the same consultation)
n_distinct(cons_2016) == count(distinct(cons_2016, patid, eventdate, constype, consid, staffid,  .keep_all = TRUE)) 
distinct_cons_2016 <- cons_2016 %>% distinct(patid, eventdate, constype, consid, staffid,  .keep_all = TRUE) 
nrow(cons_2016) - nrow(distinct_cons_2016) 

# Drop unneeded variables
cons_2016 <-  select(distinct_cons_2016, -c(sysdate, staffid))

# Change variables to the correct data type based on CPRD GOLD data specification & mapping
cons_2016$constype <- factor(cons_2016$constype, levels = c(0:61),
                             labels = c("Data not entered", "Clinic", "Night visit, deputising service", "Follow-up/routine visit", "Night visit, local rota", "Mail from patient", "Night visit, practice",
                                        "Out of hours, practice", "Out of hours, non-practice", "Surgery consultation", "Telephone call from a patient", "Acute visit", "Discharge details",
                                        "Letter from outpatients", "Repeat issue", "Other", "Results recording", "Mail to patient", "Emergency consultation",
                                        "Administration", "Casualty attendance", "Telephone call to a patient", "Third party consultation", "Hospital admission", "Children's home visit",
                                        "Day case report", "GOS18 report", "Home visit", "Hotel visit", "NHS direct report", "Nursing home visit",
                                        "Residential home visit", "Twilight visit", "Triage", "Walk-in centre", "Co-op telephone advice", "Co-op surgery consultation",
                                        "Co-op home visit", "Minor injury service", "Medicine management", "Community clinic", "Community nursing note", "Community nursing report",
                                        "Data transferred from other system", "Health authority entry", "Health visitor note", "Health visitor report", "Hospital inpatient report", "Initial post discharge review",
                                        "Laboratory request", "Night visit", "Radiology request", "Radiology result", "Referral letter", "Social services report",
                                        "Telephone consultation", "Template entry", "GP to GP communication transaction", "Non-consultation medication data", "Non-consultation data", "ePharmacy message", "Extended Hours"))

# Check correct coercion of variable classes
glimpse(cons_2016) 

# Check factors correctly labelled
levels(cons_2016$constype)

# Check for missing values
any(is.na(cons_2016)) 
sum(is.na(cons_2016)) 

# Select relevant variables 
cons_2016_clean <- dplyr::select(cons_2016, c(patid, constype, consid, eventdate))

rm(distinct_cons_2016, cons_2016)

# create 4 categories from 61 consultation type
cons_2016_clean <- cons_2016_clean %>%
  mutate(cons4type= constype)
levels(cons_2016_clean$cons4type) <-  c(levels(cons_2016_clean$cons4type),"Non-consultation", "Scheduled consultation", "Unscheduled/out-of-practice consultation", "Phone consultation")
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Clinic"] <- "Scheduled consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Follow-up/routine visit"] <- "Scheduled consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Surgery consultation"] <- "Scheduled consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Repeat issue"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Medicine management"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Community clinic"] <- "Scheduled consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Initial post discharge review"] <- "Scheduled consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Telephone call to a patient"] <- "Phone consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Telephone call from a patient"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Co-op telephone advice"] <- "Phone consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Telephone consultation"] <- "Phone consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Night visit, local rota"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Night visit, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Night visit, deputising service"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Out of hours, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Out of hours, non-practice"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Acute visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Emergency consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Casualty attendance"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Third party consultation"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Hospital admission"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Children's home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Hotel visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Residential home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Twilight visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Walk-in centre"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Co-op surgery consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Co-op home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Minor injury service"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Night visit"] <- "Unscheduled/out-of-practice consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Data not entered"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Mail from patient"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Discharge details"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Letter from outpatients"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Other"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Results recording"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Mail to patient"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Administration"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Day case report"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "GOS18 report"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "NHS direct report"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Triage"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Community nursing note"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Community nursing report"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Data transferred from other system"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Health authority entry"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Health visitor note"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Health visitor report"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Hospital inpatient report"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Laboratory request"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Radiology request"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Radiology result"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Referral letter"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Social services report"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Template entry"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "GP to GP communication transaction"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Non-consultation medication data"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Non-consultation data"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "ePharmacy message"] <- "Non-consultation"
cons_2016_clean$cons4type[cons_2016_clean$cons4type == "Extended Hours"] <- "Scheduled consultation"

# create combine 4 into 2 categories for consultation type
cons_2016_clean <- cons_2016_clean %>%
  mutate(cons2type= cons4type)
levels(cons_2016_clean$cons2type) <-  c(levels(cons_2016_clean$cons2type),"Direct", "Indirect")
cons_2016_clean$cons2type[cons_2016_clean$cons2type == "Non-consultation"] <- "Indirect"
cons_2016_clean$cons2type[cons_2016_clean$cons2type == "Scheduled consultation"] <- "Direct"
cons_2016_clean$cons2type[cons_2016_clean$cons2type == "Unscheduled/out-of-practice consultation"] <- "Direct"
cons_2016_clean$cons2type[cons_2016_clean$cons2type == "Phone consultation"] <- "Direct" 

# Drop unused factor levels
cons_2016_clean$cons4type <- droplevels(cons_2016_clean$cons4type)
cons_2016_clean$cons2type <- droplevels(cons_2016_clean$cons2type)

# Exclude indirect consultations
cons_2016_clean$cons2type <- as.integer(cons_2016_clean$cons2type)
sum(cons_2016_clean$cons2type==2) 
cons_2016_clean <- cons_2016_clean %>%
  filter(cons2type==1)

## Total non-progressive n,% of direct consultation events with each exclusion characteristic

# Exclude events belonging to patients not in the final included cohort
load(file = "filepath/cohort_England_2015_2020_annual_records.Rdata")
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  filter(studyyear == 2016)
before_exclusion_events <- nrow(cons_2016_clean)
cons_2016_clean <- cons_2016_clean %>%
  filter(patid %in% cohort_England_2015_2020_annual_records$patid)
total_events_2016 <- nrow(cons_2016_clean)  
before_exclusion_events - total_events_2016 
rm(cohort_England_2015_2020_annual_records)

# Select needed variables
load(file = "filepath/cohort_England_2015_2020.Rdata")
data_start_end <- cohort_England_2015_2020 %>%
  select(c(patid, pracid, data_start, data_end))
cons_2016_clean <- cons_2016_clean %>%
  left_join(data_start_end, by = c("patid" = "patid"))
cons_2016_clean <- select(cons_2016_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate, data_start, data_end))

# 1574612 events for the included cohort 

# Events outside of data_start and data_end
outside_data_start_end <- sum(cons_2016_clean$eventdate < cons_2016_clean$data_start | cons_2016_clean$eventdate > cons_2016_clean$data_end) 
outside_data_start_end_percent <- outside_data_start_end/total_events_2016 * 100

# Save excluded events for later assessment of bias
cons_2016_excluded <- cons_2016_clean %>%
  filter(eventdate < data_start | eventdate > data_end)
save(cons_2016_excluded, file = "filepath") 

# Exclude events that are outside data_start and data_end
cons_2016_clean <- cons_2016_clean %>%
  filter(eventdate >= data_start) %>%
  filter(eventdate <= data_end)

# Select needed variables
cons_2016_clean <- select(cons_2016_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate))

# create face-to-face and phone categories
cons_2016_clean <- cons_2016_clean %>%
  mutate(facetoface = cons4type) %>%
  mutate(phone = cons4type)
cons_2016_clean$facetoface <- as.integer(cons_2016_clean$facetoface)
cons_2016_clean$phone <- as.integer(cons_2016_clean$phone)
cons_2016_clean <- cons_2016_clean %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(1,4), 0)) %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(2,3), 1)) %>%
  mutate(phone=replace(phone, phone %in% c(1,2,3),0)) %>%
  mutate(phone=replace(phone, phone %in% c(4),1))

# Drop unused factor levels
cons_2016_clean <- droplevels(cons_2016_clean)

# Check factors correctly labelled
levels(cons_2016_clean$constype)

# Check correct coercion of variable classes
glimpse(cons_2016_clean)

# Check for NAs
sum(is.na(cons_2016_clean)) 

# Export as .Rdata
save(cons_2016_clean, file = "filepath") 

rm(list=ls())
gc()


## cons_2017  ------------------------------------------------------------

cons_2017 <- read_csv("filepath") 

# consultations 2017
glimpse(cons_2017)
names(cons_2017)
summary(cons_2017)
head(cons_2017)
tail(cons_2017)

# Check for duplicates observations (whole row)
n_distinct(cons_2017) == count(cons_2017) 

# Remove duplicates of combined patid, eventdate, constype, consid, staffid (likely clinical coding errors as consid links events at the same consultation)
n_distinct(cons_2017) == count(distinct(cons_2017, patid, eventdate, constype, consid, staffid,  .keep_all = TRUE)) 
distinct_cons_2017 <- cons_2017 %>% distinct(patid, eventdate, constype, consid, staffid,  .keep_all = TRUE) 
nrow(cons_2017) - nrow(distinct_cons_2017)

# Drop unneeded variables
cons_2017 <-  select(distinct_cons_2017, -c(sysdate, staffid))

# Change variables to the correct data type based on CPRD GOLD data specification & mapping
cons_2017$constype <- factor(cons_2017$constype, levels = c(0:61),
                             labels = c("Data not entered", "Clinic", "Night visit, deputising service", "Follow-up/routine visit", "Night visit, local rota", "Mail from patient", "Night visit, practice",
                                        "Out of hours, practice", "Out of hours, non-practice", "Surgery consultation", "Telephone call from a patient", "Acute visit", "Discharge details",
                                        "Letter from outpatients", "Repeat issue", "Other", "Results recording", "Mail to patient", "Emergency consultation",
                                        "Administration", "Casualty attendance", "Telephone call to a patient", "Third party consultation", "Hospital admission", "Children's home visit",
                                        "Day case report", "GOS18 report", "Home visit", "Hotel visit", "NHS direct report", "Nursing home visit",
                                        "Residential home visit", "Twilight visit", "Triage", "Walk-in centre", "Co-op telephone advice", "Co-op surgery consultation",
                                        "Co-op home visit", "Minor injury service", "Medicine management", "Community clinic", "Community nursing note", "Community nursing report",
                                        "Data transferred from other system", "Health authority entry", "Health visitor note", "Health visitor report", "Hospital inpatient report", "Initial post discharge review",
                                        "Laboratory request", "Night visit", "Radiology request", "Radiology result", "Referral letter", "Social services report",
                                        "Telephone consultation", "Template entry", "GP to GP communication transaction", "Non-consultation medication data", "Non-consultation data", "ePharmacy message", "Extended Hours"))

# Check correct coercion of variable classes
glimpse(cons_2017) 

# Check factors correctly labelled
levels(cons_2017$constype)

# Check for missing values
any(is.na(cons_2017)) 
sum(is.na(cons_2017)) 

# Select relevant variables 
cons_2017_clean <- dplyr::select(cons_2017, c(patid, constype, consid, eventdate))

rm(distinct_cons_2017, cons_2017)

# create 4 categories from 61 consultation type
cons_2017_clean <- cons_2017_clean %>%
  mutate(cons4type= constype)
levels(cons_2017_clean$cons4type) <-  c(levels(cons_2017_clean$cons4type),"Non-consultation", "Scheduled consultation", "Unscheduled/out-of-practice consultation", "Phone consultation")
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Clinic"] <- "Scheduled consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Follow-up/routine visit"] <- "Scheduled consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Surgery consultation"] <- "Scheduled consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Repeat issue"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Medicine management"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Community clinic"] <- "Scheduled consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Initial post discharge review"] <- "Scheduled consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Telephone call to a patient"] <- "Phone consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Telephone call from a patient"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Co-op telephone advice"] <- "Phone consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Telephone consultation"] <- "Phone consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Night visit, local rota"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Night visit, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Night visit, deputising service"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Out of hours, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Out of hours, non-practice"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Acute visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Emergency consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Casualty attendance"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Third party consultation"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Hospital admission"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Children's home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Hotel visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Residential home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Twilight visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Walk-in centre"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Co-op surgery consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Co-op home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Minor injury service"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Night visit"] <- "Unscheduled/out-of-practice consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Data not entered"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Mail from patient"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Discharge details"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Letter from outpatients"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Other"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Results recording"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Mail to patient"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Administration"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Day case report"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "GOS18 report"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "NHS direct report"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Triage"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Community nursing note"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Community nursing report"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Data transferred from other system"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Health authority entry"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Health visitor note"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Health visitor report"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Hospital inpatient report"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Laboratory request"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Radiology request"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Radiology result"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Referral letter"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Social services report"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Template entry"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "GP to GP communication transaction"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Non-consultation medication data"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Non-consultation data"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "ePharmacy message"] <- "Non-consultation"
cons_2017_clean$cons4type[cons_2017_clean$cons4type == "Extended Hours"] <- "Scheduled consultation"

# create combine 4 into 2 categories for consultation type
cons_2017_clean <- cons_2017_clean %>%
  mutate(cons2type= cons4type)
levels(cons_2017_clean$cons2type) <-  c(levels(cons_2017_clean$cons2type),"Direct", "Indirect")
cons_2017_clean$cons2type[cons_2017_clean$cons2type == "Non-consultation"] <- "Indirect"
cons_2017_clean$cons2type[cons_2017_clean$cons2type == "Scheduled consultation"] <- "Direct"
cons_2017_clean$cons2type[cons_2017_clean$cons2type == "Unscheduled/out-of-practice consultation"] <- "Direct"
cons_2017_clean$cons2type[cons_2017_clean$cons2type == "Phone consultation"] <- "Direct" 

# Drop unused factor levels
cons_2017_clean$cons4type <- droplevels(cons_2017_clean$cons4type)
cons_2017_clean$cons2type <- droplevels(cons_2017_clean$cons2type)

# Exclude indirect consultations
cons_2017_clean$cons2type <- as.integer(cons_2017_clean$cons2type)
sum(cons_2017_clean$cons2type==2) 
cons_2017_clean <- cons_2017_clean %>%
  filter(cons2type==1)

## Total non-progressive n,% of direct consultation events with each exclusion characteristic

# Exclude events belonging to patients not in the final included cohort
load(file = "filepath/cohort_England_2015_2020_annual_records.Rdata")
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  filter(studyyear == 2017)
before_exclusion_events <- nrow(cons_2017_clean)
cons_2017_clean <- cons_2017_clean %>%
  filter(patid %in% cohort_England_2015_2020_annual_records$patid)
total_events_2017 <- nrow(cons_2017_clean)  
before_exclusion_events - total_events_2017 
rm(cohort_England_2015_2020_annual_records)

# Select needed variables
load(file = "filepath/cohort_England_2015_2020.Rdata")
data_start_end <- cohort_England_2015_2020 %>%
  select(c(patid, pracid, data_start, data_end))
cons_2017_clean <- cons_2017_clean %>%
  left_join(data_start_end, by = c("patid" = "patid"))
cons_2017_clean <- select(cons_2017_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate, data_start, data_end))

# 1574612 events for the included cohort 

# Events outside of data_start and data_end
outside_data_start_end <- sum(cons_2017_clean$eventdate < cons_2017_clean$data_start | cons_2017_clean$eventdate > cons_2017_clean$data_end) 
outside_data_start_end_percent <- outside_data_start_end/total_events_2017 * 100

# Save excluded events for later assessment of bias
cons_2017_excluded <- cons_2017_clean %>%
  filter(eventdate < data_start | eventdate > data_end)
save(cons_2017_excluded, file = "filepath") 

# Exclude events that are outside data_start and data_end
cons_2017_clean <- cons_2017_clean %>%
  filter(eventdate >= data_start) %>%
  filter(eventdate <= data_end)

# Select needed variables
cons_2017_clean <- select(cons_2017_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate))

# create face-to-face and phone categories
cons_2017_clean <- cons_2017_clean %>%
  mutate(facetoface = cons4type) %>%
  mutate(phone = cons4type)
cons_2017_clean$facetoface <- as.integer(cons_2017_clean$facetoface)
cons_2017_clean$phone <- as.integer(cons_2017_clean$phone)
cons_2017_clean <- cons_2017_clean %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(1,4), 0)) %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(2,3), 1)) %>%
  mutate(phone=replace(phone, phone %in% c(1,2,3),0)) %>%
  mutate(phone=replace(phone, phone %in% c(4),1))

# Drop unused factor levels
cons_2017_clean <- droplevels(cons_2017_clean)

# Check factors correctly labelled
levels(cons_2017_clean$constype)

# Check correct coercion of variable classes
glimpse(cons_2017_clean) 

# Check for NAs
sum(is.na(cons_2017_clean)) 

# Export as .Rdata
save(cons_2017_clean, file = "filepath") 

rm(list=ls())
gc()



## cons_2018  ------------------------------------------------------------

cons_2018 <- read_csv("filepath") 

# consultations 2018
glimpse(cons_2018)
names(cons_2018)
summary(cons_2018)
head(cons_2018)
tail(cons_2018)

# Check for duplicates observations (whole row)
n_distinct(cons_2018) == count(cons_2018) 

# Remove duplicates of combined patid, eventdate, constype, consid, staffid (likely clinical coding errors as consid links events at the same consultation)
n_distinct(cons_2018) == count(distinct(cons_2018, patid, eventdate, constype, consid, staffid,  .keep_all = TRUE)) 
distinct_cons_2018 <- cons_2018 %>% distinct(patid, eventdate, constype, consid, staffid,  .keep_all = TRUE) 
nrow(cons_2018) - nrow(distinct_cons_2018) 

# Drop unneeded variables
cons_2018 <-  select(distinct_cons_2018, -c(sysdate, staffid))

# Change variables to the correct data type based on CPRD GOLD data specification & mapping
cons_2018$constype <- factor(cons_2018$constype, levels = c(0:61),
                             labels = c("Data not entered", "Clinic", "Night visit, deputising service", "Follow-up/routine visit", "Night visit, local rota", "Mail from patient", "Night visit, practice",
                                        "Out of hours, practice", "Out of hours, non-practice", "Surgery consultation", "Telephone call from a patient", "Acute visit", "Discharge details",
                                        "Letter from outpatients", "Repeat issue", "Other", "Results recording", "Mail to patient", "Emergency consultation",
                                        "Administration", "Casualty attendance", "Telephone call to a patient", "Third party consultation", "Hospital admission", "Children's home visit",
                                        "Day case report", "GOS18 report", "Home visit", "Hotel visit", "NHS direct report", "Nursing home visit",
                                        "Residential home visit", "Twilight visit", "Triage", "Walk-in centre", "Co-op telephone advice", "Co-op surgery consultation",
                                        "Co-op home visit", "Minor injury service", "Medicine management", "Community clinic", "Community nursing note", "Community nursing report",
                                        "Data transferred from other system", "Health authority entry", "Health visitor note", "Health visitor report", "Hospital inpatient report", "Initial post discharge review",
                                        "Laboratory request", "Night visit", "Radiology request", "Radiology result", "Referral letter", "Social services report",
                                        "Telephone consultation", "Template entry", "GP to GP communication transaction", "Non-consultation medication data", "Non-consultation data", "ePharmacy message", "Extended Hours"))

# Check correct coercion of variable classes
glimpse(cons_2018) 

# Check factors correctly labelled
levels(cons_2018$constype)

# Check for missing values
any(is.na(cons_2018)) 
sum(is.na(cons_2018)) 

# Select relevant variables 
cons_2018_clean <- dplyr::select(cons_2018, c(patid, constype, consid, eventdate))

rm(distinct_cons_2018, cons_2018)

# create 4 categories from 61 consultation type
cons_2018_clean <- cons_2018_clean %>%
  mutate(cons4type= constype)
levels(cons_2018_clean$cons4type) <-  c(levels(cons_2018_clean$cons4type),"Non-consultation", "Scheduled consultation", "Unscheduled/out-of-practice consultation", "Phone consultation")
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Clinic"] <- "Scheduled consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Follow-up/routine visit"] <- "Scheduled consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Surgery consultation"] <- "Scheduled consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Repeat issue"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Medicine management"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Community clinic"] <- "Scheduled consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Initial post discharge review"] <- "Scheduled consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Telephone call to a patient"] <- "Phone consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Telephone call from a patient"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Co-op telephone advice"] <- "Phone consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Telephone consultation"] <- "Phone consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Night visit, local rota"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Night visit, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Night visit, deputising service"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Out of hours, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Out of hours, non-practice"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Acute visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Emergency consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Casualty attendance"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Third party consultation"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Hospital admission"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Children's home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Hotel visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Residential home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Twilight visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Walk-in centre"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Co-op surgery consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Co-op home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Minor injury service"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Night visit"] <- "Unscheduled/out-of-practice consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Data not entered"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Mail from patient"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Discharge details"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Letter from outpatients"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Other"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Results recording"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Mail to patient"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Administration"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Day case report"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "GOS18 report"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "NHS direct report"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Triage"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Community nursing note"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Community nursing report"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Data transferred from other system"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Health authority entry"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Health visitor note"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Health visitor report"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Hospital inpatient report"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Laboratory request"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Radiology request"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Radiology result"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Referral letter"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Social services report"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Template entry"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "GP to GP communication transaction"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Non-consultation medication data"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Non-consultation data"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "ePharmacy message"] <- "Non-consultation"
cons_2018_clean$cons4type[cons_2018_clean$cons4type == "Extended Hours"] <- "Scheduled consultation"

# create combine 4 into 2 categories for consultation type
cons_2018_clean <- cons_2018_clean %>%
  mutate(cons2type= cons4type)
levels(cons_2018_clean$cons2type) <-  c(levels(cons_2018_clean$cons2type),"Direct", "Indirect")
cons_2018_clean$cons2type[cons_2018_clean$cons2type == "Non-consultation"] <- "Indirect"
cons_2018_clean$cons2type[cons_2018_clean$cons2type == "Scheduled consultation"] <- "Direct"
cons_2018_clean$cons2type[cons_2018_clean$cons2type == "Unscheduled/out-of-practice consultation"] <- "Direct"
cons_2018_clean$cons2type[cons_2018_clean$cons2type == "Phone consultation"] <- "Direct" 

# Drop unused factor levels
cons_2018_clean$cons4type <- droplevels(cons_2018_clean$cons4type)
cons_2018_clean$cons2type <- droplevels(cons_2018_clean$cons2type)

# Exclude indirect consultations
cons_2018_clean$cons2type <- as.integer(cons_2018_clean$cons2type)
sum(cons_2018_clean$cons2type==2) 
cons_2018_clean <- cons_2018_clean %>%
  filter(cons2type==1)

## Total non-progressive n,% of direct consultation events with each exclusion characteristic

# Exclude events belonging to patients not in the final included cohort
load(file = "filepath")
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  filter(studyyear == 2018)
before_exclusion_events <- nrow(cons_2018_clean)
cons_2018_clean <- cons_2018_clean %>%
  filter(patid %in% cohort_England_2015_2020_annual_records$patid)
total_events_2018 <- nrow(cons_2018_clean)  
before_exclusion_events - total_events_2018 
rm(cohort_England_2015_2020_annual_records)

# Select needed variables
load(file = "filepath")
data_start_end <- cohort_England_2015_2020 %>%
  select(c(patid, pracid, data_start, data_end))
cons_2018_clean <- cons_2018_clean %>%
  left_join(data_start_end, by = c("patid" = "patid"))
cons_2018_clean <- select(cons_2018_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate, data_start, data_end))

# 1012598 events for the included cohort 

# Events outside of data_start and data_end
outside_data_start_end <- sum(cons_2018_clean$eventdate < cons_2018_clean$data_start | cons_2018_clean$eventdate > cons_2018_clean$data_end) 
outside_data_start_end_percent <- outside_data_start_end/total_events_2018 

# Save excluded events for later assessment of bias
cons_2018_excluded <- cons_2018_clean %>%
  filter(eventdate < data_start | eventdate > data_end)
save(cons_2018_excluded, file = "filepath/cons_2018_excluded.Rdata") 

# Exclude events that are outside data_start and data_end
cons_2018_clean <- cons_2018_clean %>%
  filter(eventdate >= data_start) %>%
  filter(eventdate <= data_end)

# Select needed variables
cons_2018_clean <- select(cons_2018_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate))

# create face-to-face and phone categories
cons_2018_clean <- cons_2018_clean %>%
  mutate(facetoface = cons4type) %>%
  mutate(phone = cons4type)
cons_2018_clean$facetoface <- as.integer(cons_2018_clean$facetoface)
cons_2018_clean$phone <- as.integer(cons_2018_clean$phone)
cons_2018_clean <- cons_2018_clean %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(1,4), 0)) %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(2,3), 1)) %>%
  mutate(phone=replace(phone, phone %in% c(1,2,3),0)) %>%
  mutate(phone=replace(phone, phone %in% c(4),1))

# Drop unused factor levels
cons_2018_clean <- droplevels(cons_2018_clean)

# Check factors correctly labelled
levels(cons_2018_clean$constype)

# Check correct coercion of variable classes
glimpse(cons_2018_clean) 

# Check for NAs
sum(is.na(cons_2018_clean)) 

# Export as .Rdata
save(cons_2018_clean, file = "filepath") 

rm(list=ls())
gc()



## cons_2019  ------------------------------------------------------------

cons_2019 <- read_csv("filepath") 

# consultations 2019
glimpse(cons_2019)
names(cons_2019)
summary(cons_2019)
head(cons_2019)
tail(cons_2019)

# Check for duplicates observations (whole row)
n_distinct(cons_2019) == count(cons_2019) 

# Remove duplicates of combined patid, eventdate, constype, consid, staffid (likely clinical coding errors as consid links events at the same consultation)
n_distinct(cons_2019) == count(distinct(cons_2019, patid, eventdate, constype, consid, staffid,  .keep_all = TRUE)) 
distinct_cons_2019 <- cons_2019 %>% distinct(patid, eventdate, constype, consid, staffid,  .keep_all = TRUE) 
nrow(cons_2019) - nrow(distinct_cons_2019) 

# Drop unneeded variables
cons_2019 <-  select(distinct_cons_2019, -c(sysdate, staffid))

# Change variables to the correct data type based on CPRD GOLD data specification & mapping
cons_2019$constype <- factor(cons_2019$constype, levels = c(0:61),
                             labels = c("Data not entered", "Clinic", "Night visit, deputising service", "Follow-up/routine visit", "Night visit, local rota", "Mail from patient", "Night visit, practice",
                                        "Out of hours, practice", "Out of hours, non-practice", "Surgery consultation", "Telephone call from a patient", "Acute visit", "Discharge details",
                                        "Letter from outpatients", "Repeat issue", "Other", "Results recording", "Mail to patient", "Emergency consultation",
                                        "Administration", "Casualty attendance", "Telephone call to a patient", "Third party consultation", "Hospital admission", "Children's home visit",
                                        "Day case report", "GOS18 report", "Home visit", "Hotel visit", "NHS direct report", "Nursing home visit",
                                        "Residential home visit", "Twilight visit", "Triage", "Walk-in centre", "Co-op telephone advice", "Co-op surgery consultation",
                                        "Co-op home visit", "Minor injury service", "Medicine management", "Community clinic", "Community nursing note", "Community nursing report",
                                        "Data transferred from other system", "Health authority entry", "Health visitor note", "Health visitor report", "Hospital inpatient report", "Initial post discharge review",
                                        "Laboratory request", "Night visit", "Radiology request", "Radiology result", "Referral letter", "Social services report",
                                        "Telephone consultation", "Template entry", "GP to GP communication transaction", "Non-consultation medication data", "Non-consultation data", "ePharmacy message", "Extended Hours"))

# Check correct coercion of variable classes
glimpse(cons_2019) 

# Check factors correctly labelled
levels(cons_2019$constype)

# Check for missing values
any(is.na(cons_2019)) 
sum(is.na(cons_2019)) 

# Select relevant variables 
cons_2019_clean <- dplyr::select(cons_2019, c(patid, constype, consid, eventdate))

rm(distinct_cons_2019, cons_2019)

# create 4 categories from 61 consultation type
cons_2019_clean <- cons_2019_clean %>%
  mutate(cons4type= constype)
levels(cons_2019_clean$cons4type) <-  c(levels(cons_2019_clean$cons4type),"Non-consultation", "Scheduled consultation", "Unscheduled/out-of-practice consultation", "Phone consultation")
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Clinic"] <- "Scheduled consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Follow-up/routine visit"] <- "Scheduled consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Surgery consultation"] <- "Scheduled consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Repeat issue"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Medicine management"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Community clinic"] <- "Scheduled consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Initial post discharge review"] <- "Scheduled consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Telephone call to a patient"] <- "Phone consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Telephone call from a patient"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Co-op telephone advice"] <- "Phone consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Telephone consultation"] <- "Phone consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Night visit, local rota"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Night visit, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Night visit, deputising service"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Out of hours, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Out of hours, non-practice"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Acute visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Emergency consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Casualty attendance"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Third party consultation"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Hospital admission"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Children's home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Hotel visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Residential home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Twilight visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Walk-in centre"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Co-op surgery consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Co-op home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Minor injury service"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Night visit"] <- "Unscheduled/out-of-practice consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Data not entered"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Mail from patient"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Discharge details"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Letter from outpatients"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Other"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Results recording"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Mail to patient"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Administration"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Day case report"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "GOS18 report"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "NHS direct report"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Triage"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Community nursing note"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Community nursing report"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Data transferred from other system"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Health authority entry"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Health visitor note"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Health visitor report"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Hospital inpatient report"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Laboratory request"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Radiology request"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Radiology result"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Referral letter"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Social services report"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Template entry"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "GP to GP communication transaction"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Non-consultation medication data"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Non-consultation data"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "ePharmacy message"] <- "Non-consultation"
cons_2019_clean$cons4type[cons_2019_clean$cons4type == "Extended Hours"] <- "Scheduled consultation"

# create combine 4 into 2 categories for consultation type
cons_2019_clean <- cons_2019_clean %>%
  mutate(cons2type= cons4type)
levels(cons_2019_clean$cons2type) <-  c(levels(cons_2019_clean$cons2type),"Direct", "Indirect")
cons_2019_clean$cons2type[cons_2019_clean$cons2type == "Non-consultation"] <- "Indirect"
cons_2019_clean$cons2type[cons_2019_clean$cons2type == "Scheduled consultation"] <- "Direct"
cons_2019_clean$cons2type[cons_2019_clean$cons2type == "Unscheduled/out-of-practice consultation"] <- "Direct"
cons_2019_clean$cons2type[cons_2019_clean$cons2type == "Phone consultation"] <- "Direct" 

# Drop unused factor levels
cons_2019_clean$cons4type <- droplevels(cons_2019_clean$cons4type)
cons_2019_clean$cons2type <- droplevels(cons_2019_clean$cons2type)

# Exclude indirect consultations
cons_2019_clean$cons2type <- as.integer(cons_2019_clean$cons2type)
sum(cons_2019_clean$cons2type==2) 
cons_2019_clean <- cons_2019_clean %>%
  filter(cons2type==1)

## Total non-progressive n,% of direct consultation events with each exclusion characteristic

# Exclude events belonging to patients not in the final included cohort
load(file = "filepath")
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  filter(studyyear == 2019)
before_exclusion_events <- nrow(cons_2019_clean)
cons_2019_clean <- cons_2019_clean %>%
  filter(patid %in% cohort_England_2015_2020_annual_records$patid)
total_events_2019 <- nrow(cons_2019_clean)  
before_exclusion_events - total_events_2019 
rm(cohort_England_2015_2020_annual_records)

# Select needed variables
load(file = "filepath")
data_start_end <- cohort_England_2015_2020 %>%
  select(c(patid, pracid, data_start, data_end))
cons_2019_clean <- cons_2019_clean %>%
  left_join(data_start_end, by = c("patid" = "patid"))
cons_2019_clean <- select(cons_2019_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate, data_start, data_end))

# 880091 events for the included cohort 

# Events outside of data_start and data_end
outside_data_start_end <- sum(cons_2019_clean$eventdate < cons_2019_clean$data_start | cons_2019_clean$eventdate > cons_2019_clean$data_end) 
outside_data_start_end_percent <- outside_data_start_end/total_events_2019 

# Save excluded events for later assessment of bias
cons_2019_excluded <- cons_2019_clean %>%
  filter(eventdate < data_start | eventdate > data_end)
save(cons_2019_excluded, file = "filepath/cons_2019_excluded.Rdata") 

# Exclude events that are outside data_start and data_end
cons_2019_clean <- cons_2019_clean %>%
  filter(eventdate >= data_start) %>%
  filter(eventdate <= data_end)

# Select needed variables
cons_2019_clean <- select(cons_2019_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate))

# create face-to-face and phone categories
cons_2019_clean <- cons_2019_clean %>%
  mutate(facetoface = cons4type) %>%
  mutate(phone = cons4type)
cons_2019_clean$facetoface <- as.integer(cons_2019_clean$facetoface)
cons_2019_clean$phone <- as.integer(cons_2019_clean$phone)
cons_2019_clean <- cons_2019_clean %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(1,4), 0)) %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(2,3), 1)) %>%
  mutate(phone=replace(phone, phone %in% c(1,2,3),0)) %>%
  mutate(phone=replace(phone, phone %in% c(4),1))

# Drop unused factor levels
cons_2019_clean <- droplevels(cons_2019_clean)

# Check factors correctly labelled
levels(cons_2019_clean$constype)

# Check correct coercion of variable classes
glimpse(cons_2019_clean) 

# Check for NAs
sum(is.na(cons_2019_clean)) 

# Export as .Rdata
save(cons_2019_clean, file = "filepath") 

rm(list=ls())
gc()



## cons_2020  ------------------------------------------------------------

cons_2020 <- read_csv("filepath") 

# consultations 2020
glimpse(cons_2020)
names(cons_2020)
summary(cons_2020)
head(cons_2020)
tail(cons_2020)

# Check for duplicates observations (whole row)
n_distinct(cons_2020) == count(cons_2020) 

# Remove duplicates of combined patid, eventdate, constype, consid, staffid (likely clinical coding errors as consid links events at the same consultation)
n_distinct(cons_2020) == count(distinct(cons_2020, patid, eventdate, constype, consid, staffid,  .keep_all = TRUE)) 
distinct_cons_2020 <- cons_2020 %>% distinct(patid, eventdate, constype, consid, staffid,  .keep_all = TRUE) 
nrow(cons_2020) - nrow(distinct_cons_2020)  

# Drop unneeded variables
cons_2020 <-  select(distinct_cons_2020, -c(sysdate, staffid))

# Change variables to the correct data type based on CPRD GOLD data specification & mapping
cons_2020$constype <- factor(cons_2020$constype, levels = c(0:61),
                             labels = c("Data not entered", "Clinic", "Night visit, deputising service", "Follow-up/routine visit", "Night visit, local rota", "Mail from patient", "Night visit, practice",
                                        "Out of hours, practice", "Out of hours, non-practice", "Surgery consultation", "Telephone call from a patient", "Acute visit", "Discharge details",
                                        "Letter from outpatients", "Repeat issue", "Other", "Results recording", "Mail to patient", "Emergency consultation",
                                        "Administration", "Casualty attendance", "Telephone call to a patient", "Third party consultation", "Hospital admission", "Children's home visit",
                                        "Day case report", "GOS18 report", "Home visit", "Hotel visit", "NHS direct report", "Nursing home visit",
                                        "Residential home visit", "Twilight visit", "Triage", "Walk-in centre", "Co-op telephone advice", "Co-op surgery consultation",
                                        "Co-op home visit", "Minor injury service", "Medicine management", "Community clinic", "Community nursing note", "Community nursing report",
                                        "Data transferred from other system", "Health authority entry", "Health visitor note", "Health visitor report", "Hospital inpatient report", "Initial post discharge review",
                                        "Laboratory request", "Night visit", "Radiology request", "Radiology result", "Referral letter", "Social services report",
                                        "Telephone consultation", "Template entry", "GP to GP communication transaction", "Non-consultation medication data", "Non-consultation data", "ePharmacy message", "Extended Hours"))

# Check correct coercion of variable classes
glimpse(cons_2020) 

# Check factors correctly labelled
levels(cons_2020$constype)

# Check for missing values
any(is.na(cons_2020)) 
sum(is.na(cons_2020)) 

# Select relevant variables 
cons_2020_clean <- dplyr::select(cons_2020, c(patid, constype, consid, eventdate))

rm(distinct_cons_2020, cons_2020)

# create 4 categories from 61 consultation type
cons_2020_clean <- cons_2020_clean %>%
  mutate(cons4type= constype)
levels(cons_2020_clean$cons4type) <-  c(levels(cons_2020_clean$cons4type),"Non-consultation", "Scheduled consultation", "Unscheduled/out-of-practice consultation", "Phone consultation")
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Clinic"] <- "Scheduled consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Follow-up/routine visit"] <- "Scheduled consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Surgery consultation"] <- "Scheduled consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Repeat issue"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Medicine management"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Community clinic"] <- "Scheduled consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Initial post discharge review"] <- "Scheduled consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Telephone call to a patient"] <- "Phone consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Telephone call from a patient"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Co-op telephone advice"] <- "Phone consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Telephone consultation"] <- "Phone consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Night visit, local rota"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Night visit, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Night visit, deputising service"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Out of hours, practice"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Out of hours, non-practice"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Acute visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Emergency consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Casualty attendance"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Third party consultation"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Hospital admission"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Children's home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Hotel visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Residential home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Nursing home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Twilight visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Walk-in centre"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Co-op surgery consultation"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Co-op home visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Minor injury service"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Night visit"] <- "Unscheduled/out-of-practice consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Data not entered"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Mail from patient"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Discharge details"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Letter from outpatients"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Other"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Results recording"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Mail to patient"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Administration"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Day case report"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "GOS18 report"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "NHS direct report"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Triage"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Community nursing note"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Community nursing report"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Data transferred from other system"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Health authority entry"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Health visitor note"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Health visitor report"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Hospital inpatient report"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Laboratory request"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Radiology request"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Radiology result"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Referral letter"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Social services report"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Template entry"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "GP to GP communication transaction"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Non-consultation medication data"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Non-consultation data"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "ePharmacy message"] <- "Non-consultation"
cons_2020_clean$cons4type[cons_2020_clean$cons4type == "Extended Hours"] <- "Scheduled consultation"

# create combine 4 into 2 categories for consultation type
cons_2020_clean <- cons_2020_clean %>%
  mutate(cons2type= cons4type)
levels(cons_2020_clean$cons2type) <-  c(levels(cons_2020_clean$cons2type),"Direct", "Indirect")
cons_2020_clean$cons2type[cons_2020_clean$cons2type == "Non-consultation"] <- "Indirect"
cons_2020_clean$cons2type[cons_2020_clean$cons2type == "Scheduled consultation"] <- "Direct"
cons_2020_clean$cons2type[cons_2020_clean$cons2type == "Unscheduled/out-of-practice consultation"] <- "Direct"
cons_2020_clean$cons2type[cons_2020_clean$cons2type == "Phone consultation"] <- "Direct" 

# Drop unused factor levels
cons_2020_clean$cons4type <- droplevels(cons_2020_clean$cons4type)
cons_2020_clean$cons2type <- droplevels(cons_2020_clean$cons2type)

# Exclude indirect consultations
cons_2020_clean$cons2type <- as.integer(cons_2020_clean$cons2type)
sum(cons_2020_clean$cons2type==2) 
cons_2020_clean <- cons_2020_clean %>%
  filter(cons2type==1)

## Total non-progressive n,% of direct consultation events with each exclusion characteristic

# Exclude events belonging to patients not in the final included cohort
load(file = "filepath")
cohort_England_2015_2020_annual_records <- cohort_England_2015_2020_annual_records %>%
  filter(studyyear == 2020)
before_exclusion_events <- nrow(cons_2020_clean)
cons_2020_clean <- cons_2020_clean %>%
  filter(patid %in% cohort_England_2015_2020_annual_records$patid)
total_events_2020 <- nrow(cons_2020_clean)  
before_exclusion_events - total_events_2020 
rm(cohort_England_2015_2020_annual_records)

# Select needed variables
load(file = "filepath")
data_start_end <- cohort_England_2015_2020 %>%
  select(c(patid, pracid, data_start, data_end))
cons_2020_clean <- cons_2020_clean %>%
  left_join(data_start_end, by = c("patid" = "patid"))
cons_2020_clean <- select(cons_2020_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate, data_start, data_end))

# Events outside of data_start and data_end
outside_data_start_end <- sum(cons_2020_clean$eventdate < cons_2020_clean$data_start | cons_2020_clean$eventdate > cons_2020_clean$data_end) 
outside_data_start_end_percent <- outside_data_start_end/total_events_2020 

# Save excluded events for later assessment of bias
cons_2020_excluded <- cons_2020_clean %>%
  filter(eventdate < data_start | eventdate > data_end)
save(cons_2020_excluded, file = "filepath") 

# Exclude events that are outside data_start and data_end
cons_2020_clean <- cons_2020_clean %>%
  filter(eventdate >= data_start) %>%
  filter(eventdate <= data_end)

# Select needed variables
cons_2020_clean <- select(cons_2020_clean, c(patid, pracid, constype, cons4type, cons2type,
                                             consid, eventdate))

# create face-to-face and phone categories
cons_2020_clean <- cons_2020_clean %>%
  mutate(facetoface = cons4type) %>%
  mutate(phone = cons4type)
cons_2020_clean$facetoface <- as.integer(cons_2020_clean$facetoface)
cons_2020_clean$phone <- as.integer(cons_2020_clean$phone)
cons_2020_clean <- cons_2020_clean %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(1,4), 0)) %>%
  mutate(facetoface=replace(facetoface, facetoface %in% c(2,3), 1)) %>%
  mutate(phone=replace(phone, phone %in% c(1,2,3),0)) %>%
  mutate(phone=replace(phone, phone %in% c(4),1))

# Drop unused factor levels
cons_2020_clean <- droplevels(cons_2020_clean)

# Check factors correctly labelled
levels(cons_2020_clean$constype)

# Check correct coercion of variable classes
glimpse(cons_2020_clean) 

# Check for NAs
sum(is.na(cons_2020_clean)) 

# Export as .Rdata
save(cons_2020_clean, file = "filepath") 

rm(list=ls())
gc()



######################## UNMATCHED ANNUAL CONSULATION COUNT DATASETS ############################## ---------------------------------------------------------


# 4_Create annual conscount variables (2015-2020) -----------------------------------------------------

# cons_2015_clean ---

load(file = "filepath")

# One row per patient
cons_2015_final <- cons_2015_clean %>%
  group_by(patid) %>%
  summarise(conscount=sum(cons2type),
            phone=sum(phone),
            facetoface=sum(facetoface))

# Create studyyear variable
cons_2015_final$studyyear <- 2015

glimpse(cons_2015_final)

n_distinct(cons_2015_final$patid) 

# Drop unused factor levels
cons_2015_final <- droplevels(cons_2015_final)

# Save data
save(cons_2015_final, file = "filepath")

# Remove uneeded data from environment
rm(list = ls(all.names = TRUE))
gc()

# cons_2016_clean ---

load(file = "filepath")

# One row per patient
cons_2016_final <- cons_2016_clean %>%
  group_by(patid) %>%
  summarise(conscount=sum(cons2type),
            phone=sum(phone),
            facetoface=sum(facetoface))

# Create studyyear variable
cons_2016_final$studyyear <- 2016

glimpse(cons_2016_final)

n_distinct(cons_2016_final$patid) 

# Drop unused factor levels
cons_2016_final <- droplevels(cons_2016_final)

# Save data
save(cons_2016_final, file = "filepath")

# Remove uneeded data from environment
rm(list = ls(all.names = TRUE))
gc()

# cons_2017_clean ---

load(file = "filepath")

# One row per patient
cons_2017_final <- cons_2017_clean %>%
  group_by(patid) %>%
  summarise(conscount=sum(cons2type),
            phone=sum(phone),
            facetoface=sum(facetoface))

# Create studyyear variable
cons_2017_final$studyyear <- 2017

glimpse(cons_2017_final)

n_distinct(cons_2017_final$patid) 

# Drop unused factor levels
cons_2017_final <- droplevels(cons_2017_final)

# Save data
save(cons_2017_final, file = "filepath")

# Remove uneeded data from environment
rm(list = ls(all.names = TRUE))
gc()

# cons_2018_clean ---

load(file = "filepath")

# One row per patient
cons_2018_final <- cons_2018_clean %>%
  group_by(patid) %>%
  summarise(conscount=sum(cons2type),
            phone=sum(phone),
            facetoface=sum(facetoface))

# Create studyyear variable
cons_2018_final$studyyear <- 2018

glimpse(cons_2018_final)

n_distinct(cons_2018_final$patid) 

# Drop unused factor levels
cons_2018_final <- droplevels(cons_2018_final)

# Save data
save(cons_2018_final, file = "filepath")

# Remove uneeded data from environment
rm(list = ls(all.names = TRUE))
gc()

# cons_2019_clean ---

load(file = "filepath")

# One row per patient
cons_2019_final <- cons_2019_clean %>%
  group_by(patid) %>%
  summarise(conscount=sum(cons2type),
            phone=sum(phone),
            facetoface=sum(facetoface))

# Create studyyear variable
cons_2019_final$studyyear <- 2019

glimpse(cons_2019_final)

n_distinct(cons_2019_final$patid)  

# Drop unused factor levels
cons_2019_final <- droplevels(cons_2019_final)

# Save data
save(cons_2019_final, file = "filepath")

# Remove uneeded data from environment
rm(list = ls(all.names = TRUE))
gc()

# cons_2020_clean ---

load(file = "filepath")

# One row per patient
cons_2020_final <- cons_2020_clean %>%
  group_by(patid) %>%
  summarise(conscount=sum(cons2type),
            phone=sum(phone),
            facetoface=sum(facetoface))

# Create studyyear variable
cons_2020_final$studyyear <- 2020

glimpse(cons_2020_final)

n_distinct(cons_2020_final$patid) 

# Drop unused factor levels
cons_2020_final <- droplevels(cons_2020_final)

# Save data
save(cons_2020_final, file = "filepath")

# Remove uneeded data from environment
rm(list = ls(all.names = TRUE))
gc()

# 5_Join all annual consultation event files -------------------------------------

# Load files
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")

# Join all annual files into full file
cons_allyears <- rbind(cons_2015_final,cons_2016_final,cons_2017_final,cons_2018_final,cons_2019_final,cons_2020_final)

# Coerce to correct variable types and check levels
glimpse(cons_allyears)

# Check number of consultations
sum(cons_allyears$conscount) 
any(is.na(cons_allyears))

# Check number of distinct patients in cons_allyears
n_distinct(cons_allyears$patid) 

# Save data
save(cons_allyears, file = "filepath")

# Remove unneeded data from environment
rm(cons_2015_final, cons_2016_final, cons_2017_final,cons_2018_final,cons_2019_final,
   cons_2020_final)

# 6_Join consultations onto annual records file ------

# Load cohort file and join
load(file = "filepath")
cohort_England_2015_2020_annual_conscounts <- cohort_England_2015_2020_annual_records %>%
  left_join(cons_allyears, by = c("patid" = "patid", "studyyear" = "studyyear"))

# Check distinct number of patients matches the cohort file
load(file = "filepath")
n_distinct(cohort_England_2015_2020_annual_conscounts$patid) == n_distinct(cohort_England_2015_2020$patid) 

# Replace NAs in cons count variables with 0
cohort_England_2015_2020_annual_conscounts$conscount <- cohort_England_2015_2020_annual_conscounts$conscount %>% replace_na(0)
sum(is.na(cohort_England_2015_2020_annual_conscounts$conscount)) 
glimpse(cohort_England_2015_2020_annual_conscounts)

# Check consocunt numbers match
sum(cohort_England_2015_2020_annual_conscounts$conscount)
sum(cohort_England_2015_2020_annual_conscounts$conscount) == sum(cons_allyears$conscount) 

save(cohort_England_2015_2020_annual_conscounts, file = "filepath") 

######################## MATCHED ANNUAL CONSULATION COUNT DATASETS  ############################## ---------------------------------------------------------

# 7_Sensitivity analysis 1: Create exact-matched annual cohort dataset (matched on age_data_start, year_data_start and prac_region) -----

load(file = "filepath")
load(file = "filepath")

matched_cohort_England_2015_2020_annual_conscounts_test_ds <- cohort_England_2015_2020_annual_conscounts %>%
  semi_join(matched_cohort_England_2015_2020_test_ds, by = 'patid')

load(file = "filepath")

# check number of rows
nrow(matched_cohort_England_2015_2020_annual_records_test_ds) == nrow(matched_cohort_England_2015_2020_annual_conscounts_test_ds) 

# Check conscounts and distinct patients
sum(matched_cohort_England_2015_2020_annual_conscounts_test_ds$conscount)
n_distinct(matched_cohort_England_2015_2020_annual_conscounts_test_ds$patid) ==n_distinct(matched_cohort_England_2015_2020_test_ds$patid) 
save(matched_cohort_England_2015_2020_annual_conscounts_test_ds, file = "filepath/matched_cohort_England_2015_2020_annual_conscounts_test_ds.Rdata") 

# 8_Sensitivity analysis 2: Create exact-matched annual cohort dataset (matched on follow-up and prac_region) -----

load(file = "filepath")
load(file = "filepath")

matched_cohort_England_2015_2020_annual_conscounts_fu <- cohort_England_2015_2020_annual_conscounts %>%
  semi_join(matched_cohort_England_2015_2020_fu, by = 'patid')

load(file = "filepath")

# check number of rows
nrow(matched_cohort_England_2015_2020_annual_records_fu) == nrow(matched_cohort_England_2015_2020_annual_conscounts_fu) 

# Check conscounts and distinct patients
sum(matched_cohort_England_2015_2020_annual_conscounts_fu$conscount)
n_distinct(matched_cohort_England_2015_2020_annual_conscounts_fu$patid) ==n_distinct(matched_cohort_England_2015_2020_fu$patid) 

save(matched_cohort_England_2015_2020_annual_conscounts_fu, file = "filepath") 

# 9_Not used in paper ----
## 9a_Create exact-matched annual cohort dataset (matched on age at cohort entry, year at cohort entry and prac_region) -----

load(file = "filepath")
load(file = "filepath")

matched_cohort_England_2015_2020_annual_conscounts <- cohort_England_2015_2020_annual_conscounts %>%
  semi_join(matched_cohort_England_2015_2020, by = 'patid')

load(file = "filepath")

# check number of rows
nrow(matched_cohort_England_2015_2020_annual_records) == nrow(matched_cohort_England_2015_2020_annual_conscounts) # T

# Create new age variables
# Create year of cohort entry variable
matched_cohort_England_2015_2020_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  mutate(year_cohort_entry = year(cohort_entry))

# Create age of cohort entry variables
calc_age <- function(birthDate, refDate = Sys.Date()) {
  require(lubridate)
  period <- as.period(interval(birthDate, refDate),unit = "year")
  period$year
}

matched_cohort_England_2015_2020_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  mutate(age_cohort_entry = calc_age(dob, cohort_entry))

matched_cohort_England_2015_2020_annual_conscounts <-matched_cohort_England_2015_2020_annual_conscounts %>%
  mutate(cohort_entry_agecat = age_cohort_entry)

matched_cohort_England_2015_2020_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, cohort_entry_agecat <6, 1)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 6, 10), 2)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 11, 15), 3)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 16, 19), 4)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 20, 24), 5)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 25, 29), 6)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 30, 34), 7)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 35, 39), 8)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 40, 44), 9)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 45, 49), 10)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 50, 54), 11)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 55, 59), 12)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 60, 64), 13)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 65, 69), 14)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 70, 74), 15)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 75, 79), 16)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 80, 84), 17)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 85, 89), 18)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 90, 94), 19)) %>%
  mutate(cohort_entry_agecat=replace(cohort_entry_agecat, between (cohort_entry_agecat, 95, 99), 20)) 

matched_cohort_England_2015_2020_annual_conscounts$cohort_entry_agecat <- factor(matched_cohort_England_2015_2020_annual_conscounts$cohort_entry_agecat, 
                                                                                 levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
                                                                                 labels = c("0-5 years", "6-10 years", "11-15 years",
                                                                                            "16-19 years", "20-24 years", "25-29 years",
                                                                                            "30-34 years", "35-39 years", "40-44 years",
                                                                                            "45-49 years", "50-54 years", "55-59 years",
                                                                                            "60-64 years", "65-69 years", "70-74 years",
                                                                                            "75-79 years", "80-84 years", "85-89 years",
                                                                                            "90-94 years", "95-99 years"))

# Create age subcohort of cohort entry variable

matched_cohort_England_2015_2020_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  mutate(cohort_entry_age_subcohort = age_cohort_entry)

matched_cohort_England_2015_2020_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  mutate(cohort_entry_age_subcohort=replace(cohort_entry_age_subcohort, cohort_entry_age_subcohort <16, 1)) %>%
  mutate(cohort_entry_age_subcohort=replace(cohort_entry_age_subcohort, between (cohort_entry_age_subcohort, 16, 24), 2)) %>%
  mutate(cohort_entry_age_subcohort=replace(cohort_entry_age_subcohort, between (cohort_entry_age_subcohort, 25, 34), 3)) %>%
  mutate(cohort_entry_age_subcohort=replace(cohort_entry_age_subcohort, between (cohort_entry_age_subcohort, 35, 49), 4)) %>%
  mutate(cohort_entry_age_subcohort=replace(cohort_entry_age_subcohort, between (cohort_entry_age_subcohort, 50, 64), 5)) %>%
  mutate(cohort_entry_age_subcohort=replace(cohort_entry_age_subcohort, cohort_entry_age_subcohort >= 65, 6)) 

matched_cohort_England_2015_2020_annual_conscounts$cohort_entry_age_subcohort <- factor(matched_cohort_England_2015_2020_annual_conscounts$cohort_entry_age_subcohort,levels = c(1,2,3,4,5,6),
                                                                                        labels = c("0-15 years", "16-24 years", "25-34 years",
                                                                                                   "35-49 years","50-64 years", ">=65 years"))

glimpse(matched_cohort_England_2015_2020_annual_conscounts)

# Check conscounts and distinct patients
sum(matched_cohort_England_2015_2020_annual_conscounts$conscount)
n_distinct(matched_cohort_England_2015_2020_annual_conscounts$patid) ==n_distinct(matched_cohort_England_2015_2020$patid) # T

save(matched_cohort_England_2015_2020_annual_conscounts, file = "filepath") 



