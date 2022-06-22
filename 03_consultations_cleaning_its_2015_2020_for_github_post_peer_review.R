# 0_Description ---------------------------------------------------------------------------

# Cleaning of CPRD GOLD consultations for weekly consultation outcomes 2015-20 (ITS analysis)
# Date started: 25/03/2021
# Author(s): Claire Zhang / Yamina Boukari
# QC (date): Yamina Boukari (27/02/2022)

# 1_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc)

# 2_Set working directory ------------------------------------------------------------------

setwd("filepath")

######################### CREATE WEEKLY CONSULTATION COUNTS ######################## ---------------------------------------------------------
# 3_Create weekly conscount variables (2015-2020) -----------------------------------------------------

load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")

# Join all files into full file
cons_allyears_clean <- rbind(cons_2015_clean, cons_2016_clean, cons_2017_clean, 
                             cons_2018_clean, cons_2019_clean, cons_2020_clean)

summary(cons_allyears_clean$eventdate)

rm(cons_2015_clean, cons_2016_clean, cons_2017_clean, cons_2018_clean, cons_2019_clean, cons_2020_clean)

# Create week variable (weeks start on Sunday, end on Saturday)
cons_2015_2020_weekly_conscounts <- cons_allyears_clean %>%
  mutate(week = cut.Date(eventdate, breaks = "1 week", labels = FALSE, start.on.monday = FALSE)) %>%
  arrange(eventdate)

rm(cons_allyears_clean)

# One row per patient
cons_2015_2020_weekly_conscounts <- cons_2015_2020_weekly_conscounts %>%
  group_by(patid, week) %>%
  summarise(conscount=sum(cons2type),
            phone=sum(phone),
            facetoface=sum(facetoface))

glimpse(cons_2015_2020_weekly_conscounts)

n_distinct(cons_2015_2020_weekly_conscounts$patid) 

# Drop unused factor levels
cons_2015_2020_weekly_conscounts <- droplevels(cons_2015_2020_weekly_conscounts)

# Save data
save(cons_2015_2020_weekly_conscounts, file = "filepath")

# 4_Join weekly conscount variables onto cohort weekly records file  ------

load(file = "filepath")

# Load cohort file and join
load(file = "filepath")

cohort_England_2015_2020_weekly_conscounts <- matched_its_weekly_records %>%
  left_join(cons_2015_2020_weekly_conscounts, by = c("patid" = "patid", "studyweek" = "week"))

rm(matched_its_weekly_records,cons_2015_2020_weekly_conscounts)

# Replace NAs in cons count variables with 0
cohort_England_2015_2020_weekly_conscounts$conscount <- cohort_England_2015_2020_weekly_conscounts$conscount %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts$phone <- cohort_England_2015_2020_weekly_conscounts$phone %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts$facetoface <- cohort_England_2015_2020_weekly_conscounts$facetoface %>% replace_na(0)
sum(is.na(cohort_England_2015_2020_weekly_conscounts$conscount))
glimpse(cohort_England_2015_2020_weekly_conscounts)
sum(cohort_England_2015_2020_weekly_conscounts$conscount)

# Save data
save(cohort_England_2015_2020_weekly_conscounts, file = "filepath")

rm(list = ls()[!(ls() %in% c('cohort_England_2015_2020_weekly_conscounts_England'))])
gc()

#################### CREATE WEEKLY AGGREGATE DATASETS FOR ITS ANALYSIS #################### ---------------------------------------------------------

# 5_Create aggregate weekly data (migrant status)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts <- cohort_England_2015_2020_weekly_conscounts %>%
  group_by(migrant_status, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            pweeks=sum(pweeks),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts <- left_join(England_2015_2020_aggregate_weekly_conscounts, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts <- England_2015_2020_aggregate_weekly_conscounts %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts <- England_2015_2020_aggregate_weekly_conscounts %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Check no pweeks are greater than popsize
any(England_2015_2020_aggregate_weekly_conscounts$pweeks > England_2015_2020_aggregate_weekly_conscounts$popsize)

# Save
save(England_2015_2020_aggregate_weekly_conscounts, file = "filepath")

rm(list=ls())
gc()

# 6_Create aggregate weekly data (migrant status age subcohort)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_agesubcohort <- cohort_England_2015_2020_weekly_conscounts %>%
  group_by(migrant_status, studyweek, age_subcohort) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_agesubcohort <- left_join(England_2015_2020_aggregate_weekly_conscounts_agesubcohort, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_agesubcohort$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_agesubcohort <- England_2015_2020_aggregate_weekly_conscounts_agesubcohort %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_agesubcohort <- England_2015_2020_aggregate_weekly_conscounts_agesubcohort %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_agesubcohort, file = "filepath")

rm(list=ls())
gc()

# 7_Create aggregate weekly data (ethnicity - 6 categories)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_ethnicity <- cohort_England_2015_2020_weekly_conscounts %>%
  group_by(ethnicat6, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_ethnicity <- left_join(England_2015_2020_aggregate_weekly_conscounts_ethnicity, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_ethnicity$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_ethnicity %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_ethnicity %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_ethnicity, file = "filepath")

rm(list=ls())
gc()

# 8_Create aggregate weekly data (migrant_status and ethnicity 6 categories)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- cohort_England_2015_2020_weekly_conscounts %>%
  group_by(migrant_status, studyweek, ethnicat6) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- left_join(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, file = "filepath")

rm(list=ls())
gc()

# 9_Create aggregate weekly data (migrant_status and ethnicity 18 categories)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat <- cohort_England_2015_2020_weekly_conscounts %>%
  group_by(migrant_status, studyweek, ethnicat) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat <- left_join(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat, file = "filepath")

rm(list=ls())
gc()


# 10_Create aggregate weekly data (certainty of migration)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_migcertainty <- cohort_England_2015_2020_weekly_conscounts %>%
  group_by(migcertainty, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_migcertainty <- left_join(England_2015_2020_aggregate_weekly_conscounts_migcertainty, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_migcertainty$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_migcertainty <- England_2015_2020_aggregate_weekly_conscounts_migcertainty %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_migcertainty <- England_2015_2020_aggregate_weekly_conscounts_migcertainty %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_migcertainty, file = "filepath")

rm(list=ls())
gc()

# 11_Create aggregate weekly data (London only, migrant status)  ------

load(file = "filepath")

London_2015_2020_aggregate_weekly_conscounts <- cohort_England_2015_2020_weekly_conscounts %>%
  filter(prac_region == "London") 

London_2015_2020_aggregate_weekly_conscounts <- London_2015_2020_aggregate_weekly_conscounts %>%
  group_by(migrant_status, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

London_2015_2020_aggregate_weekly_conscounts <- left_join(London_2015_2020_aggregate_weekly_conscounts, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
London_2015_2020_aggregate_weekly_conscounts$lockdown1 <- 0
London_2015_2020_aggregate_weekly_conscounts <- London_2015_2020_aggregate_weekly_conscounts %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

London_2015_2020_aggregate_weekly_conscounts <- London_2015_2020_aggregate_weekly_conscounts %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(London_2015_2020_aggregate_weekly_conscounts, file = "filepath")

rm(list=ls())
gc()

# 12_Create aggregate weekly data (certainty of migration by age_subcohort)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort <- cohort_England_2015_2020_weekly_conscounts %>%
  group_by(migcertainty, age_subcohort, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort <- left_join(England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort, file = "filepath")

rm(list=ls())
gc()

## ITS SENSITIVITY ANALYSIS 1 - alternative cohort matched on data_start variables --------------------

# 13_Join weekly conscount variables onto cohort weekly records file  ------

load(file = "filepath")

# Load cohort file and join
load(file = "filepath")

cohort_England_2015_2020_weekly_conscounts_ds <- matched_its_weekly_records_ds %>%
  left_join(cons_2015_2020_weekly_conscounts, by = c("patid" = "patid", "studyweek" = "week"))

rm(matched_its_weekly_records_ds,cons_2015_2020_weekly_conscounts)

# Replace NAs in cons count variables with 0
cohort_England_2015_2020_weekly_conscounts_ds$conscount <- cohort_England_2015_2020_weekly_conscounts_ds$conscount %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts_ds$phone <- cohort_England_2015_2020_weekly_conscounts_ds$phone %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts_ds$facetoface <- cohort_England_2015_2020_weekly_conscounts_ds$facetoface %>% replace_na(0)
sum(is.na(cohort_England_2015_2020_weekly_conscounts_ds$conscount))
glimpse(cohort_England_2015_2020_weekly_conscounts_ds)

# Save data
save(cohort_England_2015_2020_weekly_conscounts_ds, file = "filepath")

rm(list = ls())
gc()

# 14_Create aggregate weekly data (migrant status)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_ds <- cohort_England_2015_2020_weekly_conscounts_ds %>%
  group_by(migrant_status, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            pweeks=sum(pweeks),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Check that no pweeks are greater than popsize
England_2015_2020_aggregate_weekly_conscounts_ds$pweeks <= England_2015_2020_aggregate_weekly_conscounts_ds$popsize

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_ds <- left_join(England_2015_2020_aggregate_weekly_conscounts_ds, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_ds$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_ds <- England_2015_2020_aggregate_weekly_conscounts_ds %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_ds <- England_2015_2020_aggregate_weekly_conscounts_ds %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_ds, file = "filepath")

rm(list=ls())
gc()

# 15_Not included in paper ---
# 15a_Create aggregate weekly data (migrant status age subcohort)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds <- cohort_England_2015_2020_weekly_conscounts_ds %>%
  group_by(migrant_status, studyweek, age_subcohort) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds <- left_join(England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds <- England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds <- England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_agesubcohort_ds, file = "filepath")

rm(list=ls())
gc()

# 15b_Create aggregate weekly data (ethnicity)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds <- cohort_England_2015_2020_weekly_conscounts_ds %>%
  group_by(ethnicat6, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds <- left_join(England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds <- England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds <- England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_ethnicity_ds, file = "filepath")

rm(list=ls())
gc()

# 15c_Create aggregate weekly data (migrant_status and ethnicity)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds <- cohort_England_2015_2020_weekly_conscounts_ds %>%
  group_by(migrant_status, studyweek, ethnicat6) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds <- left_join(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ds, file = "filepath")

rm(list=ls())
gc()

# 15d_Create aggregate weekly data (certainty of migration)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds <- cohort_England_2015_2020_weekly_conscounts_ds %>%
  group_by(migcertainty, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds <- left_join(England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_migcertainty_ds, file = "filepath")

rm(list=ls())
gc()

# 15f_Create aggregate weekly data (London only, migrant status)  ------

load(file = "filepath")

London_2015_2020_aggregate_weekly_conscounts_ds <- cohort_England_2015_2020_weekly_conscounts_ds %>%
  filter(prac_region == "London") 

London_2015_2020_aggregate_weekly_conscounts_ds <- London_2015_2020_aggregate_weekly_conscounts_ds %>%
  group_by(migrant_status, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

London_2015_2020_aggregate_weekly_conscounts_ds <- left_join(London_2015_2020_aggregate_weekly_conscounts_ds, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
London_2015_2020_aggregate_weekly_conscounts_ds$lockdown1 <- 0
London_2015_2020_aggregate_weekly_conscounts_ds <- London_2015_2020_aggregate_weekly_conscounts_ds %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

London_2015_2020_aggregate_weekly_conscounts_ds <- London_2015_2020_aggregate_weekly_conscounts_ds %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(London_2015_2020_aggregate_weekly_conscounts_ds, file = "filepath")

rm(list=ls())
gc()

# 15g_Create aggregate weekly data (certainty of migration by age_subcohort)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds <- cohort_England_2015_2020_weekly_conscounts_ds %>%
  group_by(migcertainty, age_subcohort, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort <- left_join(England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort_ds, file = "filepath")

rm(list=ls())
gc()
## ITS SENSITIVITY ANALYSIS 2 - alternative cohort matched on age and gender only --------------------

# 16_Join weekly conscount variables onto cohort weekly records file  ------

load(file = "filepath")

# Load cohort file and join
load(file = "filepath")

cohort_England_2015_2020_weekly_conscounts_SA2 <- matched_its_weekly_records_SA2 %>%
  left_join(cons_2015_2020_weekly_conscounts, by = c("patid" = "patid", "studyweek" = "week"))

rm(matched_its_weekly_records_ds,cons_2015_2020_weekly_conscounts)

# Replace NAs in cons count variables with 0
cohort_England_2015_2020_weekly_conscounts_SA2$conscount <- cohort_England_2015_2020_weekly_conscounts_SA2$conscount %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts_SA2$phone <- cohort_England_2015_2020_weekly_conscounts_SA2$phone %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts_SA2$facetoface <- cohort_England_2015_2020_weekly_conscounts_SA2$facetoface %>% replace_na(0)
sum(is.na(cohort_England_2015_2020_weekly_conscounts_SA2$conscount))
glimpse(cohort_England_2015_2020_weekly_conscounts_SA2)

# Save data
save(cohort_England_2015_2020_weekly_conscounts_SA2, file = "filepath")

rm(list = ls())
gc()

# 17_Create aggregate weekly data (migrant status)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_SA2 <- cohort_England_2015_2020_weekly_conscounts_SA2 %>%
  group_by(migrant_status, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            pweeks=sum(pweeks),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Check that no pweeks are greater than popsize
England_2015_2020_aggregate_weekly_conscounts_SA2$pweeks <= England_2015_2020_aggregate_weekly_conscounts_SA2$popsize

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_SA2 <- left_join(England_2015_2020_aggregate_weekly_conscounts_SA2, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_SA2$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_SA2 <- England_2015_2020_aggregate_weekly_conscounts_SA2 %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_SA2 <- England_2015_2020_aggregate_weekly_conscounts_SA2 %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_SA2, file = "filepath")

rm(list=ls())
gc()

# 15_Not included in paper ---
## ITS SENSITIVITY ANALYSIS 3 - alternative cohort matched on age, prac_region and gender only --------------------

# 18_Join weekly conscount variables onto cohort weekly records file  ------

load(file = "filepath")

# Load cohort file and join
load(file = "filepath")

cohort_England_2015_2020_weekly_conscounts_SA3 <- matched_its_weekly_records_SA3 %>%
  left_join(cons_2015_2020_weekly_conscounts, by = c("patid" = "patid", "studyweek" = "week"))

rm(matched_its_weekly_records_SA3,cons_2015_2020_weekly_conscounts)

# Replace NAs in cons count variables with 0
cohort_England_2015_2020_weekly_conscounts_SA3$conscount <- cohort_England_2015_2020_weekly_conscounts_SA3$conscount %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts_SA3$phone <- cohort_England_2015_2020_weekly_conscounts_SA3$phone %>% replace_na(0)
cohort_England_2015_2020_weekly_conscounts_SA3$facetoface <- cohort_England_2015_2020_weekly_conscounts_SA3$facetoface %>% replace_na(0)
sum(is.na(cohort_England_2015_2020_weekly_conscounts_SA3$conscount))
glimpse(cohort_England_2015_2020_weekly_conscounts_SA3)

# Save data
save(cohort_England_2015_2020_weekly_conscounts_SA3, file = "filepath")

rm(list = ls())
gc()

# 19_Create aggregate weekly data (migrant status)  ------

load(file = "filepath")

England_2015_2020_aggregate_weekly_conscounts_SA3 <- cohort_England_2015_2020_weekly_conscounts_SA3 %>%
  group_by(migrant_status, studyweek) %>%
  summarise(conscount=sum(conscount),
            phone=sum(phone),
            facetoface=sum(facetoface),
            pdays=sum(pdays),
            pyears=sum(pyears),
            pweeks=sum(pweeks),
            popsize=sum(n_distinct(patid)),
            studymonth=first(studymonth))

# Check that no pweeks are greater than popsize
England_2015_2020_aggregate_weekly_conscounts_SA3$pweeks <= England_2015_2020_aggregate_weekly_conscounts_SA3$popsize

# Add study year
load(file = "filepath")
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>%
  mutate(studyyear = year(date))
date_to_week_conversion_2015_2020 <- dplyr::select(date_to_week_conversion_2015_2020, -c(date))
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020 %>% distinct()
date_to_week_conversion_2015_2020 <- date_to_week_conversion_2015_2020[-c(54, 160, 213, 266),]# only take the year at the start of the weeks

England_2015_2020_aggregate_weekly_conscounts_SA3 <- left_join(England_2015_2020_aggregate_weekly_conscounts_SA3, date_to_week_conversion_2015_2020, by=c("studyweek" = "week"))

# Create binary lockdown 1 variable
England_2015_2020_aggregate_weekly_conscounts_SA3$lockdown1 <- 0
England_2015_2020_aggregate_weekly_conscounts_SA3 <- England_2015_2020_aggregate_weekly_conscounts_SA3 %>%
  mutate(lockdown1 = replace(lockdown1, studyweek >= 275, 1))

# Add a date variable for plotting purposes 
load(file = "filepath")
first_day_week <- date_to_week_conversion_2015_2020 %>%
  group_by(week) %>%
  filter(row_number()==1)

England_2015_2020_aggregate_weekly_conscounts_SA3 <- England_2015_2020_aggregate_weekly_conscounts_SA3 %>%
  left_join(first_day_week, by = c('studyweek' = 'week'))

# Save
save(England_2015_2020_aggregate_weekly_conscounts_SA3, file = "filepath")

rm(list=ls())
gc()
