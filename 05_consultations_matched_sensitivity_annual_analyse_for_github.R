## 0_Description ---------------------------------------------------------------------------

# Migrants' primary care utilisation before and during the COVID-19 pandemic in England: An interrupted time series
# Annual pre-pandemic sensitivity analysis 
# Date started: 25/03/2021
# Author(s): Claire Zhang / Yamina Boukari
# QA (date): Yamina Boukari (09/03/2022)

## 1_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc, epitools, data.table, epiR, MASS, psych, forestplot, dsr, gmodels, Epi, msm)

## 2_Set working directory ------------------------------------------------------------------

setwd("filepath")

## 3_Functions ------------------------------------------------------

function_to_calculate_IR <- function(dataset, output_name) {
  output <- pois.exact(dataset$conscount, dataset$pyears, conf.level=0.95)
  output <- output %>%
    mutate(incidence_rate1=rate*1) # format as per 1 person-year (change number if you want per 10,000 person-years etc)
  output <- output %>%
    mutate(lower1=lower*1)
  output <- output %>%
    mutate(upper1=upper*1)
  output <- output %>% 
    rename(incidence_rate = rate)
  output <- output[, 3:9]
  output <- bind_cols(dataset, output) 
  output <- output %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower1, output$upper1, sep ="-")
  output$ir_ci <- paste(output$incidence_rate1, output$ci, sep =",")
  output <-  dplyr::select(output, -c(incidence_rate, lower, upper, conf.level, incidence_rate1, lower1,upper1,ci))
  output <- output %>% 
    rename(events = conscount) %>%
    rename(person_years = pyears)
  assign(x = output_name, value = output, envir = globalenv())
} # per 1 person-year

extract_glm_results_allages <- function (x, output_name, dataset, variable) {
  names <- names(coef(x))
  names <- sub(".*)", "", names) # extracts all text after ')' from names to avoid having to recode later (excluding the intercept, the below lines of code are deal with this)
  if (is.factor(eval(substitute(variable), dataset)) == F) {
    variable <- as.factor(eval(substitute(variable), dataset))
  } # if a variable isn't a factor (which is needed for the below line), e.g. studyyear, it converts it to a factor
  names[1] <- levels(eval(substitute(variable), dataset))[1] # takes the reference level of a specified factor variable and replaces it with 'Intercept' in names
  estimate<- exp(coef(x))
  confint.default <- exp(confint.default(x))
  p <- coef(summary(x))[,4]
  output <- cbind(names,estimate,confint.default,p)
  output <- as_tibble(output)
  output$estimate <- as.numeric(output$estimate)
  output <- output %>% rename(lower = "2.5 %") %>% rename(upper = "97.5 %") 
  output$upper <- as.numeric(output$upper)
  output$lower <- as.numeric(output$lower)
  output$p <- as.numeric(output$p)
  output <- output %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$estimate[1] <- 1.00
  output$lower[1] <- 1.00
  output$upper[1] <- 1.00
  output$ci[1] <- "1.00-1.00"
  output$irr_ci <- paste(output$estimate, output$ci, sep =",")
  output_table <- dplyr::select(output,names, irr_ci )
  output_table <- output_table %>%
    mutate(age_subcohort = 'All_ages')
  assign(x = output_name, value = output, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
  assign(x = paste0(output_name,'_table'), value = output_table, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
}

extract_glm_results_agesubcohorts <- function (x, output_name, subcohort, dataset, variable) {
  names <- names(coef(x))
  names <- sub(".*)", "", names) # extracts all text after ')' from names to avoid having to recode later (excluding the intercept, the below lines of code are deal with this)
  if (is.factor(eval(substitute(variable), dataset)) == F) {
    variable <- as.factor(eval(substitute(variable), dataset))
  } # if a variable isn't a factor (which is needed for the below line), e.g. studyyear, it converts it to a factor)
  names[1] <- levels(eval(substitute(variable), dataset))[1] # takes the reference level of a specified factor variable and replaces it with 'Intercept' in names
  estimate<- exp(coef(x))
  confint.default <- exp(confint.default(x))
  p <- coef(summary(x))[,4]
  output <- cbind(names,estimate,confint.default,p)
  output <- as_tibble(output)
  output$estimate <- as.numeric(output$estimate)
  output <- output %>% rename(lower = "2.5 %") %>% rename(upper = "97.5 %") 
  output$upper <- as.numeric(output$upper)
  output$lower <- as.numeric(output$lower)
  output$p <- as.numeric(output$p)
  output <- output %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$estimate[1] <- 1.00
  output$lower[1] <- 1.00
  output$upper[1] <- 1.00
  output$ci[1] <- "1.00-1.00"
  output$irr_ci <- paste(output$estimate, output$ci, sep =",")
  output <- output %>%
    mutate(age_subcohort = subcohort)
  output_table <- dplyr::select(output, names, irr_ci, age_subcohort)
  assign(x = output_name, value = output, envir = globalenv()) 
  assign(x = paste0(output_name,'_table'), value = output_table, envir = globalenv()) 
}

# 4_Sensitivity analysis 1 - MATCHED ON: age_data_start, year_data_start and prac_region -------

## Load dataset and select variables ----

load(file = "filepath")

# Matched on age_data_start, year_data_start and prac_region 

matched_cohort_England_2015_2020_annual_conscounts_test_ds <- dplyr::select(matched_cohort_England_2015_2020_annual_conscounts_test_ds, 
                                                                            c(patid, pyears, conscount, migrant_status, migcertainty,
                                                                              age_subcohort, studyyear, prac_region, gender,
                                                                              ethnicat6, imd, studyyear_agecat))

# Relevel
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, gender <- relevel (gender, ref="Male")) 
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, imd <- relevel (imd, ref="IMD 1"))
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, prac_region <- relevel (prac_region, ref="London")) 
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, ethnicat6 <- relevel (ethnicat6, ref="White British"))


## Summary outcome measures (not included in paper) ------------------------------------------------------------------

### Overall (i.e. full period)

## All patients

# All patients
conscount_summary_overall_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% 
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migrant status
conscount_summary_mig_status_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migcertainty
conscount_summary_migcertainty_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(migcertainty) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Join for overall all patients + by migrant status
conscount_summary_overall_1 <- full_join(conscount_summary_overall_allpatients, conscount_summary_mig_status_allpatients, 
                                         by = c("n_events"="n_events", "mean" = "mean",
                                                "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals")) %>%
  full_join(conscount_summary_migcertainty_allpatients, by = c("n_events"="n_events", "mean" = "mean",
                                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                                                               "migrant_status" = "migcertainty","n_individuals"="n_individuals"))

conscount_summary_overall_1$migrant_status <- fct_explicit_na(conscount_summary_overall_1$migrant_status, na_level = "All_patients")
conscount_summary_overall_1$age <- 'all_ages'
conscount_summary_overall_1 <- conscount_summary_overall_1 %>% relocate(age, migrant_status)

## Age subcohorts

# Age_subcohorts
conscount_summary_overall_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(age_subcohort) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migrant status
conscount_summary_mig_status_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(age_subcohort, migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migcertainty
conscount_summary_migcertainty_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(age_subcohort,migcertainty) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Join for overall age_subcohorts + by migrant status
conscount_summary_overall_2 <- full_join(conscount_summary_overall_agesubcohorts, conscount_summary_mig_status_agesubcohorts, 
                                         by = c('age_subcohort' = 'age_subcohort', "n_events"="n_events", "mean" = "mean",
                                                "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals")) %>%
  full_join(conscount_summary_migcertainty_agesubcohorts, by = c('age_subcohort' = 'age_subcohort',"migrant_status" = "migcertainty","n_events"="n_events", "mean" = "mean",
                                                                 "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals"))
conscount_summary_overall_2$migrant_status <- fct_explicit_na(conscount_summary_overall_2$migrant_status, na_level = "All_patients")
conscount_summary_overall_2 <- conscount_summary_overall_2 %>% relocate(age_subcohort, migrant_status) %>%
  rename(age = age_subcohort) %>%
  arrange(age)

# Join all patients to the age_subcohorts
matched_conscount_summary_fullperiod_ds <- bind_rows(conscount_summary_overall_1,conscount_summary_overall_2)

# save file
write_csv(matched_conscount_summary_fullperiod_ds, "filepath" )

### Annual

## All patients

# All patients
conscount_annual_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migrant status
conscount_annual_mig_status_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migcertainty 
conscount_annual_migcertainty_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(migcertainty, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))


# Join to give annual for allpatients + migrant status
conscount_summary_annual_1 <- full_join(conscount_annual_allpatients, conscount_annual_mig_status_allpatients, 
                                        by = c("studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max")) %>%
  full_join(conscount_annual_migcertainty_allpatients,
            by = c("studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                   "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                   "migrant_status" = "migcertainty"))

conscount_summary_annual_1$migrant_status <- fct_explicit_na(conscount_summary_annual_1$migrant_status, na_level = "All_patients") 
conscount_summary_annual_1$age <- 'all_ages'
conscount_summary_annual_1 <- conscount_summary_annual_1 %>% relocate(age, migrant_status)

## age_subcohorts

# age_subcohorts
conscount_annual_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(age_subcohort, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migrant status
conscount_annual_mig_status_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(age_subcohort,migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migcertainty 
conscount_annual_migcertainty_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% group_by(age_subcohort, migcertainty, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))


## Join to give annual for age_subcohorts + migrant status
conscount_summary_annual_2 <- full_join(conscount_annual_agesubcohorts, conscount_annual_mig_status_agesubcohorts, 
                                        by = c('age_subcohort' = 'age_subcohort', "studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max")) %>%
  full_join(conscount_annual_migcertainty_agesubcohorts,
            by = c('age_subcohort' = 'age_subcohort',"studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                   "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                   "migrant_status" = "migcertainty"))

conscount_summary_annual_2$migrant_status <- fct_explicit_na(conscount_summary_annual_2$migrant_status, na_level = "All_patients") 
conscount_summary_annual_2 <- conscount_summary_annual_2 %>% relocate(age_subcohort, migrant_status) %>%
  rename(age = age_subcohort) %>%
  arrange(age)

# Join all patients to the age_subcohorts
matched_conscount_summary_annual_ds <- bind_rows(conscount_summary_annual_1,conscount_summary_annual_2)

## Export  summary measures
write_csv(matched_conscount_summary_annual_ds, "filepath" )

## IRs by migrant_status + studyyear (not included in paper) -----

# IRs by migrant_status and studyyear

## All ages

pyears_allpatients_overall_migstatus_studyyear <- aggregate(pyears ~ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_allpatients_overall_migstatus_studyyear <- pyears_allpatients_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus_studyyear <- aggregate(conscount ~ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_allpatients_overall_migstatus_studyyear <- inner_join(conscount_allpatients_overall_migstatus_studyyear,pyears_allpatients_overall_migstatus_studyyear, by= c("migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migstatus_studyyear, output_name = 'IR_allpatients_overall_migstatus_studyyear')
IR_allpatients_overall_migstatus_studyyear$age_subcohort <- 'All_ages'

## age_subcohorts

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus_studyyear <- aggregate(pyears ~ age_subcohort + migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_agesubcohorts_overall_migstatus_studyyear <- pyears_agesubcohorts_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus_studyyear <- aggregate(conscount ~ age_subcohort+ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_agesubcohorts_overall_migstatus_studyyear <- inner_join(conscount_agesubcohorts_overall_migstatus_studyyear,pyears_agesubcohorts_overall_migstatus_studyyear, 
                                                           by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus_studyyear, output_name = 'IR_agesubcohorts_overall_migstatus_studyyear')

## Combine

IR_migstatus_studyyear <- full_join(IR_allpatients_overall_migstatus_studyyear, IR_agesubcohorts_overall_migstatus_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migrant_status' = 'migrant_status',
                                                                                                                                     'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))

## Reorder results (can change depending on how we want to present it)
IR_migstatus_studyyear <- arrange(IR_migstatus_studyyear, migrant_status, age_subcohort) %>%
  relocate(age_subcohort)

# IR_migstatus_studyyear <- arrange(IR_migstatus_studyyear, studyyear) %>%
#   relocate(age_subcohort)

write_csv(IR_migstatus_studyyear, "filepath")

## IR & univariable IRR by migrant status ----------------------------------------------------------------------

# Filter out 2020 for modelling
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>%
  filter(studyyear != 2020)

## All patients (i.e. all ages)

# IR by migrant status
pyears_allpatients_overall_migstatus <- aggregate(pyears ~ migrant_status, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_allpatients_overall_migstatus <- pyears_allpatients_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus <- aggregate(conscount ~ migrant_status, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_allpatients_overall_migstatus <- inner_join(conscount_allpatients_overall_migstatus,pyears_allpatients_overall_migstatus, by= c("migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_allpatients_overall_migstatus, output_name = 'IR_allpatients_overall_migstatus')
IR_allpatients_overall_migstatus$age_subcohort <- 'All_ages'

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = matched_cohort_England_2015_2020_annual_conscounts_test_ds)

extract_glm_results_allages(x, 'glm_mig', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)
round(ci.lin(x,Exp=T),3)

save(glm_mig, file="filepath")

# Join glm + IR + univariable_mig 
univariable_migrant_status <- full_join(IR_allpatients_overall_migstatus, glm_mig_table, by = c("age_subcohort" = "age_subcohort","migrant_status" = "names")) %>%
  relocate(age_subcohort)


## Age_subcohorts (not included in paper)

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus <- aggregate(pyears ~ age_subcohort + migrant_status, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_agesubcohorts_overall_migstatus <- pyears_agesubcohorts_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus <- aggregate(conscount ~ age_subcohort+ migrant_status, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_agesubcohorts_overall_migstatus <- inner_join(conscount_agesubcohorts_overall_migstatus,pyears_agesubcohorts_overall_migstatus, 
                                                 by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus, output_name = 'IR_agesubcohorts_overall_migstatus')

# IRR by migrant status (univariable glm)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)
round(ci.lin(x,Exp=T),3)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_50to64", "50-64 years",matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 65 years and over  
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## NOTE IRR is comparing migrants at X age to non-migrants of the SAME age 

# Join age_subcohorts glm+IR and all_ages glm+IR 
glm_mig_combined_table <- bind_rows(glm_mig_0to15_table,glm_mig_16to24_table,glm_mig_25to34_table,glm_mig_35to49_table,glm_mig_50to64_table,glm_mig_65plus_table)

glm_mig_combined_table$names <- as.factor(glm_mig_combined_table$names)

univariable_migrant_status_agesubcohorts <- left_join(IR_agesubcohorts_overall_migstatus, glm_mig_combined_table, 
                                                      by = c("migrant_status" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migrant_status <- bind_rows(univariable_migrant_status, univariable_migrant_status_agesubcohorts)  
univariable_migrant_status$age_subcohort <- as.factor(univariable_migrant_status$age_subcohort)
univariable_migrant_status$age_subcohort <- factor(univariable_migrant_status$age_subcohort, 
                                                   levels = c('All_ages', '0-15 years','16-24 years', '25-34 years', '35-49 years','50-64 years', '>=65 years'))
## Reorder results (can change depending on how we want to present it)
univariable_migrant_status_ds <- arrange(univariable_migrant_status, migrant_status, age_subcohort)

write_csv(univariable_migrant_status_ds, "filepath")

## Multivariable model by migrant_status with IMD adjustment  ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)),
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds)
extract_glm_results_allages(x, 'multivariable_matched_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)
round(ci.lin(x,Exp=T),3)

## age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)
round(ci.lin(x,Exp=T),3)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)
round(ci.lin(x,Exp=T),3)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)
round(ci.lin(x,Exp=T),3)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)
round(ci.lin(x,Exp=T),5)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)
round(ci.lin(x,Exp=T),3)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)
round(ci.lin(x,Exp=T),3)

# Join all ages and age_subcohorts glm results
multivariable_matched_migrant_status_ds <- bind_rows(multivariable_matched_allages_table, glm_multivariable_matched_0to15_table,glm_multivariable_matched_16to24_table,
                                                  glm_multivariable_matched_25to34_table,glm_multivariable_matched_35to49_table,glm_multivariable_matched_50to64_table,glm_multivariable_matched_65plus_table)

write_csv(multivariable_matched_migrant_status_ds, "filepath")

# Save .Rdata for forestplots
save(multivariable_matched_allages, file="filepath")
save(glm_multivariable_matched_0to15, file="filepath")
save(glm_multivariable_matched_16to24, file="filepath")
save(glm_multivariable_matched_25to34, file="filepath")
save(glm_multivariable_matched_35to49, file="filepath")
save(glm_multivariable_matched_50to64, file="filepath")
save(glm_multivariable_matched_65plus, file="filepath")

## Multivariable model by migrant_status without IMD adjustment  ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)),
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds)
extract_glm_results_allages(x, 'multivariable_matched_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)
round(ci.lin(x,Exp=T),3)

## age_subcohorts (not included in paper)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_matched_migrant_status_noIMD_ds <- bind_rows(multivariable_matched_allages_table, glm_multivariable_matched_0to15_table,glm_multivariable_matched_16to24_table,
                                                  glm_multivariable_matched_25to34_table,glm_multivariable_matched_35to49_table,glm_multivariable_matched_50to64_table,glm_multivariable_matched_65plus_table)

write_csv(multivariable_matched_migrant_status_noIMD_ds, "filepath")

# Save .Rdata for forestplots
save(multivariable_matched_allages, file="filepath")
save(glm_multivariable_matched_0to15, file="filepath")
save(glm_multivariable_matched_16to24, file="filepath")
save(glm_multivariable_matched_25to34, file="filepath")
save(glm_multivariable_matched_35to49, file="filepath")
save(glm_multivariable_matched_50to64, file="filepath")
save(glm_multivariable_matched_65plus, file="filepath")

## Analyses not included in paper ------
### Certainty of migration status  ----

# Reload to include 2020
load(file = "filepath")

matched_cohort_England_2015_2020_annual_conscounts_test_ds <- dplyr::select(matched_cohort_England_2015_2020_annual_conscounts_test_ds, 
                                                                    c(patid, pyears, conscount, migrant_status,
                                                                      age_subcohort, imd, ethnicat6, gender,
                                                                      studyyear, studyyear_agecat, prac_region,
                                                                      migcertainty, cohort_entry))

# Relevel
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, gender <- relevel (gender, ref="Male")) 
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, imd <- relevel (imd, ref="IMD 1"))
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, prac_region <- relevel (prac_region, ref="London")) 
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- within(matched_cohort_England_2015_2020_annual_conscounts_test_ds, ethnicat6 <- relevel (ethnicat6, ref="White British"))

# IRs by migcertainty + studyyear

## All ages

pyears_allpatients_overall_migcertainty_studyyear <- aggregate(pyears ~ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_allpatients_overall_migcertainty_studyyear <- pyears_allpatients_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty_studyyear <- aggregate(conscount ~ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_allpatients_overall_migcertainty_studyyear <- inner_join(conscount_allpatients_overall_migcertainty_studyyear,pyears_allpatients_overall_migcertainty_studyyear, by= c("migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migcertainty_studyyear, output_name = 'IR_allpatients_overall_migcertainty_studyyear')
IR_allpatients_overall_migcertainty_studyyear$age_subcohort <- 'All_ages'

## Age_subcohorts

# IR by migcertainty
pyears_agesubcohorts_overall_migcertainty_studyyear <- aggregate(pyears ~ age_subcohort + migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_agesubcohorts_overall_migcertainty_studyyear <- pyears_agesubcohorts_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty_studyyear <- aggregate(conscount ~ age_subcohort+ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_agesubcohorts_overall_migcertainty_studyyear <- inner_join(conscount_agesubcohorts_overall_migcertainty_studyyear,pyears_agesubcohorts_overall_migcertainty_studyyear, 
                                                              by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty_studyyear, output_name = 'IR_agesubcohorts_overall_migcertainty_studyyear')

## Combine

IR_migcertainty_studyyear <- full_join(IR_allpatients_overall_migcertainty_studyyear, IR_agesubcohorts_overall_migcertainty_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migcertainty' = 'migcertainty',
                                                                                                                                              'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))
## Reorder results (can change depending on how we want to present it)
IR_migcertainty_studyyear_ds <- arrange(IR_migcertainty_studyyear, migcertainty, age_subcohort) %>%
  relocate(age_subcohort)

write_csv(IR_migcertainty_studyyear_ds, "results")

# Filter out 2020 for modelling
matched_cohort_England_2015_2020_annual_conscounts_test_ds <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>%
  filter(studyyear != 2020)

## Univariable model

## All patients (i.e. all ages)

# IR 
pyears_allpatients_overall_migcertainty <- aggregate(pyears ~ migcertainty, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_allpatients_overall_migcertainty <- pyears_allpatients_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty <- aggregate(conscount ~ migcertainty, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_allpatients_overall_migcertainty <- inner_join(conscount_allpatients_overall_migcertainty,pyears_allpatients_overall_migcertainty, by= c("migcertainty" = "migcertainty"))

function_to_calculate_IR(IR_allpatients_overall_migcertainty, output_name = 'IR_allpatients_overall_migcertainty')
IR_allpatients_overall_migcertainty$age_subcohort <- 'All_ages'

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), data = matched_cohort_England_2015_2020_annual_conscounts_test_ds)
extract_glm_results_allages(x, 'glm_migcertainty', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

save(glm_migcertainty, file="filepath")

# Join glm + IR + univariable_mig 
univariable_migcertainty <- full_join(IR_allpatients_overall_migcertainty, glm_migcertainty_table, by = c("age_subcohort" = "age_subcohort","migcertainty" = "names")) %>%
  relocate(age_subcohort)

## Age_subcohorts

# IR
pyears_agesubcohorts_overall_migcertainty <- aggregate(pyears ~ age_subcohort + migcertainty, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
pyears_agesubcohorts_overall_migcertainty <- pyears_agesubcohorts_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty <- aggregate(conscount ~ age_subcohort+ migcertainty, matched_cohort_England_2015_2020_annual_conscounts_test_ds, sum) 
IR_agesubcohorts_overall_migcertainty <- inner_join(conscount_agesubcohorts_overall_migcertainty,pyears_agesubcohorts_overall_migcertainty, 
                                                    by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty"))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty, output_name = 'IR_agesubcohorts_overall_migcertainty')

# IRR (univariable glm)
### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

## 65 years and over NOTE  - haven't run this yet as the sample dataset doesn't have any >65s in
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)


# Join age_subcohorts glm+IR and all_ages glm+IR
glm_migcertainty_combined_table <- bind_rows(glm_migcertainty_0to15_table,glm_migcertainty_16to24_table,glm_migcertainty_25to34_table,glm_migcertainty_35to49_table,glm_migcertainty_50to64_table,glm_migcertainty_65plus_table)

univariable_migcertainty_agesubcohorts <- left_join(IR_agesubcohorts_overall_migcertainty, glm_migcertainty_combined_table, 
                                                    by = c("migcertainty" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migcertainty <- bind_rows(univariable_migcertainty, univariable_migcertainty_agesubcohorts)  
univariable_migcertainty$age_subcohort <- as.factor(univariable_migcertainty$age_subcohort)
univariable_migcertainty$age_subcohort <- factor(univariable_migcertainty$age_subcohort, 
                                                 levels = c('All_ages', '0-15 years','16-24 years', '25-34 years', '35-49 years','50-64 years', '>=65 years'))
univariable_migcertainty_ds <- arrange(univariable_migcertainty, age_subcohort)


write_csv(univariable_migcertainty_ds, "results")

## Multivariable model controlling for year + gender + as.factor(studyyear_agecat) + as.factor(imd) + prac_region

## All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds)
extract_glm_results_allages(x, 'multivariable_migcertainty_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds[matched_cohort_England_2015_2020_annual_conscounts_test_ds$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migcertainty)

# Join all ages and age_subcohorts glm results
multivariable_migcertainty_ds <- bind_rows(multivariable_migcertainty_allages_table, glm_multivariable_migcertainty_0to15_table,glm_multivariable_migcertainty_16to24_table,
                                        glm_multivariable_migcertainty_25to34_table,glm_multivariable_migcertainty_35to49_table,glm_multivariable_migcertainty_50to64_table,glm_multivariable_migcertainty_65plus_table)


write_csv(multivariable_migcertainty_ds, "filepath")

# Save .Rdata for forestplots
save(multivariable_migcertainty_allages, file="filepath")
save(glm_multivariable_migcertainty_0to15, file="filepath")
save(glm_multivariable_migcertainty_16to24, file="filepath")
save(glm_multivariable_migcertainty_25to34, file="filepath")
save(glm_multivariable_migcertainty_35to49, file="filepath")
save(glm_multivariable_migcertainty_50to64, file="filepath")
save(glm_multivariable_migcertainty_65plus, file="filepath")


### Ethnicity (ethnicat6) - interaction term, additive and multiplicative effects ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status)*as.factor(ethnicat6) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = matched_cohort_England_2015_2020_annual_conscounts_test_ds)
extract_glm_results_allages(x, 'multivariable_ethnicat6interaction_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# Save .Rdata for forestplots
save(multivariable_ethnicat6interaction_allages, file="filepath")
write_csv(multivariable_ethnicat6interaction_allages, "filepath")

# RR (95% CI) for non-migrants in each ethnicity compared to White British non-migrants
# For table
white_non_migrant <- data.frame(Ethnicity = 'White British', RR_nonmigrantsVsWBNM = '1.0')
non_migrants_by_ethnicity <- multivariable_ethnicat6interaction_allages[3:8,] %>%
  mutate(ci = paste(lower, upper, sep ="-"),
         RR_nonmigrantsVsWBNM = paste0(estimate, ' (', ci, ")")) %>%
  mutate(Ethnicity = c('White non-British', 'Mixed', 'Asian', 'Black', 'Other', 'Unknown')) %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, estimate, upper, lower)
non_migrants_by_ethnicity <- bind_rows(white_non_migrant, non_migrants_by_ethnicity)

# Calculate interaction effects - migrants in each ethnicity compared to White British non-migrant (WBNM) reference group

get_RR_and_95_CI <- function(dataframe) {
  output <- dataframe %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$Lower.CI, output$Upper.CI, sep ="-")
  output$RR <- paste0(output$Estimate, ' (', output$ci, ")")
  #output <- output %>%
  #  dplyr:: select(RR)
  output_name <- deparse(substitute(dataframe))
  assign(x = output_name, value = output, envir = globalenv())
}

## White British migrants 
white_migrant_vsWBNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_migrant_vsWBNM)
save(white_migrant_vsWBNM, file="filepath")
white_migrant <- white_migrant_vsWBNM %>%
  dplyr::select(RR)

# White non-British migrants 
white_nonbritish_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)White non-British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_migrant_vsWBNM)
save(white_nonbritish_migrant_vsWBNM, file="filepath")
white_nonbritish_migrant <- white_nonbritish_migrant_vsWBNM %>%
  dplyr::select(RR)

## Black migrants 
black_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_migrant_vsWBNM)
save(black_migrant_vsWBNM, file="filepath")
black_migrant <- black_migrant_vsWBNM %>%
  dplyr::select(RR)

## Asian
asian_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Asian/Asian British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_migrant_vsWBNM)
save(asian_migrant_vsWBNM, file="filepath")
asian_migrant <- asian_migrant_vsWBNM %>%
  dplyr::select(RR)

## Mixed
mixed_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_migrant_vsWBNM)
save(mixed_migrant_vsWBNM, file="filepath")
mixed_migrant <- mixed_migrant_vsWBNM %>%
  dplyr::select(RR)

## Other
other_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Other ethnic group' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_migrant_vsWBNM)
save(other_migrant_vsWBNM, file="filepath")
other_migrant <- other_migrant_vsWBNM %>%
  dplyr::select(RR)

## Unknown
unknown_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Unknown' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_migrant_vsWBNM)
save(unknown_migrant_vsWBNM, file="filepath")
unknown_migrant <- unknown_migrant_vsWBNM %>%
  dplyr::select(RR)

# Calculate interaction effects - migrants vs non-migrants in each ethnicity strata
## White
white_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_strata_vsNM)
save(white_strata_vsNM, file="filepath")
white_strata <- white_strata_vsNM %>%
  dplyr::select(RR)

## White non-British
white_nonbritish_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_strata_vsNM)
save(white_nonbritish_strata_vsNM, file="filepath")
white_nonbritish_strata <- white_nonbritish_strata_vsNM %>%
  dplyr::select(RR)

## Black
black_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_strata_vsNM)
save(black_strata_vsNM, file="filepath")
black_strata <- black_strata_vsNM %>%
  dplyr::select(RR)

## Asian
asian_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_strata_vsNM)
save(asian_strata_vsNM, file="filepath")
asian_strata <- asian_strata_vsNM %>%
  dplyr::select(RR)

## Mixed
mixed_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_strata_vsNM)
save(mixed_strata_vsNM, file="filepath")
mixed_strata <- mixed_strata_vsNM %>%
  dplyr::select(RR)

## Other
other_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_strata_vsNM)
save(other_strata_vsNM, file="filepath")
other_strata <- other_strata_vsNM %>%
  dplyr::select(RR)

## Unknown
unknown_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_strata_vsNM)
save(unknown_strata_vsNM, file="filepath")
unknown_strata <- unknown_strata_vsNM %>%
  dplyr::select(RR)

# Combine results to make table
Ethnicity <- data.frame(c('White British', 'White non-British', 'Black', 'Asian', 'Mixed', 'Other', 'Unknown'))
migrants_by_ethnicity <- bind_rows(white_migrant, white_nonbritish_migrant, black_migrant, asian_migrant, mixed_migrant, other_migrant, unknown_migrant)
Ethniciy_strata <- bind_rows(white_strata, white_nonbritish_strata, black_strata, asian_strata, mixed_strata, other_strata, unknown_strata)
matched_ethnicity_table <- bind_cols(Ethnicity, migrants_by_ethnicity, Ethniciy_strata)
row.names(matched_ethnicity_table) <- NULL

matched_ethnicity_table <- matched_ethnicity_table %>%
  rename(
    Ethnicity = c..White.British....White.non.British....Black....Asian....Mixed...,
    RR_migrantsVsWBNM = RR...2,
    RR_MvsNM_each_ethnicity = RR...3
  ) %>%
  right_join(non_migrants_by_ethnicity, by = 'Ethnicity') %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, RR_migrantsVsWBNM, RR_MvsNM_each_ethnicity)

# Save dataframe for knitting in RMD file
save(matched_ethnicity_table, file="filepath")

### Additive and multiplicative interaction effects for each ethnicity 
# Based on Mathur MB & VanderWeele TJ (2018). R function for additive interaction measures. Epidemiology 29(1), e5-e6

# White non-British
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[3]
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_white_non_british <- data.frame(ethnicity = 'white_non_british', 
                                                    RERI = RERI,
                                                    RERI_upper_CI = upper_CI,
                                                    RERI_lower_CI = lower_CI,
                                                    multiplicative_effect = RR_interaction,
                                                    mult_lower_CI = beta_interaction_CI[1],
                                                    mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_white_non_british) <- NULL

# Black 
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[6]
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_black <- data.frame(ethnicity = 'black', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_black) <- NULL

# Asian
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[5] # Asian
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_asian <- data.frame(ethnicity = 'asian', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_asian) <- NULL

# Mixed
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[4] # Mixed
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_mixed <- data.frame(ethnicity = 'mixed', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_mixed) <- NULL

# Other
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[7] # Other
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_other <- data.frame(ethnicity = 'other', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_other) <- NULL

# Unknown
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[8] # Unknown
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_unknown <- data.frame(ethnicity = 'unknown', 
                                          RERI = RERI,
                                          RERI_upper_CI = upper_CI,
                                          RERI_lower_CI = lower_CI,
                                          multiplicative_effect = RR_interaction,
                                          mult_lower_CI = beta_interaction_CI[1],
                                          mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_unknown) <- NULL

# Combine RERIs and save
int_effects_all_ethnicities_ds <- bind_rows(interaction_effects_white_non_british, interaction_effects_black, 
                                         interaction_effects_asian, interaction_effects_mixed, 
                                         interaction_effects_other, interaction_effects_unknown)

write_csv(int_effects_all_ethnicities_ds, "filepath")

### London only  ----

## Multivariable model by migrant_status  ---

exact_match_London_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>%
  filter(prac_region == "London")

## All patients (i.e. all ages)

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = exact_match_London_annual_conscounts)
extract_glm_results_allages(x, 'London_glm_mig', exact_match_London_annual_conscounts, migrant_status)

write_csv(London_glm_mig_table, "filepath")


# IRR (multivariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts)
extract_glm_results_allages(x, 'London_multivariable_allages', exact_match_London_annual_conscounts, migrant_status)

## Age_subcohorts

# IRR (multivariable glm)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_0to15", "0-15 years", exact_match_London_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_16to24", "16-24 years", exact_match_London_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_25to34", "25-34 years", exact_match_London_annual_conscounts, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_35to49", "35-49 years", exact_match_London_annual_conscounts, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_50to64", "50-64 years", exact_match_London_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_65plus", ">=65 years", exact_match_London_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
London_multivariable_migrant_status_ds <- bind_rows(London_multivariable_allages_table, London_glm_multivariable_0to15_table,London_glm_multivariable_16to24_table,
                                                 London_glm_multivariable_25to34_table,London_glm_multivariable_35to49_table,London_glm_multivariable_50to64_table,London_glm_multivariable_65plus_table)


write_csv(London_multivariable_migrant_status_ds, "filepath")

### Study year  ----

## Multivariable IR and IRRs 

## 2015

studyyear_2015 <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% filter(studyyear == 2015)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015)
extract_glm_results_allages(x, 'multivariable_2015_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2015_migrant_status_ds <- bind_rows(multivariable_2015_allages_table, glm_2015_multivariable_0to15_table,glm_2015_multivariable_16to24_table,
                                               glm_2015_multivariable_25to34_table,glm_2015_multivariable_35to49_table,glm_2015_multivariable_50to64_table,glm_2015_multivariable_65plus_table)

write_csv(multivariable_2015_migrant_status_ds, "filepath")

## 2016

studyyear_2016 <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% filter(studyyear == 2016)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016)
extract_glm_results_allages(x, 'multivariable_2016_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2016_migrant_status_ds <- bind_rows(multivariable_2016_allages_table, glm_2016_multivariable_0to15_table,glm_2016_multivariable_16to24_table,
                                               glm_2016_multivariable_25to34_table,glm_2016_multivariable_35to49_table,glm_2016_multivariable_50to64_table,glm_2016_multivariable_65plus_table)


write_csv(multivariable_2016_migrant_status_ds, "filepath")

## 2017

studyyear_2017 <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% filter(studyyear == 2017)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017)
extract_glm_results_allages(x, 'multivariable_2017_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2017_migrant_status_ds <- bind_rows(multivariable_2017_allages_table, glm_2017_multivariable_0to15_table,glm_2017_multivariable_16to24_table,
                                               glm_2017_multivariable_25to34_table,glm_2017_multivariable_35to49_table,glm_2017_multivariable_50to64_table,glm_2017_multivariable_65plus_table)


write_csv(multivariable_2017_migrant_status_ds, "filepath")

## 2018

studyyear_2018 <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% filter(studyyear == 2018)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018)
extract_glm_results_allages(x, 'multivariable_2018_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2018_migrant_status_ds <- bind_rows(multivariable_2018_allages_table, glm_2018_multivariable_0to15_table,glm_2018_multivariable_16to24_table,
                                               glm_2018_multivariable_25to34_table,glm_2018_multivariable_35to49_table,glm_2018_multivariable_50to64_table,glm_2018_multivariable_65plus_table)


write_csv(multivariable_2018_migrant_status_ds, "filepath")

## 2019

studyyear_2019 <- matched_cohort_England_2015_2020_annual_conscounts_test_ds %>% filter(studyyear == 2019)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019)
extract_glm_results_allages(x, 'multivariable_2019_allages', matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_test_ds, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2019_migrant_status_ds <- bind_rows(multivariable_2019_allages_table, glm_2019_multivariable_0to15_table,glm_2019_multivariable_16to24_table,
                                               glm_2019_multivariable_25to34_table,glm_2019_multivariable_35to49_table,glm_2019_multivariable_50to64_table,glm_2019_multivariable_65plus_table)


write_csv(multivariable_2019_migrant_status_ds, "filepath")

# 5_Sensitivity analysis 2 -MATCHED ON: pyears and prac_region -------

## Load dataset and select variables ----

load(file = "filepath")

# Matched on age_data_start, year_data_start and prac_region 

matched_cohort_England_2015_2020_annual_conscounts_fu <- dplyr::select(matched_cohort_England_2015_2020_annual_conscounts_fu, 
                                                                            c(patid, pyears, conscount, migrant_status, migcertainty,
                                                                              age_subcohort, studyyear, prac_region, gender,
                                                                              ethnicat6, imd, studyyear_agecat))

# Relevel
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, gender <- relevel (gender, ref="Male")) 
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, imd <- relevel (imd, ref="IMD 1"))
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, prac_region <- relevel (prac_region, ref="London")) 
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, ethnicat6 <- relevel (ethnicat6, ref="White British"))

## Summary outcome measures (not included in paper) ------------------------------------------------------------------

### Overall (i.e. full period)

## All patients

# All patients
conscount_summary_overall_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% 
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migrant status
conscount_summary_mig_status_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migcertainty
conscount_summary_migcertainty_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(migcertainty) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Join for overall all patients + by migrant status
conscount_summary_overall_1 <- full_join(conscount_summary_overall_allpatients, conscount_summary_mig_status_allpatients, 
                                         by = c("n_events"="n_events", "mean" = "mean",
                                                "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals")) %>%
  full_join(conscount_summary_migcertainty_allpatients, by = c("n_events"="n_events", "mean" = "mean",
                                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                                                               "migrant_status" = "migcertainty","n_individuals"="n_individuals"))

conscount_summary_overall_1$migrant_status <- fct_explicit_na(conscount_summary_overall_1$migrant_status, na_level = "All_patients")
conscount_summary_overall_1$age <- 'all_ages'
conscount_summary_overall_1 <- conscount_summary_overall_1 %>% relocate(age, migrant_status)

## Age subcohorts

# Age_subcohorts
conscount_summary_overall_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(age_subcohort) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migrant status
conscount_summary_mig_status_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(age_subcohort, migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migcertainty
conscount_summary_migcertainty_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(age_subcohort,migcertainty) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Join for overall age_subcohorts + by migrant status
conscount_summary_overall_2 <- full_join(conscount_summary_overall_agesubcohorts, conscount_summary_mig_status_agesubcohorts, 
                                         by = c('age_subcohort' = 'age_subcohort', "n_events"="n_events", "mean" = "mean",
                                                "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals")) %>%
  full_join(conscount_summary_migcertainty_agesubcohorts, by = c('age_subcohort' = 'age_subcohort',"migrant_status" = "migcertainty","n_events"="n_events", "mean" = "mean",
                                                                 "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals"))
conscount_summary_overall_2$migrant_status <- fct_explicit_na(conscount_summary_overall_2$migrant_status, na_level = "All_patients")
conscount_summary_overall_2 <- conscount_summary_overall_2 %>% relocate(age_subcohort, migrant_status) %>%
  rename(age = age_subcohort) %>%
  arrange(age)

# Join all patients to the age_subcohorts
matched_conscount_summary_fullperiod_fu <- bind_rows(conscount_summary_overall_1,conscount_summary_overall_2)

# save file
write_csv(matched_conscount_summary_fullperiod_fu, "filepath" )

### Annual

## All patients

# All patients
conscount_annual_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migrant status
conscount_annual_mig_status_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migcertainty 
conscount_annual_migcertainty_allpatients <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(migcertainty, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))


# Join to give annual for allpatients + migrant status
conscount_summary_annual_1 <- full_join(conscount_annual_allpatients, conscount_annual_mig_status_allpatients, 
                                        by = c("studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max")) %>%
  full_join(conscount_annual_migcertainty_allpatients,
            by = c("studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                   "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                   "migrant_status" = "migcertainty"))

conscount_summary_annual_1$migrant_status <- fct_explicit_na(conscount_summary_annual_1$migrant_status, na_level = "All_patients") 
conscount_summary_annual_1$age <- 'all_ages'
conscount_summary_annual_1 <- conscount_summary_annual_1 %>% relocate(age, migrant_status)

## age_subcohorts

# age_subcohorts
conscount_annual_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(age_subcohort, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migrant status
conscount_annual_mig_status_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(age_subcohort,migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migcertainty 
conscount_annual_migcertainty_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% group_by(age_subcohort, migcertainty, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))


## Join to give annual for age_subcohorts + migrant status
conscount_summary_annual_2 <- full_join(conscount_annual_agesubcohorts, conscount_annual_mig_status_agesubcohorts, 
                                        by = c('age_subcohort' = 'age_subcohort', "studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max")) %>%
  full_join(conscount_annual_migcertainty_agesubcohorts,
            by = c('age_subcohort' = 'age_subcohort',"studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                   "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                   "migrant_status" = "migcertainty"))

conscount_summary_annual_2$migrant_status <- fct_explicit_na(conscount_summary_annual_2$migrant_status, na_level = "All_patients") 
conscount_summary_annual_2 <- conscount_summary_annual_2 %>% relocate(age_subcohort, migrant_status) %>%
  rename(age = age_subcohort) %>%
  arrange(age)

# Join all patients to the age_subcohorts
matched_conscount_summary_annual_fu <- bind_rows(conscount_summary_annual_1,conscount_summary_annual_2)

## Export  summary measures
write_csv(matched_conscount_summary_annual_fu, "filepath" )
 
## IRs by migrant_status + studyyear (not included in paper) -----

# IRs by migrant_status and studyyear

## All ages

pyears_allpatients_overall_migstatus_studyyear <- aggregate(pyears ~ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_allpatients_overall_migstatus_studyyear <- pyears_allpatients_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus_studyyear <- aggregate(conscount ~ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_allpatients_overall_migstatus_studyyear <- inner_join(conscount_allpatients_overall_migstatus_studyyear,pyears_allpatients_overall_migstatus_studyyear, by= c("migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migstatus_studyyear, output_name = 'IR_allpatients_overall_migstatus_studyyear')
IR_allpatients_overall_migstatus_studyyear$age_subcohort <- 'All_ages'

## age_subcohorts

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus_studyyear <- aggregate(pyears ~ age_subcohort + migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_agesubcohorts_overall_migstatus_studyyear <- pyears_agesubcohorts_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus_studyyear <- aggregate(conscount ~ age_subcohort+ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_agesubcohorts_overall_migstatus_studyyear <- inner_join(conscount_agesubcohorts_overall_migstatus_studyyear,pyears_agesubcohorts_overall_migstatus_studyyear, 
                                                           by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus_studyyear, output_name = 'IR_agesubcohorts_overall_migstatus_studyyear')

## Combine

IR_migstatus_studyyear <- full_join(IR_allpatients_overall_migstatus_studyyear, IR_agesubcohorts_overall_migstatus_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migrant_status' = 'migrant_status',
                                                                                                                                     'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))

## Reorder results (can change depending on how we want to present it)
IR_migstatus_studyyear <- arrange(IR_migstatus_studyyear, migrant_status, age_subcohort) %>%
  relocate(age_subcohort)

# IR_migstatus_studyyear <- arrange(IR_migstatus_studyyear, studyyear) %>%
#   relocate(age_subcohort)

write_csv(IR_migstatus_studyyear, "filepath")

## IR & univariable IRR by migrant status ----------------------------------------------------------------------

# Filter out 2020 for modelling
matched_cohort_England_2015_2020_annual_conscounts_fu <- matched_cohort_England_2015_2020_annual_conscounts_fu %>%
  filter(studyyear != 2020)

## All patients (i.e. all ages)

# IR by migrant status
pyears_allpatients_overall_migstatus <- aggregate(pyears ~ migrant_status, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_allpatients_overall_migstatus <- pyears_allpatients_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus <- aggregate(conscount ~ migrant_status, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_allpatients_overall_migstatus <- inner_join(conscount_allpatients_overall_migstatus,pyears_allpatients_overall_migstatus, by= c("migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_allpatients_overall_migstatus, output_name = 'IR_allpatients_overall_migstatus')
IR_allpatients_overall_migstatus$age_subcohort <- 'All_ages'

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = matched_cohort_England_2015_2020_annual_conscounts_fu)

extract_glm_results_allages(x, 'glm_mig', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)
round(ci.lin(x,Exp=T),3)

save(glm_mig, file="filepath")

# Join glm + IR + univariable_mig 
univariable_migrant_status <- full_join(IR_allpatients_overall_migstatus, glm_mig_table, by = c("age_subcohort" = "age_subcohort","migrant_status" = "names")) %>%
  relocate(age_subcohort)


## Age_subcohorts (not included in paper)

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus <- aggregate(pyears ~ age_subcohort + migrant_status, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_agesubcohorts_overall_migstatus <- pyears_agesubcohorts_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus <- aggregate(conscount ~ age_subcohort+ migrant_status, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_agesubcohorts_overall_migstatus <- inner_join(conscount_agesubcohorts_overall_migstatus,pyears_agesubcohorts_overall_migstatus, 
                                                 by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus, output_name = 'IR_agesubcohorts_overall_migstatus')

# IRR by migrant status (univariable glm)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_50to64", "50-64 years",matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 65 years and over  
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## NOTE IRR is comparing migrants at X age to non-migrants of the SAME age 

# Join age_subcohorts glm+IR and all_ages glm+IR 
glm_mig_combined_table <- bind_rows(glm_mig_0to15_table,glm_mig_16to24_table,glm_mig_25to34_table,glm_mig_35to49_table,glm_mig_50to64_table,glm_mig_65plus_table)

glm_mig_combined_table$names <- as.factor(glm_mig_combined_table$names)

univariable_migrant_status_agesubcohorts <- left_join(IR_agesubcohorts_overall_migstatus, glm_mig_combined_table, 
                                                      by = c("migrant_status" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migrant_status <- bind_rows(univariable_migrant_status, univariable_migrant_status_agesubcohorts)  
univariable_migrant_status$age_subcohort <- as.factor(univariable_migrant_status$age_subcohort)
univariable_migrant_status$age_subcohort <- factor(univariable_migrant_status$age_subcohort, 
                                                   levels = c('All_ages', '0-15 years','16-24 years', '25-34 years', '35-49 years','50-64 years', '>=65 years'))
## Reorder results (can change depending on how we want to present it)
univariable_migrant_status_fu <- arrange(univariable_migrant_status, migrant_status, age_subcohort)

write_csv(univariable_migrant_status_fu, "filepath")

## Multivariable model by migrant_status with IMD adjustment  ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)),
            data = matched_cohort_England_2015_2020_annual_conscounts_fu)
extract_glm_results_allages(x, 'multivariable_matched_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)
round(ci.lin(x,Exp=T),3)

## age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)
round(ci.lin(x,Exp=T),3)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)
round(ci.lin(x,Exp=T),3)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)
round(ci.lin(x,Exp=T),3)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)
round(ci.lin(x,Exp=T),3)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)
round(ci.lin(x,Exp=T),3)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)
round(ci.lin(x,Exp=T),3)

# Join all ages and age_subcohorts glm results
multivariable_matched_migrant_status_fu <- bind_rows(multivariable_matched_allages_table, glm_multivariable_matched_0to15_table,glm_multivariable_matched_16to24_table,
                                                     glm_multivariable_matched_25to34_table,glm_multivariable_matched_35to49_table,glm_multivariable_matched_50to64_table,glm_multivariable_matched_65plus_table)

write_csv(multivariable_matched_migrant_status_fu, "filepath")

# Save .Rdata for forestplots
save(multivariable_matched_allages, file="filepath")
save(glm_multivariable_matched_0to15, file="filepath")
save(glm_multivariable_matched_16to24, file="filepath")
save(glm_multivariable_matched_25to34, file="filepath")
save(glm_multivariable_matched_35to49, file="filepath")
save(glm_multivariable_matched_50to64, file="filepath")
save(glm_multivariable_matched_65plus, file="filepath")


## Analyses not included in paper ------
### Multivariable model by migrant_status without IMD adjustment  ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)),
            data = matched_cohort_England_2015_2020_annual_conscounts_fu)
extract_glm_results_allages(x, 'multivariable_matched_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_matched_migrant_status_noIMD_fu <- bind_rows(multivariable_matched_allages_table, glm_multivariable_matched_0to15_table,glm_multivariable_matched_16to24_table,
                                                           glm_multivariable_matched_25to34_table,glm_multivariable_matched_35to49_table,glm_multivariable_matched_50to64_table,glm_multivariable_matched_65plus_table)

write_csv(multivariable_matched_migrant_status_noIMD_fu, "filepath")

# Save .Rdata for forestplots
save(multivariable_matched_allages, file="filepath")
save(glm_multivariable_matched_0to15, file="filepath")
save(glm_multivariable_matched_16to24, file="filepath")
save(glm_multivariable_matched_25to34, file="filepath")
save(glm_multivariable_matched_35to49, file="filepath")
save(glm_multivariable_matched_50to64, file="filepath")
save(glm_multivariable_matched_65plus, file="filepath")

### Certainty of migration status  ----

# Reload to include 2020
load(file = "filepath")

matched_cohort_England_2015_2020_annual_conscounts_fu <- dplyr::select(matched_cohort_England_2015_2020_annual_conscounts_fu, 
                                                                            c(patid, pyears, conscount, migrant_status,
                                                                              age_subcohort, imd, ethnicat6, gender,
                                                                              studyyear, studyyear_agecat, prac_region,
                                                                              migcertainty, cohort_entry))

# Relevel
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, gender <- relevel (gender, ref="Male")) 
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, imd <- relevel (imd, ref="IMD 1"))
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, prac_region <- relevel (prac_region, ref="London")) 
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
matched_cohort_England_2015_2020_annual_conscounts_fu <- within(matched_cohort_England_2015_2020_annual_conscounts_fu, ethnicat6 <- relevel (ethnicat6, ref="White British"))

# IRs by migcertainty + studyyear

## All ages

pyears_allpatients_overall_migcertainty_studyyear <- aggregate(pyears ~ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_allpatients_overall_migcertainty_studyyear <- pyears_allpatients_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty_studyyear <- aggregate(conscount ~ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_allpatients_overall_migcertainty_studyyear <- inner_join(conscount_allpatients_overall_migcertainty_studyyear,pyears_allpatients_overall_migcertainty_studyyear, by= c("migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migcertainty_studyyear, output_name = 'IR_allpatients_overall_migcertainty_studyyear')
IR_allpatients_overall_migcertainty_studyyear$age_subcohort <- 'All_ages'

## Age_subcohorts

# IR by migcertainty
pyears_agesubcohorts_overall_migcertainty_studyyear <- aggregate(pyears ~ age_subcohort + migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_agesubcohorts_overall_migcertainty_studyyear <- pyears_agesubcohorts_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty_studyyear <- aggregate(conscount ~ age_subcohort+ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_agesubcohorts_overall_migcertainty_studyyear <- inner_join(conscount_agesubcohorts_overall_migcertainty_studyyear,pyears_agesubcohorts_overall_migcertainty_studyyear, 
                                                              by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty_studyyear, output_name = 'IR_agesubcohorts_overall_migcertainty_studyyear')

## Combine

IR_migcertainty_studyyear <- full_join(IR_allpatients_overall_migcertainty_studyyear, IR_agesubcohorts_overall_migcertainty_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migcertainty' = 'migcertainty',
                                                                                                                                              'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))
## Reorder results (can change depending on how we want to present it)
IR_migcertainty_studyyear_fu <- arrange(IR_migcertainty_studyyear, migcertainty, age_subcohort) %>%
  relocate(age_subcohort)

write_csv(IR_migcertainty_studyyear_fu, "filepath")

# Filter out 2020 for modelling
matched_cohort_England_2015_2020_annual_conscounts_fu <- matched_cohort_England_2015_2020_annual_conscounts_fu %>%
  filter(studyyear != 2020)

## Univariable model

## All patients (i.e. all ages)

# IR 
pyears_allpatients_overall_migcertainty <- aggregate(pyears ~ migcertainty, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_allpatients_overall_migcertainty <- pyears_allpatients_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty <- aggregate(conscount ~ migcertainty, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_allpatients_overall_migcertainty <- inner_join(conscount_allpatients_overall_migcertainty,pyears_allpatients_overall_migcertainty, by= c("migcertainty" = "migcertainty"))

function_to_calculate_IR(IR_allpatients_overall_migcertainty, output_name = 'IR_allpatients_overall_migcertainty')
IR_allpatients_overall_migcertainty$age_subcohort <- 'All_ages'

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), data = matched_cohort_England_2015_2020_annual_conscounts_fu)
extract_glm_results_allages(x, 'glm_migcertainty', matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

save(glm_migcertainty, file="filepath")

# Join glm + IR + univariable_mig 
univariable_migcertainty <- full_join(IR_allpatients_overall_migcertainty, glm_migcertainty_table, by = c("age_subcohort" = "age_subcohort","migcertainty" = "names")) %>%
  relocate(age_subcohort)

## Age_subcohorts

# IR
pyears_agesubcohorts_overall_migcertainty <- aggregate(pyears ~ age_subcohort + migcertainty, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
pyears_agesubcohorts_overall_migcertainty <- pyears_agesubcohorts_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty <- aggregate(conscount ~ age_subcohort+ migcertainty, matched_cohort_England_2015_2020_annual_conscounts_fu, sum) 
IR_agesubcohorts_overall_migcertainty <- inner_join(conscount_agesubcohorts_overall_migcertainty,pyears_agesubcohorts_overall_migcertainty, 
                                                    by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty"))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty, output_name = 'IR_agesubcohorts_overall_migcertainty')

# IRR (univariable glm)
### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

## 65 years and over NOTE  - haven't run this yet as the sample dataset doesn't have any >65s in
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)


# Join age_subcohorts glm+IR and all_ages glm+IR
glm_migcertainty_combined_table <- bind_rows(glm_migcertainty_0to15_table,glm_migcertainty_16to24_table,glm_migcertainty_25to34_table,glm_migcertainty_35to49_table,glm_migcertainty_50to64_table,glm_migcertainty_65plus_table)

univariable_migcertainty_agesubcohorts <- left_join(IR_agesubcohorts_overall_migcertainty, glm_migcertainty_combined_table, 
                                                    by = c("migcertainty" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migcertainty <- bind_rows(univariable_migcertainty, univariable_migcertainty_agesubcohorts)  
univariable_migcertainty$age_subcohort <- as.factor(univariable_migcertainty$age_subcohort)
univariable_migcertainty$age_subcohort <- factor(univariable_migcertainty$age_subcohort, 
                                                 levels = c('All_ages', '0-15 years','16-24 years', '25-34 years', '35-49 years','50-64 years', '>=65 years'))
univariable_migcertainty_fu <- arrange(univariable_migcertainty, age_subcohort)


write_csv(univariable_migcertainty_fu, "filepath")

## Multivariable model controlling for year + gender + as.factor(studyyear_agecat) + as.factor(imd) + prac_region

## All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu)
extract_glm_results_allages(x, 'multivariable_migcertainty_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts_fu[matched_cohort_England_2015_2020_annual_conscounts_fu$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migcertainty)

# Join all ages and age_subcohorts glm results
multivariable_migcertainty_fu <- bind_rows(multivariable_migcertainty_allages_table, glm_multivariable_migcertainty_0to15_table,glm_multivariable_migcertainty_16to24_table,
                                           glm_multivariable_migcertainty_25to34_table,glm_multivariable_migcertainty_35to49_table,glm_multivariable_migcertainty_50to64_table,glm_multivariable_migcertainty_65plus_table)


write_csv(multivariable_migcertainty_fu, "filepath")

# Save .Rdata for forestplots
save(multivariable_migcertainty_allages, file="filepath")
save(glm_multivariable_migcertainty_0to15, file="filepath")
save(glm_multivariable_migcertainty_16to24, file="filepath")
save(glm_multivariable_migcertainty_25to34, file="filepath")
save(glm_multivariable_migcertainty_35to49, file="filepath")
save(glm_multivariable_migcertainty_50to64, file="filepath")
save(glm_multivariable_migcertainty_65plus, file="filepath")


### Ethnicity (ethnicat6) - interaction term, additive and multiplicative effects ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status)*as.factor(ethnicat6) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = matched_cohort_England_2015_2020_annual_conscounts_fu)
extract_glm_results_allages(x, 'multivariable_ethnicat6interaction_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# Save .Rdata for forestplots
save(multivariable_ethnicat6interaction_allages, file="filepath")
write_csv(multivariable_ethnicat6interaction_allages, "filepath")

# RR (95% CI) for non-migrants in each ethnicity compared to White British non-migrants
# For table
white_non_migrant <- data.frame(Ethnicity = 'White British', RR_nonmigrantsVsWBNM = '1.0')
non_migrants_by_ethnicity <- multivariable_ethnicat6interaction_allages[3:8,] %>%
  mutate(ci = paste(lower, upper, sep ="-"),
         RR_nonmigrantsVsWBNM = paste0(estimate, ' (', ci, ")")) %>%
  mutate(Ethnicity = c('White non-British', 'Mixed', 'Asian', 'Black', 'Other', 'Unknown')) %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, estimate, upper, lower)
non_migrants_by_ethnicity <- bind_rows(white_non_migrant, non_migrants_by_ethnicity)

# Calculate interaction effects - migrants in each ethnicity compared to White British non-migrant (WBNM) reference group

get_RR_and_95_CI <- function(dataframe) {
  output <- dataframe %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$Lower.CI, output$Upper.CI, sep ="-")
  output$RR <- paste0(output$Estimate, ' (', output$ci, ")")
  #output <- output %>%
  #  dplyr:: select(RR)
  output_name <- deparse(substitute(dataframe))
  assign(x = output_name, value = output, envir = globalenv())
}

## White British migrants 
white_migrant_vsWBNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_migrant_vsWBNM)
save(white_migrant_vsWBNM, file="filepath")
white_migrant <- white_migrant_vsWBNM %>%
  dplyr::select(RR)

# White non-British migrants 
white_nonbritish_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)White non-British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_migrant_vsWBNM)
save(white_nonbritish_migrant_vsWBNM, file="filepath")
white_nonbritish_migrant <- white_nonbritish_migrant_vsWBNM %>%
  dplyr::select(RR)

## Black migrants 
black_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_migrant_vsWBNM)
save(black_migrant_vsWBNM, file="filepath")
black_migrant <- black_migrant_vsWBNM %>%
  dplyr::select(RR)

## Asian
asian_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Asian/Asian British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_migrant_vsWBNM)
save(asian_migrant_vsWBNM, file="filepath")
asian_migrant <- asian_migrant_vsWBNM %>%
  dplyr::select(RR)

## Mixed
mixed_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_migrant_vsWBNM)
save(mixed_migrant_vsWBNM, file="filepath")
mixed_migrant <- mixed_migrant_vsWBNM %>%
  dplyr::select(RR)

## Other
other_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Other ethnic group' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_migrant_vsWBNM)
save(other_migrant_vsWBNM, file="filepath")
other_migrant <- other_migrant_vsWBNM %>%
  dplyr::select(RR)

## Unknown
unknown_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Unknown' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_migrant_vsWBNM)
save(unknown_migrant_vsWBNM, file="filepath")
unknown_migrant <- unknown_migrant_vsWBNM %>%
  dplyr::select(RR)

# Calculate interaction effects - migrants vs non-migrants in each ethnicity strata
## White
white_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_strata_vsNM)
save(white_strata_vsNM, file="filepath")
white_strata <- white_strata_vsNM %>%
  dplyr::select(RR)

## White non-British
white_nonbritish_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_strata_vsNM)
save(white_nonbritish_strata_vsNM, file="filepath")
white_nonbritish_strata <- white_nonbritish_strata_vsNM %>%
  dplyr::select(RR)

## Black
black_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_strata_vsNM)
save(black_strata_vsNM, file="filepath")
black_strata <- black_strata_vsNM %>%
  dplyr::select(RR)

## Asian
asian_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_strata_vsNM)
save(asian_strata_vsNM, file="filepath")
asian_strata <- asian_strata_vsNM %>%
  dplyr::select(RR)

## Mixed
mixed_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_strata_vsNM)
save(mixed_strata_vsNM, file="filepath")
mixed_strata <- mixed_strata_vsNM %>%
  dplyr::select(RR)

## Other
other_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_strata_vsNM)
save(other_strata_vsNM, file="filepath")
other_strata <- other_strata_vsNM %>%
  dplyr::select(RR)

## Unknown
unknown_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_strata_vsNM)
save(unknown_strata_vsNM, file="filepath")
unknown_strata <- unknown_strata_vsNM %>%
  dplyr::select(RR)

# Combine results to make table
Ethnicity <- data.frame(c('White British', 'White non-British', 'Black', 'Asian', 'Mixed', 'Other', 'Unknown'))
migrants_by_ethnicity <- bind_rows(white_migrant, white_nonbritish_migrant, black_migrant, asian_migrant, mixed_migrant, other_migrant, unknown_migrant)
Ethniciy_strata <- bind_rows(white_strata, white_nonbritish_strata, black_strata, asian_strata, mixed_strata, other_strata, unknown_strata)
matched_ethnicity_table <- bind_cols(Ethnicity, migrants_by_ethnicity, Ethniciy_strata)
row.names(matched_ethnicity_table) <- NULL

matched_ethnicity_table <- matched_ethnicity_table %>%
  rename(
    Ethnicity = c..White.British....White.non.British....Black....Asian....Mixed...,
    RR_migrantsVsWBNM = RR...2,
    RR_MvsNM_each_ethnicity = RR...3
  ) %>%
  right_join(non_migrants_by_ethnicity, by = 'Ethnicity') %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, RR_migrantsVsWBNM, RR_MvsNM_each_ethnicity)

# Save dataframe for knitting in RMD file
save(matched_ethnicity_table, file="filepath")

### Additive and multiplicative interaction effects for each ethnicity 
# Based on Mathur MB & VanderWeele TJ (2018). R function for additive interaction measures. Epidemiology 29(1), e5-e6

# White non-British
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[3]
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_white_non_british <- data.frame(ethnicity = 'white_non_british', 
                                                    RERI = RERI,
                                                    RERI_upper_CI = upper_CI,
                                                    RERI_lower_CI = lower_CI,
                                                    multiplicative_effect = RR_interaction,
                                                    mult_lower_CI = beta_interaction_CI[1],
                                                    mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_white_non_british) <- NULL

# Black 
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[6]
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_black <- data.frame(ethnicity = 'black', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_black) <- NULL

# Asian
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[5] # Asian
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_asian <- data.frame(ethnicity = 'asian', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_asian) <- NULL

# Mixed
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[4] # Mixed
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_mixed <- data.frame(ethnicity = 'mixed', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_mixed) <- NULL

# Other
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[7] # Other
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_other <- data.frame(ethnicity = 'other', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_other) <- NULL

# Unknown
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[8] # Unknown
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_unknown <- data.frame(ethnicity = 'unknown', 
                                          RERI = RERI,
                                          RERI_upper_CI = upper_CI,
                                          RERI_lower_CI = lower_CI,
                                          multiplicative_effect = RR_interaction,
                                          mult_lower_CI = beta_interaction_CI[1],
                                          mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_unknown) <- NULL

# Combine RERIs and save
int_effects_all_ethnicities_fu <- bind_rows(interaction_effects_white_non_british, interaction_effects_black, 
                                            interaction_effects_asian, interaction_effects_mixed, 
                                            interaction_effects_other, interaction_effects_unknown)

write_csv(int_effects_all_ethnicities_fu, "filepath")

### London only  ----

## Multivariable model by migrant_status  ---

exact_match_London_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts_fu %>%
  filter(prac_region == "London")

## All patients (i.e. all ages)

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = exact_match_London_annual_conscounts)
extract_glm_results_allages(x, 'London_glm_mig', exact_match_London_annual_conscounts, migrant_status)

write_csv(London_glm_mig_table, "filepath")


# IRR (multivariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts)
extract_glm_results_allages(x, 'London_multivariable_allages', exact_match_London_annual_conscounts, migrant_status)

## Age_subcohorts

# IRR (multivariable glm)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_0to15", "0-15 years", exact_match_London_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_16to24", "16-24 years", exact_match_London_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_25to34", "25-34 years", exact_match_London_annual_conscounts, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_35to49", "35-49 years", exact_match_London_annual_conscounts, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_50to64", "50-64 years", exact_match_London_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_65plus", ">=65 years", exact_match_London_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
London_multivariable_migrant_status_fu <- bind_rows(London_multivariable_allages_table, London_glm_multivariable_0to15_table,London_glm_multivariable_16to24_table,
                                                    London_glm_multivariable_25to34_table,London_glm_multivariable_35to49_table,London_glm_multivariable_50to64_table,London_glm_multivariable_65plus_table)


write_csv(London_multivariable_migrant_status_fu, "filepath")

### Study year  ----

## Multivariable IR and IRRs 

## 2015

studyyear_2015 <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% filter(studyyear == 2015)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015)
extract_glm_results_allages(x, 'multivariable_2015_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2015_migrant_status_fu <- bind_rows(multivariable_2015_allages_table, glm_2015_multivariable_0to15_table,glm_2015_multivariable_16to24_table,
                                                  glm_2015_multivariable_25to34_table,glm_2015_multivariable_35to49_table,glm_2015_multivariable_50to64_table,glm_2015_multivariable_65plus_table)

write_csv(multivariable_2015_migrant_status_fu, "filepath")

## 2016

studyyear_2016 <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% filter(studyyear == 2016)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016)
extract_glm_results_allages(x, 'multivariable_2016_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2016_migrant_status_fu <- bind_rows(multivariable_2016_allages_table, glm_2016_multivariable_0to15_table,glm_2016_multivariable_16to24_table,
                                                  glm_2016_multivariable_25to34_table,glm_2016_multivariable_35to49_table,glm_2016_multivariable_50to64_table,glm_2016_multivariable_65plus_table)


write_csv(multivariable_2016_migrant_status_fu, "filepath")

## 2017

studyyear_2017 <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% filter(studyyear == 2017)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017)
extract_glm_results_allages(x, 'multivariable_2017_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2017_migrant_status_fu <- bind_rows(multivariable_2017_allages_table, glm_2017_multivariable_0to15_table,glm_2017_multivariable_16to24_table,
                                                  glm_2017_multivariable_25to34_table,glm_2017_multivariable_35to49_table,glm_2017_multivariable_50to64_table,glm_2017_multivariable_65plus_table)


write_csv(multivariable_2017_migrant_status_fu, "filepath")

## 2018

studyyear_2018 <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% filter(studyyear == 2018)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018)
extract_glm_results_allages(x, 'multivariable_2018_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2018_migrant_status_fu <- bind_rows(multivariable_2018_allages_table, glm_2018_multivariable_0to15_table,glm_2018_multivariable_16to24_table,
                                                  glm_2018_multivariable_25to34_table,glm_2018_multivariable_35to49_table,glm_2018_multivariable_50to64_table,glm_2018_multivariable_65plus_table)


write_csv(multivariable_2018_migrant_status_fu, "filepath")

## 2019

studyyear_2019 <- matched_cohort_England_2015_2020_annual_conscounts_fu %>% filter(studyyear == 2019)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019)
extract_glm_results_allages(x, 'multivariable_2019_allages', matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts_fu, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts_fu, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts_fu, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2019_migrant_status_fu <- bind_rows(multivariable_2019_allages_table, glm_2019_multivariable_0to15_table,glm_2019_multivariable_16to24_table,
                                                  glm_2019_multivariable_25to34_table,glm_2019_multivariable_35to49_table,glm_2019_multivariable_50to64_table,glm_2019_multivariable_65plus_table)


write_csv(multivariable_2019_migrant_status_fu, "filepath")

# 6_Further analyses not included in paper -------
## Sensitivity analysis matched on age_cohort_entry, year_cohort_entry and prac_region
## Load dataset and select required variables 

# Load dataset 

load(file = "filepath")

# Matched on age_cohort_entry, year_cohort_entry and prac_region 

matched_cohort_England_2015_2020_annual_conscounts <- dplyr::select(matched_cohort_England_2015_2020_annual_conscounts, 
                                                                    c(patid, pyears, conscount, migrant_status, migcertainty,
                                                                      age_subcohort, studyyear, prac_region, gender,
                                                                      ethnicat6, imd, studyyear_agecat))

# Relevel
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, gender <- relevel (gender, ref="Male")) 
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, imd <- relevel (imd, ref="IMD 1"))
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, prac_region <- relevel (prac_region, ref="London")) 
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, ethnicat6 <- relevel (ethnicat6, ref="White British"))

## MAIN ANALYSIS ### -----

## 5_Summary outcome measures ------------------------------------------------------------------

### Overall (i.e. full period)

## All patients

# All patients
conscount_summary_overall_allpatients <- matched_cohort_England_2015_2020_annual_conscounts %>% 
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migrant status
conscount_summary_mig_status_allpatients <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migcertainty
conscount_summary_migcertainty_allpatients <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(migcertainty) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Join for overall all patients + by migrant status
conscount_summary_overall_1 <- full_join(conscount_summary_overall_allpatients, conscount_summary_mig_status_allpatients, 
                                         by = c("n_events"="n_events", "mean" = "mean",
                                                "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals")) %>%
  full_join(conscount_summary_migcertainty_allpatients, by = c("n_events"="n_events", "mean" = "mean",
                                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                                                               "migrant_status" = "migcertainty","n_individuals"="n_individuals"))

conscount_summary_overall_1$migrant_status <- fct_explicit_na(conscount_summary_overall_1$migrant_status, na_level = "All_patients")
conscount_summary_overall_1$age <- 'all_ages'
conscount_summary_overall_1 <- conscount_summary_overall_1 %>% relocate(age, migrant_status)

## Age subcohorts

# Age_subcohorts
conscount_summary_overall_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migrant status
conscount_summary_mig_status_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort, migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migcertainty
conscount_summary_migcertainty_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort,migcertainty) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Join for overall age_subcohorts + by migrant status
conscount_summary_overall_2 <- full_join(conscount_summary_overall_agesubcohorts, conscount_summary_mig_status_agesubcohorts, 
                                         by = c('age_subcohort' = 'age_subcohort', "n_events"="n_events", "mean" = "mean",
                                                "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals")) %>%
  full_join(conscount_summary_migcertainty_agesubcohorts, by = c('age_subcohort' = 'age_subcohort',"migrant_status" = "migcertainty","n_events"="n_events", "mean" = "mean",
                                                                 "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max","n_individuals"="n_individuals"))
conscount_summary_overall_2$migrant_status <- fct_explicit_na(conscount_summary_overall_2$migrant_status, na_level = "All_patients")
conscount_summary_overall_2 <- conscount_summary_overall_2 %>% relocate(age_subcohort, migrant_status) %>%
  rename(age = age_subcohort) %>%
  arrange(age)

# Join all patients to the age_subcohorts
matched_conscount_summary_fullperiod <- bind_rows(conscount_summary_overall_1,conscount_summary_overall_2)

# save file
write_csv(matched_conscount_summary_fullperiod, "filepath" )
#view(conscount_summary_overall_final)

### Annual

## All patients

# All patients
conscount_annual_allpatients <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migrant status
conscount_annual_mig_status_allpatients <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migcertainty 
conscount_annual_migcertainty_allpatients <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(migcertainty, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))


# Join to give annual for allpatients + migrant status
conscount_summary_annual_1 <- full_join(conscount_annual_allpatients, conscount_annual_mig_status_allpatients, 
                                        by = c("studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max")) %>%
  full_join(conscount_annual_migcertainty_allpatients,
            by = c("studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                   "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                   "migrant_status" = "migcertainty"))

conscount_summary_annual_1$migrant_status <- fct_explicit_na(conscount_summary_annual_1$migrant_status, na_level = "All_patients") 
conscount_summary_annual_1$age <- 'all_ages'
conscount_summary_annual_1 <- conscount_summary_annual_1 %>% relocate(age, migrant_status)

## age_subcohorts

# age_subcohorts
conscount_annual_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migrant status
conscount_annual_mig_status_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort,migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migcertainty 
conscount_annual_migcertainty_agesubcohorts <- matched_cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort, migcertainty, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))


## Join to give annual for age_subcohorts + migrant status
conscount_summary_annual_2 <- full_join(conscount_annual_agesubcohorts, conscount_annual_mig_status_agesubcohorts, 
                                        by = c('age_subcohort' = 'age_subcohort', "studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                                               "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max")) %>%
  full_join(conscount_annual_migcertainty_agesubcohorts,
            by = c('age_subcohort' = 'age_subcohort',"studyyear" = "studyyear", "n_events"="n_events", "mean" = "mean",
                   "sd"="sd", "median"="median", "iqr"="iqr", "min"="min", "max"="max",
                   "migrant_status" = "migcertainty"))

conscount_summary_annual_2$migrant_status <- fct_explicit_na(conscount_summary_annual_2$migrant_status, na_level = "All_patients") 
conscount_summary_annual_2 <- conscount_summary_annual_2 %>% relocate(age_subcohort, migrant_status) %>%
  rename(age = age_subcohort) %>%
  arrange(age)

# Join all patients to the age_subcohorts
matched_conscount_summary_annual <- bind_rows(conscount_summary_annual_1,conscount_summary_annual_2)

## Export  summary measures
write_csv(matched_conscount_summary_annual, "filepath" )

## 6_IRs by migrant_status + studyyear -----

# IRs by migrant_status and studyyear

## All ages

pyears_allpatients_overall_migstatus_studyyear <- aggregate(pyears ~ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migstatus_studyyear <- pyears_allpatients_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus_studyyear <- aggregate(conscount ~ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migstatus_studyyear <- inner_join(conscount_allpatients_overall_migstatus_studyyear,pyears_allpatients_overall_migstatus_studyyear, by= c("migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migstatus_studyyear, output_name = 'IR_allpatients_overall_migstatus_studyyear')
IR_allpatients_overall_migstatus_studyyear$age_subcohort <- 'All_ages'

## age_subcohorts

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus_studyyear <- aggregate(pyears ~ age_subcohort + migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_agesubcohorts_overall_migstatus_studyyear <- pyears_agesubcohorts_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus_studyyear <- aggregate(conscount ~ age_subcohort+ migrant_status + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_agesubcohorts_overall_migstatus_studyyear <- inner_join(conscount_agesubcohorts_overall_migstatus_studyyear,pyears_agesubcohorts_overall_migstatus_studyyear, 
                                                           by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus_studyyear, output_name = 'IR_agesubcohorts_overall_migstatus_studyyear')

## Combine

IR_migstatus_studyyear <- full_join(IR_allpatients_overall_migstatus_studyyear, IR_agesubcohorts_overall_migstatus_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migrant_status' = 'migrant_status',
                                                                                                                                     'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))

## Reorder results (can change depending on how we want to present it)
IR_migstatus_studyyear <- arrange(IR_migstatus_studyyear, migrant_status, age_subcohort) %>%
  relocate(age_subcohort)

# IR_migstatus_studyyear <- arrange(IR_migstatus_studyyear, studyyear) %>%
#   relocate(age_subcohort)

write_csv(IR_migstatus_studyyear, "filepath")

## 7_IR & univariable IRR by migrant status ----------------------------------------------------------------------

# Filter out 2020 for modelling
matched_cohort_England_2015_2020_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  filter(studyyear != 2020)

## All patients (i.e. all ages)

# IR by migrant status
pyears_allpatients_overall_migstatus <- aggregate(pyears ~ migrant_status, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migstatus <- pyears_allpatients_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus <- aggregate(conscount ~ migrant_status, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migstatus <- inner_join(conscount_allpatients_overall_migstatus,pyears_allpatients_overall_migstatus, by= c("migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_allpatients_overall_migstatus, output_name = 'IR_allpatients_overall_migstatus')
IR_allpatients_overall_migstatus$age_subcohort <- 'All_ages'

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = matched_cohort_England_2015_2020_annual_conscounts)

extract_glm_results_allages(x, 'glm_mig', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)


save(glm_mig, file="filepath")

# Join glm + IR + univariable_mig 
univariable_migrant_status <- full_join(IR_allpatients_overall_migstatus, glm_mig_table, by = c("age_subcohort" = "age_subcohort","migrant_status" = "names")) %>%
  relocate(age_subcohort)


## Age_subcohorts

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus <- aggregate(pyears ~ age_subcohort + migrant_status, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_agesubcohorts_overall_migstatus <- pyears_agesubcohorts_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus <- aggregate(conscount ~ age_subcohort+ migrant_status, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_agesubcohorts_overall_migstatus <- inner_join(conscount_agesubcohorts_overall_migstatus,pyears_agesubcohorts_overall_migstatus, 
                                                 by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus, output_name = 'IR_agesubcohorts_overall_migstatus')

# IRR by migrant status (univariable glm)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_50to64", "50-64 years",matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over  
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## NOTE IRR is comparing migrants at X age to non-migrants of the SAME age 

# Join age_subcohorts glm+IR and all_ages glm+IR 
glm_mig_combined_table <- bind_rows(glm_mig_0to15_table,glm_mig_16to24_table,glm_mig_25to34_table,glm_mig_35to49_table,glm_mig_50to64_table,glm_mig_65plus_table)

glm_mig_combined_table$names <- as.factor(glm_mig_combined_table$names)

univariable_migrant_status_agesubcohorts <- left_join(IR_agesubcohorts_overall_migstatus, glm_mig_combined_table, 
                                                      by = c("migrant_status" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migrant_status <- bind_rows(univariable_migrant_status, univariable_migrant_status_agesubcohorts)  
univariable_migrant_status$age_subcohort <- as.factor(univariable_migrant_status$age_subcohort)
univariable_migrant_status$age_subcohort <- factor(univariable_migrant_status$age_subcohort, 
                                                   levels = c('All_ages', '0-15 years','16-24 years', '25-34 years', '35-49 years','50-64 years', '>=65 years'))
## Reorder results (can change depending on how we want to present it)
univariable_migrant_status <- arrange(univariable_migrant_status, migrant_status, age_subcohort)

write_csv(univariable_migrant_status, "filepath")

## 8a_Multivariable model by migrant_status with IMD adjustment  ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)),
            data = matched_cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'multivariable_matched_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_matched_migrant_status <- bind_rows(multivariable_matched_allages_table, glm_multivariable_matched_0to15_table,glm_multivariable_matched_16to24_table,
                                                  glm_multivariable_matched_25to34_table,glm_multivariable_matched_35to49_table,glm_multivariable_matched_50to64_table,glm_multivariable_matched_65plus_table)

write_csv(multivariable_matched_migrant_status, "filepath")

# Save .Rdata for forestplots
save(multivariable_matched_allages, file="filepath")
save(glm_multivariable_matched_0to15, file="filepath")
save(glm_multivariable_matched_16to24, file="filepath")
save(glm_multivariable_matched_25to34, file="filepath")
save(glm_multivariable_matched_35to49, file="filepath")
save(glm_multivariable_matched_50to64, file="filepath")
save(glm_multivariable_matched_65plus, file="filepath")


## 8b_Multivariable model by migrant_status without IMD adjustment  ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)),
            data = matched_cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'multivariable_matched_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + studyyear_agecat  + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_matched_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_matched_migrant_status <- bind_rows(multivariable_matched_allages_table, glm_multivariable_matched_0to15_table,glm_multivariable_matched_16to24_table,
                                                  glm_multivariable_matched_25to34_table,glm_multivariable_matched_35to49_table,glm_multivariable_matched_50to64_table,glm_multivariable_matched_65plus_table)


multivariable_matched_migrant_status_noIMD <- multivariable_matched_migrant_status
write_csv(multivariable_matched_migrant_status_noIMD, "filepath")

# Save .Rdata for forestplots
save(multivariable_matched_allages, file="filepath")
save(glm_multivariable_matched_0to15, file="filepath")
save(glm_multivariable_matched_16to24, file="filepath")
save(glm_multivariable_matched_25to34, file="filepath")
save(glm_multivariable_matched_35to49, file="filepath")
save(glm_multivariable_matched_50to64, file="filepath")
save(glm_multivariable_matched_65plus, file="filepath")

# SENSITIVITY ANALYSES ## ----

## 9_Certainty of migration status  ----

# Reload to include 2020
load(file = "filepath")

matched_cohort_England_2015_2020_annual_conscounts <- dplyr::select(matched_cohort_England_2015_2020_annual_conscounts, 
                                                                    c(patid, pyears, conscount, migrant_status,
                                                                      age_subcohort, imd, ethnicat6, gender,
                                                                      studyyear, studyyear_agecat, prac_region,
                                                                      migcertainty, cohort_entry))

# Relevel
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, gender <- relevel (gender, ref="Male")) 
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, imd <- relevel (imd, ref="IMD 1"))
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, prac_region <- relevel (prac_region, ref="London")) 
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
matched_cohort_England_2015_2020_annual_conscounts <- within(matched_cohort_England_2015_2020_annual_conscounts, ethnicat6 <- relevel (ethnicat6, ref="White British"))

# IRs by migcertainty + studyyear

## All ages

pyears_allpatients_overall_migcertainty_studyyear <- aggregate(pyears ~ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migcertainty_studyyear <- pyears_allpatients_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty_studyyear <- aggregate(conscount ~ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migcertainty_studyyear <- inner_join(conscount_allpatients_overall_migcertainty_studyyear,pyears_allpatients_overall_migcertainty_studyyear, by= c("migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migcertainty_studyyear, output_name = 'IR_allpatients_overall_migcertainty_studyyear')
IR_allpatients_overall_migcertainty_studyyear$age_subcohort <- 'All_ages'

## Age_subcohorts

# IR by migcertainty
pyears_agesubcohorts_overall_migcertainty_studyyear <- aggregate(pyears ~ age_subcohort + migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_agesubcohorts_overall_migcertainty_studyyear <- pyears_agesubcohorts_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty_studyyear <- aggregate(conscount ~ age_subcohort+ migcertainty + studyyear, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_agesubcohorts_overall_migcertainty_studyyear <- inner_join(conscount_agesubcohorts_overall_migcertainty_studyyear,pyears_agesubcohorts_overall_migcertainty_studyyear, 
                                                              by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty_studyyear, output_name = 'IR_agesubcohorts_overall_migcertainty_studyyear')

## Combine

IR_migcertainty_studyyear <- full_join(IR_allpatients_overall_migcertainty_studyyear, IR_agesubcohorts_overall_migcertainty_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migcertainty' = 'migcertainty',
                                                                                                                                              'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))
## Reorder results (can change depending on how we want to present it)
IR_migcertainty_studyyear <- arrange(IR_migcertainty_studyyear, migcertainty, age_subcohort) %>%
  relocate(age_subcohort)

# IR_migcertainty_studyyear <- arrange(IR_migcertainty_studyyear, studyyear) %>%
#   relocate(age_subcohort)


write_csv(IR_migcertainty_studyyear, "filepath")

# Filter out 2020 for modelling
matched_cohort_England_2015_2020_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  filter(studyyear != 2020)

## Univariable model

## All patients (i.e. all ages)

# IR 
pyears_allpatients_overall_migcertainty <- aggregate(pyears ~ migcertainty, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migcertainty <- pyears_allpatients_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty <- aggregate(conscount ~ migcertainty, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migcertainty <- inner_join(conscount_allpatients_overall_migcertainty,pyears_allpatients_overall_migcertainty, by= c("migcertainty" = "migcertainty"))

function_to_calculate_IR(IR_allpatients_overall_migcertainty, output_name = 'IR_allpatients_overall_migcertainty')
IR_allpatients_overall_migcertainty$age_subcohort <- 'All_ages'

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), data = matched_cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'glm_migcertainty', matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

save(glm_migcertainty, file="filepath")

# Join glm + IR + univariable_mig 
univariable_migcertainty <- full_join(IR_allpatients_overall_migcertainty, glm_migcertainty_table, by = c("age_subcohort" = "age_subcohort","migcertainty" = "names")) %>%
  relocate(age_subcohort)

## Age_subcohorts

# IR
pyears_agesubcohorts_overall_migcertainty <- aggregate(pyears ~ age_subcohort + migcertainty, matched_cohort_England_2015_2020_annual_conscounts, sum) 
pyears_agesubcohorts_overall_migcertainty <- pyears_agesubcohorts_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty <- aggregate(conscount ~ age_subcohort+ migcertainty, matched_cohort_England_2015_2020_annual_conscounts, sum) 
IR_agesubcohorts_overall_migcertainty <- inner_join(conscount_agesubcohorts_overall_migcertainty,pyears_agesubcohorts_overall_migcertainty, 
                                                    by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty"))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty, output_name = 'IR_agesubcohorts_overall_migcertainty')

# IRR (univariable glm)
### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

## 65 years and over NOTE  - haven't run this yet as the sample dataset doesn't have any >65s in
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)


# Join age_subcohorts glm+IR and all_ages glm+IR
glm_migcertainty_combined_table <- bind_rows(glm_migcertainty_0to15_table,glm_migcertainty_16to24_table,glm_migcertainty_25to34_table,glm_migcertainty_35to49_table,glm_migcertainty_50to64_table,glm_migcertainty_65plus_table)

univariable_migcertainty_agesubcohorts <- left_join(IR_agesubcohorts_overall_migcertainty, glm_migcertainty_combined_table, 
                                                    by = c("migcertainty" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migcertainty <- bind_rows(univariable_migcertainty, univariable_migcertainty_agesubcohorts)  
univariable_migcertainty$age_subcohort <- as.factor(univariable_migcertainty$age_subcohort)
univariable_migcertainty$age_subcohort <- factor(univariable_migcertainty$age_subcohort, 
                                                 levels = c('All_ages', '0-15 years','16-24 years', '25-34 years', '35-49 years','50-64 years', '>=65 years'))
univariable_migcertainty <- arrange(univariable_migcertainty, age_subcohort)


write_csv(univariable_migcertainty, "filepath")

## Multivariable model controlling for year + gender + as.factor(studyyear_agecat) + as.factor(imd) + prac_region

## All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'multivariable_migcertainty_allages', matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = matched_cohort_England_2015_2020_annual_conscounts[matched_cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migcertainty)

# Join all ages and age_subcohorts glm results
multivariable_migcertainty <- bind_rows(multivariable_migcertainty_allages_table, glm_multivariable_migcertainty_0to15_table,glm_multivariable_migcertainty_16to24_table,
                                        glm_multivariable_migcertainty_25to34_table,glm_multivariable_migcertainty_35to49_table,glm_multivariable_migcertainty_50to64_table,glm_multivariable_migcertainty_65plus_table)


write_csv(multivariable_migcertainty, "filepath")

# Save .Rdata for forestplots
save(multivariable_migcertainty_allages, file="filepath")
save(glm_multivariable_migcertainty_0to15, file="filepath")
save(glm_multivariable_migcertainty_16to24, file="filepath")
save(glm_multivariable_migcertainty_25to34, file="filepath")
save(glm_multivariable_migcertainty_35to49, file="filepath")
save(glm_multivariable_migcertainty_50to64, file="filepath")
save(glm_multivariable_migcertainty_65plus, file="filepath")


## 10_Ethnicity (ethnicat6) - interaction term, additive and multiplicative effects ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status)*as.factor(ethnicat6) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = matched_cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'multivariable_ethnicat6interaction_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Save .Rdata for forestplots
save(multivariable_ethnicat6interaction_allages, file="filepath")
write_csv(multivariable_ethnicat6interaction_allages, "filepath")

# RR (95% CI) for non-migrants in each ethnicity compared to White British non-migrants
# For table
white_non_migrant <- data.frame(Ethnicity = 'White British', RR_nonmigrantsVsWBNM = '1.0')
non_migrants_by_ethnicity <- multivariable_ethnicat6interaction_allages[3:8,] %>%
  mutate(ci = paste(lower, upper, sep ="-"),
         RR_nonmigrantsVsWBNM = paste0(estimate, ' (', ci, ")")) %>%
  mutate(Ethnicity = c('White non-British', 'Mixed', 'Asian', 'Black', 'Other', 'Unknown')) %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, estimate, upper, lower)
non_migrants_by_ethnicity <- bind_rows(white_non_migrant, non_migrants_by_ethnicity)

# Calculate interaction effects - migrants in each ethnicity compared to White British non-migrant (WBNM) reference group

get_RR_and_95_CI <- function(dataframe) {
  output <- dataframe %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$Lower.CI, output$Upper.CI, sep ="-")
  output$RR <- paste0(output$Estimate, ' (', output$ci, ")")
  #output <- output %>%
  #  dplyr:: select(RR)
  output_name <- deparse(substitute(dataframe))
  assign(x = output_name, value = output, envir = globalenv())
}

## White British migrants 
white_migrant_vsWBNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_migrant_vsWBNM)
save(white_migrant_vsWBNM, file="filepath")
white_migrant <- white_migrant_vsWBNM %>%
  dplyr::select(RR)

# White non-British migrants 
white_nonbritish_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)White non-British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_migrant_vsWBNM)
save(white_nonbritish_migrant_vsWBNM, file="filepath")
white_nonbritish_migrant <- white_nonbritish_migrant_vsWBNM %>%
  dplyr::select(RR)

## Black migrants 
black_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_migrant_vsWBNM)
save(black_migrant_vsWBNM, file="filepath")
black_migrant <- black_migrant_vsWBNM %>%
  dplyr::select(RR)

## Asian
asian_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Asian/Asian British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_migrant_vsWBNM)
save(asian_migrant_vsWBNM, file="filepath")
asian_migrant <- asian_migrant_vsWBNM %>%
  dplyr::select(RR)

## Mixed
mixed_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_migrant_vsWBNM)
save(mixed_migrant_vsWBNM, file="filepath")
mixed_migrant <- mixed_migrant_vsWBNM %>%
  dplyr::select(RR)

## Other
other_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Other ethnic group' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_migrant_vsWBNM)
save(other_migrant_vsWBNM, file="filepath")
other_migrant <- other_migrant_vsWBNM %>%
  dplyr::select(RR)

## Unknown
unknown_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Unknown' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_migrant_vsWBNM)
save(unknown_migrant_vsWBNM, file="filepath")
unknown_migrant <- unknown_migrant_vsWBNM %>%
  dplyr::select(RR)

# Calculate interaction effects - migrants vs non-migrants in each ethnicity strata
## White
white_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_strata_vsNM)
save(white_strata_vsNM, file="filepath")
white_strata <- white_strata_vsNM %>%
  dplyr::select(RR)

## White non-British
white_nonbritish_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_strata_vsNM)
save(white_nonbritish_strata_vsNM, file="filepath")
white_nonbritish_strata <- white_nonbritish_strata_vsNM %>%
  dplyr::select(RR)

## Black
black_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_strata_vsNM)
save(black_strata_vsNM, file="filepath")
black_strata <- black_strata_vsNM %>%
  dplyr::select(RR)

## Asian
asian_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_strata_vsNM)
save(asian_strata_vsNM, file="filepath")
asian_strata <- asian_strata_vsNM %>%
  dplyr::select(RR)

## Mixed
mixed_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_strata_vsNM)
save(mixed_strata_vsNM, file="filepath")
mixed_strata <- mixed_strata_vsNM %>%
  dplyr::select(RR)

## Other
other_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_strata_vsNM)
save(other_strata_vsNM, file="filepath")
other_strata <- other_strata_vsNM %>%
  dplyr::select(RR)

## Unknown
unknown_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_strata_vsNM)
save(unknown_strata_vsNM, file="filepath")
unknown_strata <- unknown_strata_vsNM %>%
  dplyr::select(RR)

# Combine results to make table
Ethnicity <- data.frame(c('White British', 'White non-British', 'Black', 'Asian', 'Mixed', 'Other', 'Unknown'))
migrants_by_ethnicity <- bind_rows(white_migrant, white_nonbritish_migrant, black_migrant, asian_migrant, mixed_migrant, other_migrant, unknown_migrant)
Ethniciy_strata <- bind_rows(white_strata, white_nonbritish_strata, black_strata, asian_strata, mixed_strata, other_strata, unknown_strata)
matched_ethnicity_table <- bind_cols(Ethnicity, migrants_by_ethnicity, Ethniciy_strata)
row.names(ethnicity_table) <- NULL

matched_ethnicity_table <- matched_ethnicity_table %>%
  rename(
    Ethnicity = c..White.British....White.non.British....Black....Asian....Mixed...,
    RR_migrantsVsWBNM = RR...2,
    RR_MvsNM_each_ethnicity = RR...3
  ) %>%
  right_join(non_migrants_by_ethnicity, by = 'Ethnicity') %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, RR_migrantsVsWBNM, RR_MvsNM_each_ethnicity)

# Save dataframe for knitting in RMD file
save(matched_ethnicity_table, file="filepath")

### Additive and multiplicative interaction effects for each ethnicity 
# Based on Mathur MB & VanderWeele TJ (2018). R function for additive interaction measures. Epidemiology 29(1), e5-e6

# White non-British
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[3]
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_white_non_british <- data.frame(ethnicity = 'white_non_british', 
                                                    RERI = RERI,
                                                    RERI_upper_CI = upper_CI,
                                                    RERI_lower_CI = lower_CI,
                                                    multiplicative_effect = RR_interaction,
                                                    mult_lower_CI = beta_interaction_CI[1],
                                                    mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_white_non_british) <- NULL

# Black 
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[6]
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_black <- data.frame(ethnicity = 'black', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_black) <- NULL

# Asian
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[5] # Asian
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_asian <- data.frame(ethnicity = 'asian', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_asian) <- NULL

# Mixed
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[4] # Mixed
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_mixed <- data.frame(ethnicity = 'mixed', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_mixed) <- NULL

# Other
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[7] # Other
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_other <- data.frame(ethnicity = 'other', 
                                        RERI = RERI,
                                        RERI_upper_CI = upper_CI,
                                        RERI_lower_CI = lower_CI,
                                        multiplicative_effect = RR_interaction,
                                        mult_lower_CI = beta_interaction_CI[1],
                                        mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_other) <- NULL

# Unknown
# get coefficients and RRs
exposure1 <- names(coef(x))[2]
exposure2 <- names(coef(x))[8] # Unknown
interaction <- paste(exposure1, exposure2, sep = ':')

beta_migrant_status <- coef(x)[exposure1]
beta_ethnicity <- coef(x)[exposure2]
beta_interaction <- coef(x)[interaction]
beta_interaction_CI <- exp(confint.default(x)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_ethnicity <- exp(beta_ethnicity)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_ethnicity <- exp(beta_migrant_status + beta_ethnicity + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_ethnicity - RR_migrant_status - RR_ethnicity + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(x))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_ethnicity, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine
interaction_effects_unknown <- data.frame(ethnicity = 'unknown', 
                                          RERI = RERI,
                                          RERI_upper_CI = upper_CI,
                                          RERI_lower_CI = lower_CI,
                                          multiplicative_effect = RR_interaction,
                                          mult_lower_CI = beta_interaction_CI[1],
                                          mult_upper_CI = beta_interaction_CI[2])
row.names(interaction_effects_unknown) <- NULL

# Combine RERIs and save
int_effects_all_ethnicities <- bind_rows(interaction_effects_white_non_british, interaction_effects_black, 
                                         interaction_effects_asian, interaction_effects_mixed, 
                                         interaction_effects_other, interaction_effects_unknown)

write_csv(int_effects_all_ethnicities, "filepath")

## 11_London only  ----

## Multivariable model by migrant_status  ---

exact_match_London_annual_conscounts <- matched_cohort_England_2015_2020_annual_conscounts %>%
  filter(prac_region == "London")

## All patients (i.e. all ages)

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = exact_match_London_annual_conscounts)
extract_glm_results_allages(x, 'London_glm_mig', exact_match_London_annual_conscounts, migrant_status)

write_csv(London_glm_mig_table, "filepath")


# IRR (multivariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts)
extract_glm_results_allages(x, 'London_multivariable_allages', exact_match_London_annual_conscounts, migrant_status)

## Age_subcohorts

# IRR (multivariable glm)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_0to15", "0-15 years", exact_match_London_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_16to24", "16-24 years", exact_match_London_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_25to34", "25-34 years", exact_match_London_annual_conscounts, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_35to49", "35-49 years", exact_match_London_annual_conscounts, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_50to64", "50-64 years", exact_match_London_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) , 
            data = exact_match_London_annual_conscounts[exact_match_London_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_65plus", ">=65 years", exact_match_London_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
London_multivariable_migrant_status <- bind_rows(London_multivariable_allages_table, London_glm_multivariable_0to15_table,London_glm_multivariable_16to24_table,
                                                 London_glm_multivariable_25to34_table,London_glm_multivariable_35to49_table,London_glm_multivariable_50to64_table,London_glm_multivariable_65plus_table)


write_csv(London_multivariable_migrant_status, "filepath")

## 12_Study year  ----

## Multivariable IR and IRRs 

## 2015

studyyear_2015 <- matched_cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2015)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015)
extract_glm_results_allages(x, 'multivariable_2015_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2015[studyyear_2015$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2015_migrant_status <- bind_rows(multivariable_2015_allages_table, glm_2015_multivariable_0to15_table,glm_2015_multivariable_16to24_table,
                                               glm_2015_multivariable_25to34_table,glm_2015_multivariable_35to49_table,glm_2015_multivariable_50to64_table,glm_2015_multivariable_65plus_table)

write_csv(multivariable_2015_migrant_status, "filepath")

## 2016

studyyear_2016 <- matched_cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2016)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016)
extract_glm_results_allages(x, 'multivariable_2016_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2016[studyyear_2016$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2016_migrant_status <- bind_rows(multivariable_2016_allages_table, glm_2016_multivariable_0to15_table,glm_2016_multivariable_16to24_table,
                                               glm_2016_multivariable_25to34_table,glm_2016_multivariable_35to49_table,glm_2016_multivariable_50to64_table,glm_2016_multivariable_65plus_table)


write_csv(multivariable_2016_migrant_status, "filepath")

## 2017

studyyear_2017 <- matched_cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2017)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017)
extract_glm_results_allages(x, 'multivariable_2017_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2017[studyyear_2017$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2017_migrant_status <- bind_rows(multivariable_2017_allages_table, glm_2017_multivariable_0to15_table,glm_2017_multivariable_16to24_table,
                                               glm_2017_multivariable_25to34_table,glm_2017_multivariable_35to49_table,glm_2017_multivariable_50to64_table,glm_2017_multivariable_65plus_table)


write_csv(multivariable_2017_migrant_status, "filepath")

## 2018

studyyear_2018 <- matched_cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2018)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018)
extract_glm_results_allages(x, 'multivariable_2018_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2018[studyyear_2018$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2018_migrant_status <- bind_rows(multivariable_2018_allages_table, glm_2018_multivariable_0to15_table,glm_2018_multivariable_16to24_table,
                                               glm_2018_multivariable_25to34_table,glm_2018_multivariable_35to49_table,glm_2018_multivariable_50to64_table,glm_2018_multivariable_65plus_table)


write_csv(multivariable_2018_migrant_status, "filepath")

## 2019

studyyear_2019 <- matched_cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2019)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019)
extract_glm_results_allages(x, 'multivariable_2019_allages', matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_0to15", "0-15 years", matched_cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_16to24", "16-24 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_25to34", "25-34 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_35to49", "35-49 years", matched_cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_50to64", "50-64 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) +  
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = studyyear_2019[studyyear_2019$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_65plus", ">=65 years", matched_cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_2019_migrant_status <- bind_rows(multivariable_2019_allages_table, glm_2019_multivariable_0to15_table,glm_2019_multivariable_16to24_table,
                                               glm_2019_multivariable_25to34_table,glm_2019_multivariable_35to49_table,glm_2019_multivariable_50to64_table,glm_2019_multivariable_65plus_table)


write_csv(multivariable_2019_migrant_status, "filepath")
