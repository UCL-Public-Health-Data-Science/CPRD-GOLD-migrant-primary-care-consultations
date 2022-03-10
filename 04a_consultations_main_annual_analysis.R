## 0_Description ---------------------------------------------------------------------------

# Migrants' primary care utilisation before and during the COVID-19 pandemic in England: An interrupted time series
# Annual pre-pandemic main analysis 
# Date started: 30/11/2020
# Author(s): Claire Zhang / Yamina Boukari / Neha Pathak 
# QA (date): Yamina Boukari (04/03/2022)

## 1_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc, epitools, data.table, epiR, MASS, psych, forestplot, dsr, gmodels, Epi, msm, car)

## 2_Set working directory ------------------------------------------------------------------

setwd("S:/CALIBER_19_062R/01_gold/01_all_patients")

## 3_Load datasets and functions ------------------------------------------------------

# Datasets 

load(file = "primary_care_cleaned_files/02_Cons_Rdata/cohort_England_2015_2020_annual_conscounts.Rdata")

# Functions

function_to_calculate_IR <- function(dataset, output_name) {
  output <- pois.exact(dataset$conscount, dataset$pyears, conf.level=0.95)
  output <- output %>%
    dplyr::mutate(incidence_rate1=rate*1)
  output <- output %>%
    dplyr::mutate(lower1=lower*1)
  output <- output %>%
    dplyr::mutate(upper1=upper*1)
  output <- output %>% 
    dplyr::rename(incidence_rate = rate)
  output <- output[, 3:9]
  output <- bind_cols(dataset, output) 
  output <- output %>% dplyr::mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower1, output$upper1, sep ="-")
  output$ir_ci <- paste0(output$incidence_rate, ' (', output$ci, ')')
  output <-  dplyr::select(output, -c(incidence_rate, lower, upper, conf.level, incidence_rate1, lower1,upper1,ci))
  output <- output %>% 
    dplyr::rename(events = conscount) %>%
    dplyr::rename(person_years = pyears)
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
  output <- tibble::as_tibble(output)
  output$estimate <- as.numeric(output$estimate)
  output <- output %>% dplyr::rename(lower = "2.5 %") %>% dplyr::rename(upper = "97.5 %") 
  output$upper <- as.numeric(output$upper)
  output$lower <- as.numeric(output$lower)
  output$p <- as.numeric(output$p)
  output <- output %>% dplyr::mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$estimate[1] <- 1.00
  output$lower[1] <- 1.00
  output$upper[1] <- 1.00
  output$ci[1] <- "1.00-1.00"
  output$irr_ci <- paste(output$estimate, output$ci, sep =",")
  output_table <- dplyr::select(output,names, irr_ci )
  output_table <- output_table %>%
    dplyr::mutate(age_subcohort = 'All_ages')
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
  output <- tibble::as_tibble(output)
  output$estimate <- as.numeric(output$estimate)
  output <- output %>% dplyr::rename(lower = "2.5 %") %>% dplyr::rename(upper = "97.5 %") 
  output$upper <- as.numeric(output$upper)
  output$lower <- as.numeric(output$lower)
  output$p <- as.numeric(output$p)
  output <- output %>% dplyr::mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$estimate[1] <- 1.00
  output$lower[1] <- 1.00
  output$upper[1] <- 1.00
  output$ci[1] <- "1.00-1.00"
  output$irr_ci <- paste(output$estimate, output$ci, sep =",")
  output <- output %>%
    dplyr::mutate(age_subcohort = subcohort)
  output_table <- dplyr::select(output, names, irr_ci, age_subcohort)
  assign(x = output_name, value = output, envir = globalenv()) 
  assign(x = paste0(output_name,'_table'), value = output_table, envir = globalenv()) 
}

## 4_Select needed variables -------------------------------------------------

cohort_England_2015_2020_annual_conscounts <- dplyr::select(cohort_England_2015_2020_annual_conscounts, 
                                                          c(patid, pyears, conscount, migrant_status,
                                                            age_subcohort, imd, ethnicat6, gender,
                                                            studyyear, studyyear_agecat, prac_region,
                                                            migcertainty, cohort_entry))

# Relevel
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, gender <- relevel (gender, ref="Male")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, imd <- relevel (imd, ref="IMD 1"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, prac_region <- relevel (prac_region, ref="London")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, ethnicat6 <- relevel (ethnicat6, ref="White British"))

########################### MAIN ANALYSIS####################################### -----

## 5_Summary outcome measures ------------------------------------------------------------------

### Overall (i.e. full study period)

## All patients

# All patients
conscount_summary_overall_allpatients <- cohort_England_2015_2020_annual_conscounts %>% 
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migrant status
conscount_summary_mig_status_allpatients <- cohort_England_2015_2020_annual_conscounts %>% group_by(migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# All patients by migcertainty
conscount_summary_migcertainty_allpatients <- cohort_England_2015_2020_annual_conscounts %>% group_by(migcertainty) %>%
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
conscount_summary_overall_agesubcohorts <- cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migrant status
conscount_summary_mig_status_agesubcohorts <- cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort, migrant_status) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount),n_individuals = n_distinct(patid))

# Age_subcohorts by migcertainty
conscount_summary_migcertainty_agesubcohorts <- cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort,migcertainty) %>%
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
conscount_summary_fullperiod <- bind_rows(conscount_summary_overall_1,conscount_summary_overall_2)

# save file
write_csv(conscount_summary_fullperiod, "results/01_Consultations/full_cohort/conscount_summary_fullperiod.csv" )

### Annual

## All patients

# All patients
conscount_annual_allpatients <- cohort_England_2015_2020_annual_conscounts %>% group_by(studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migrant status
conscount_annual_mig_status_allpatients <- cohort_England_2015_2020_annual_conscounts %>% group_by(migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

# All patients by migcertainty 
conscount_annual_migcertainty_allpatients <- cohort_England_2015_2020_annual_conscounts %>% group_by(migcertainty, studyyear) %>%
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
conscount_annual_agesubcohorts <- cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migrant status
conscount_annual_mig_status_agesubcohorts <- cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort,migrant_status, studyyear) %>%
  summarise(n_events = sum(conscount), mean = mean(conscount), sd = sd(conscount), median = median(conscount), iqr = IQR(conscount), 
            min = min(conscount), max = max(conscount))

## age_subcohorts by migcertainty 
conscount_annual_migcertainty_agesubcohorts <- cohort_England_2015_2020_annual_conscounts %>% group_by(age_subcohort, migcertainty, studyyear) %>%
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
conscount_summary_annual <- bind_rows(conscount_summary_annual_1,conscount_summary_annual_2)

## Export  summary measures
write_csv(conscount_summary_annual, "results/01_Consultations/full_cohort/conscount_summary_annual.csv" )


## 6_IRs by migrant_status + studyyear -----

# IRs by migrant_status and studyyear

## All ages

pyears_allpatients_overall_migstatus_studyyear <- aggregate(pyears ~ migrant_status + studyyear, cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migstatus_studyyear <- pyears_allpatients_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus_studyyear <- aggregate(conscount ~ migrant_status + studyyear, cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migstatus_studyyear <- inner_join(conscount_allpatients_overall_migstatus_studyyear,pyears_allpatients_overall_migstatus_studyyear, by= c("migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migstatus_studyyear, output_name = 'IR_allpatients_overall_migstatus_studyyear')
# Save for tables
save(IR_allpatients_overall_migstatus_studyyear, file = 'S:/CALIBER_19_062R/01_gold/01_all_patients/results/01_Consultations/tables/IR_allpatients_overall_migstatus_studyyear.Rdata')
IR_allpatients_overall_migstatus_studyyear$age_subcohort <- 'All_ages'

## age_subcohorts

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus_studyyear <- aggregate(pyears ~ age_subcohort + migrant_status + studyyear, cohort_England_2015_2020_annual_conscounts, sum) 
pyears_agesubcohorts_overall_migstatus_studyyear <- pyears_agesubcohorts_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus_studyyear <- aggregate(conscount ~ age_subcohort+ migrant_status + studyyear, cohort_England_2015_2020_annual_conscounts, sum) 
IR_agesubcohorts_overall_migstatus_studyyear <- inner_join(conscount_agesubcohorts_overall_migstatus_studyyear,pyears_agesubcohorts_overall_migstatus_studyyear, 
                                                           by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus_studyyear, output_name = 'IR_agesubcohorts_overall_migstatus_studyyear')

## Combine

IR_migstatus_studyyear <- full_join(IR_allpatients_overall_migstatus_studyyear, IR_agesubcohorts_overall_migstatus_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migrant_status' = 'migrant_status',
                                                                                                                                     'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))

## Reorder results
IR_migstatus_studyyear <- arrange(IR_migstatus_studyyear, migrant_status, age_subcohort) %>%
  relocate(age_subcohort)

write_csv(IR_migstatus_studyyear, "results/01_Consultations/full_cohort/IR_migstatus_studyyear.csv")

# IMPORTANT STEP: Filter out 2020 for modelling ------
cohort_England_2015_2020_annual_conscounts <- cohort_England_2015_2020_annual_conscounts %>%
  filter(studyyear != 2020)

## 7_IR & univariable IRR by migrant status ----------------------------------------------------------------------

## All patients (i.e. all ages)

# IR by migrant status (i.e. not including 2020, which was filtered out in the step above)
pyears_allpatients_overall_migstatus <- aggregate(pyears ~ migrant_status, cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migstatus <- pyears_allpatients_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus <- aggregate(conscount ~ migrant_status, cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migstatus <- inner_join(conscount_allpatients_overall_migstatus,pyears_allpatients_overall_migstatus, by= c("migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_allpatients_overall_migstatus, output_name = 'IR_allpatients_overall_migstatus')
save(IR_allpatients_overall_migstatus, file = 'S:/CALIBER_19_062R/01_gold/01_all_patients/results/01_Consultations/tables/IR_allpatients_overall_migstatus.Rdata')
IR_allpatients_overall_migstatus$age_subcohort <- 'All_ages'

# Relevel
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, migrant_status <- relevel (migrant_status, ref="Non-migrant") )
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, age_subcohort <- relevel (age_subcohort, ref="0-15 years") )

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = cohort_England_2015_2020_annual_conscounts)

extract_glm_results_allages(x, 'glm_mig', cohort_England_2015_2020_annual_conscounts, migrant_status)
round(ci.lin(x,Exp=T),3)

save(glm_mig, file="results/01_Consultations/full_cohort/glm_mig.Rdata")

# Join glm + IR + univariable_mig 
univariable_migrant_status <- full_join(IR_allpatients_overall_migstatus, glm_mig_table, by = c("age_subcohort" = "age_subcohort","migrant_status" = "names")) %>%
  relocate(age_subcohort)

## Age_subcohorts - not included in paper 

# IR by migrant_status
pyears_agesubcohorts_overall_migstatus <- aggregate(pyears ~ age_subcohort + migrant_status, cohort_England_2015_2020_annual_conscounts, sum) 
pyears_agesubcohorts_overall_migstatus <- pyears_agesubcohorts_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migstatus <- aggregate(conscount ~ age_subcohort+ migrant_status, cohort_England_2015_2020_annual_conscounts, sum) 
IR_agesubcohorts_overall_migstatus <- inner_join(conscount_agesubcohorts_overall_migstatus,pyears_agesubcohorts_overall_migstatus, 
                                                 by= c("age_subcohort" = 'age_subcohort',"migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_agesubcohorts_overall_migstatus, output_name = 'IR_agesubcohorts_overall_migstatus')

# IRR by migrant status (univariable glm) 

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_50to64", "50-64 years",cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over  
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_mig_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## NOTE IRR is comparing migrants at X age to non-migrants of the SAME age 

# Join age_subcohorts glm+IR and all_ages glm+IR 
glm_mig_combined_table <- bind_rows(glm_mig_0to15_table,glm_mig_16to24_table,glm_mig_25to34_table,glm_mig_35to49_table,glm_mig_50to64_table,glm_mig_65plus_table)

glm_mig_combined_table$names <- as.factor(glm_mig_combined_table$names)

univariable_migrant_status_agesubcohorts <- left_join(IR_agesubcohorts_overall_migstatus, glm_mig_combined_table, 
                                                      by = c("migrant_status" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migrant_status <- bind_rows(univariable_migrant_status, univariable_migrant_status_agesubcohorts)  
univariable_migrant_status$age_subcohort <- as.factor(univariable_migrant_status$age_subcohort)
univariable_migrant_status$age_subcohort <- factor(univariable_migrant_status$age_subcohort, 
                                                   levels = c('All_ages', '0-15 years','16-24 years', '25-34 years','35-49 years', '50-64 years', '>=65 years'))
## Reorder results (can change depending on how we want to present it)
univariable_migrant_status <- arrange(univariable_migrant_status, migrant_status, age_subcohort)

write_csv(univariable_migrant_status, "results/01_Consultations/full_cohort/univariable_migrant_status.csv")


## 8a_Multivariable model by migrant_status with IMD adjustment  ----

## All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts)

extract_glm_results_allages(x, 'multivariable_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)
round(ci.lin(x,Exp=T),3)

# Check the variance inflation factor (to check multicolinearity)
car::vif(x)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)
round(ci.lin(x,Exp=T),3)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)
round(ci.lin(x,Exp=T),3)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)
round(ci.lin(x,Exp=T),3)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)
round(ci.lin(x,Exp=T),3)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migrant_status)
round(ci.lin(x,Exp=T),3)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)
round(ci.lin(x,Exp=T),3)

# Join all ages and age_subcohorts glm results
multivariable_migrant_status <- dplyr::bind_rows(multivariable_allages_table, glm_multivariable_0to15_table,glm_multivariable_16to24_table,
                                                 glm_multivariable_25to34_table,glm_multivariable_35to49_table,glm_multivariable_50to64_table,glm_multivariable_65plus_table)

write_csv(multivariable_migrant_status, "results/01_Consultations/full_cohort/multivariable_migrant_status.csv")

# Save .Rdata for forestplots
save(multivariable_allages, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs/multivariable_allages.Rdata")
save(glm_multivariable_0to15, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs/glm_multivariable_0to15.Rdata")
save(glm_multivariable_16to24, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs/glm_multivariable_16to24.Rdata")
save(glm_multivariable_25to34, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs/glm_multivariable_25to34.Rdata")
save(glm_multivariable_35to49, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs/glm_multivariable_35to49.Rdata")
save(glm_multivariable_50to64, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs/glm_multivariable_50to64.Rdata")
save(glm_multivariable_65plus, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs/glm_multivariable_65plus.Rdata")


## 8b_Multivariable model by migrant_status without IMD adjustment  ----

## All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts)

extract_glm_results_allages(x, 'multivariable_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)
round(ci.lin(x,Exp=T),3)

## Age_subcohorts (not included in paper)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) + 
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(prac_region) + offset(log(pyears)) , 
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

# Join all ages and age_subcohorts glm results
multivariable_migrant_status <- dplyr::bind_rows(multivariable_allages_table, glm_multivariable_0to15_table,glm_multivariable_16to24_table,
                                                 glm_multivariable_25to34_table,glm_multivariable_35to49_table,glm_multivariable_50to64_table,glm_multivariable_65plus_table)

multivariable_migrant_status_noIMD <- multivariable_migrant_status
write_csv(multivariable_migrant_status, "results/01_Consultations/full_cohort/multivariable_migrant_status_noIMD.csv")

# Save .Rdata for forestplots
save(multivariable_allages, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs_noIMD/multivariable_allages.Rdata")
save(glm_multivariable_0to15, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs_noIMD/glm_multivariable_0to15.Rdata")
save(glm_multivariable_16to24, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs_noIMD/glm_multivariable_16to24.Rdata")
save(glm_multivariable_25to34, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs_noIMD/glm_multivariable_25to34.Rdata")
save(glm_multivariable_35to49, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs_noIMD/glm_multivariable_35to49.Rdata")
save(glm_multivariable_50to64, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs_noIMD/glm_multivariable_50to64.Rdata")
save(glm_multivariable_65plus, file="results/01_Consultations/full_cohort/multivariable_migrant_status_forestplot_inputs_noIMD/glm_multivariable_65plus.Rdata")


## 9_London only  ----

## Multivariable model by migrant_status ---

London_annual_conscounts_2015_2020 <- cohort_England_2015_2020_annual_conscounts %>%
  filter(prac_region == "London")

## All patients (i.e. all ages)

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migrant_status) + offset(log(pyears)), data = London_annual_conscounts_2015_2020)
extract_glm_results_allages(x, 'London_glm_mig', London_annual_conscounts_2015_2020, migrant_status)
round(ci.lin(x,Exp=T),3)

write_csv(London_glm_mig_table, "results/01_Consultations/full_cohort/London_univariable_migrant_status.csv")

# IRR (multivariable glm) without IMD
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020)
extract_glm_results_allages(x, 'London_multivariable_allages_noIMD', London_annual_conscounts_2015_2020, migrant_status)
round(ci.lin(x,Exp=T),3)

# IRR (multivariable glm) with IMD
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020)
extract_glm_results_allages(x, 'London_multivariable_allages', London_annual_conscounts_2015_2020, migrant_status)
round(ci.lin(x,Exp=T),3)

## Age_subcohorts

# IRR (multivariable glm)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020[London_annual_conscounts_2015_2020$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_0to15", "0-15 years", London_annual_conscounts_2015_2020, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020[London_annual_conscounts_2015_2020$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_16to24", "16-24 years", London_annual_conscounts_2015_2020, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020[London_annual_conscounts_2015_2020$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_25to34", "25-34 years", London_annual_conscounts_2015_2020, migrant_status)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020[London_annual_conscounts_2015_2020$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_35to49", "35-49 years", London_annual_conscounts_2015_2020, migrant_status)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020[London_annual_conscounts_2015_2020$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_50to64", "50-64 years", London_annual_conscounts_2015_2020, migrant_status)

## 65 years and over
x <- glm.nb(conscount ~ as.factor(migrant_status) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + offset(log(pyears)) ,
            data = London_annual_conscounts_2015_2020[London_annual_conscounts_2015_2020$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "London_glm_multivariable_65plus", ">=65 years", London_annual_conscounts_2015_2020, migrant_status)

# Join all ages and age_subcohorts glm results
London_multivariable_migrant_status <- bind_rows(London_multivariable_allages_table, London_glm_multivariable_0to15_table,London_glm_multivariable_16to24_table,
                                                 London_glm_multivariable_25to34_table,London_glm_multivariable_35to49_table, London_glm_multivariable_50to64_table,London_glm_multivariable_65plus_table)

write_csv(London_multivariable_migrant_status, "results/01_Consultations/full_cohort/London_multivariable_migrant_status.csv")

# Save .Rdata for forestplots

save(London_glm_mig, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_glm_mig.Rdata")
save(London_multivariable_allages, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_multivariable_allages.Rdata")
save(London_multivariable_allages_noIMD, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs_noIMD/London_multivariable_allages_noIMD.Rdata")
save(London_glm_multivariable_0to15, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_glm_multivariable_0to15.Rdata")
save(London_glm_multivariable_16to24, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_glm_multivariable_16to24.Rdata")
save(London_glm_multivariable_25to34, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_glm_multivariable_25to34.Rdata")
save(London_glm_multivariable_35to49, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_glm_multivariable_35to49.Rdata")
save(London_glm_multivariable_50to64, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_glm_multivariable_50to64.Rdata")
save(London_glm_multivariable_65plus, file="results/01_Consultations/full_cohort/London_multivariable_migrant_status_forestplot_inputs/London_glm_multivariable_65plus.Rdata")


## 10_Ethnicity (ethnicat6) - interaction term, additive and multiplicative effects ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status)*as.factor(ethnicat6) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'multivariable_ethnicat6interaction_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)
round(ci.lin(x,Exp=T),3)

# Save .Rdata for forestplots
save(multivariable_ethnicat6interaction_allages, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/multivariable_ethnicat6interaction_allages.Rdata")
write_csv(multivariable_ethnicat6interaction_allages, "results/01_Consultations/full_cohort/multivariable_ethnicat6interaction_allages.csv")

# RR (95% CI) for non-migrants in each ethnicity compared to White British non-migrants
# For table
white_non_migrant <- data.frame(Ethnicity = 'White British', RR_nonmigrantsVsWBNM = '1.0')
non_migrants_by_ethnicity <- multivariable_ethnicat6interaction_allages[3:8,] %>%
  mutate(ci = paste(lower, upper, sep ="-"),
         RR_nonmigrantsVsWBNM = paste0(estimate, ' (', ci, ")")) %>%
  mutate(Ethnicity = c('White non-British', 'Mixed', 'Asian', 'Black', 'Other', 'Unknown')) %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, estimate, upper, lower)
non_migrants_by_ethnicity <- bind_rows(white_non_migrant, non_migrants_by_ethnicity)
# Add row for All (to match table columns further down)
all_non_migrants <- data.frame(Ethnicity = 'All', RR_nonmigrantsVsWBNM = '-')
non_migrants_by_ethnicity <- bind_rows(non_migrants_by_ethnicity, all_non_migrants)

# Calculate interaction effects - migrants in each ethnicity compared to White British non-migrant (WBNM) reference group

get_RR_and_95_CI <- function(dataframe) {
  output <- dataframe %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$Lower.CI, output$Upper.CI, sep ="-")
  output$RR <- paste0(output$Estimate, ' (', output$ci, ")")
  output_name <- deparse(substitute(dataframe))
  assign(x = output_name, value = output, envir = globalenv())
}

## White British migrants 
white_migrant_vsWBNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_migrant_vsWBNM)
save(white_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/white_migrant_vsWBNM.Rdata")
white_migrant <- white_migrant_vsWBNM %>%
  dplyr::select(RR)

# White non-British migrants 
white_nonbritish_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)White non-British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_migrant_vsWBNM)
save(white_nonbritish_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/white_nonbritish_migrant_vsWBNM.Rdata")
white_nonbritish_migrant <- white_nonbritish_migrant_vsWBNM %>%
  dplyr::select(RR)

## Asian
asian_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Asian/Asian British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_migrant_vsWBNM)
save(asian_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/asian_migrant_vsWBNM.Rdata")
asian_migrant <- asian_migrant_vsWBNM %>%
  dplyr::select(RR)

## Black migrants 
black_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_migrant_vsWBNM)
save(black_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/black_migrant_vsWBNM.Rdata")
black_migrant <- black_migrant_vsWBNM %>%
  dplyr::select(RR)

## Mixed
mixed_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_migrant_vsWBNM)
save(mixed_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/mixed_migrant_vsWBNM.Rdata")
mixed_migrant <- mixed_migrant_vsWBNM %>%
  dplyr::select(RR)

## Other
other_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Other ethnic group' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_migrant_vsWBNM)
save(other_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/other_migrant_vsWBNM.Rdata")
other_migrant <- other_migrant_vsWBNM %>%
  dplyr::select(RR)

## Unknown
unknown_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)Unknown' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_migrant_vsWBNM)
save(unknown_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/unknown_migrant_vsWBNM.Rdata")
unknown_migrant <- unknown_migrant_vsWBNM %>%
  dplyr::select(RR)

# All 
all_migrant_vsWBNM <- exp(estimable(x, c('as.factor(ethnicat6)White non-British' = 1,
                                         'as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1,
                                         'as.factor(ethnicat6)Asian/Asian British' = 1,
                                         'as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1, 
                                         'as.factor(ethnicat6)Other ethnic group' = 1, 
                                         'as.factor(ethnicat6)Unknown' = 1, 
                                         'as.factor(migrant_status)Migrant' = 1,
                                         'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1,
                                         'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1,
                                         'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1,
                                         'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1,
                                         'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1,
                                         'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(all_migrant_vsWBNM)
save(all_migrant_vsWBNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/all_migrant_vsWBNM.Rdata")
all_migrant <- all_migrant_vsWBNM %>%
  dplyr::select(RR)

# Calculate interaction effects - migrants vs non-migrants in each ethnicity strata
## White
white_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(white_strata_vsNM)
save(white_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/white_strata_vsNM.Rdata")
white_strata <- white_strata_vsNM %>%
  dplyr::select(RR)

## White non-British
white_nonbritish_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)White non-British' = 1), conf=.95))
get_RR_and_95_CI(white_nonbritish_strata_vsNM)
save(white_nonbritish_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/white_nonbritish_strata_vsNM.Rdata")
white_nonbritish_strata <- white_nonbritish_strata_vsNM %>%
  dplyr::select(RR)

## Black
black_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Black/African/Caribbean/Black British' = 1), conf=.95))
get_RR_and_95_CI(black_strata_vsNM)
save(black_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/black_strata_vsNM.Rdata")
black_strata <- black_strata_vsNM %>%
  dplyr::select(RR)

## Asian
asian_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Asian/Asian British' = 1), conf=.95))
get_RR_and_95_CI(asian_strata_vsNM)
save(asian_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/asian_strata_vsNM.Rdata")
asian_strata <- asian_strata_vsNM %>%
  dplyr::select(RR)

## Mixed
mixed_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Mixed/Multiple ethnic groups' = 1), conf=.95))
get_RR_and_95_CI(mixed_strata_vsNM)
save(mixed_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/mixed_strata_vsNM.Rdata")
mixed_strata <- mixed_strata_vsNM %>%
  dplyr::select(RR)

## Other
other_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_strata_vsNM)
save(other_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/other_strata_vsNM.Rdata")
other_strata <- other_strata_vsNM %>%
  dplyr::select(RR)

## Unknown
unknown_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat6)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_strata_vsNM)
save(unknown_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/unknown_strata_vsNM.Rdata")
unknown_strata <- unknown_strata_vsNM %>%
  dplyr::select(RR)

# Empty row to match the vsWBNM table
all_strata <- data.frame(RR = '-')

# Combine results to make table
Ethnicity <- data.frame(c('White British', 'White non-British', 'Mixed', 'Asian', 'Black', 'Other', 'Unknown', 'All'))
migrants_by_ethnicity <- bind_rows(white_migrant, white_nonbritish_migrant, mixed_migrant, asian_migrant, black_migrant, other_migrant, unknown_migrant, all_migrant)
Ethniciy_strata <- bind_rows(white_strata, white_nonbritish_strata, mixed_strata, asian_strata, black_strata, other_strata, unknown_strata, all_strata)
ethnicity_table <- bind_cols(Ethnicity, migrants_by_ethnicity, Ethniciy_strata)
row.names(ethnicity_table) <- NULL

ethnicity_table <- ethnicity_table %>%
  rename(
    Ethnicity = c..White.British....White.non.British....Mixed....Asian....Black...,
    RR_migrantsVsWBNM = RR...2,
    RR_MvsNM_each_ethnicity = RR...3
  ) %>%
  right_join(non_migrants_by_ethnicity, by = 'Ethnicity') %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, RR_migrantsVsWBNM, RR_MvsNM_each_ethnicity)

# Save dataframe for knitting in RMD file
save(ethnicity_table, file="results/01_Consultations/full_cohort/ethnicity_interactions_table.Rdata")

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
int_effects_all_ethnicities <- bind_rows(interaction_effects_white_non_british, interaction_effects_mixed,
                                         interaction_effects_asian, interaction_effects_black,
                                         interaction_effects_other, interaction_effects_unknown)

write_csv(int_effects_all_ethnicities, "results/01_Consultations/full_cohort/int_effects_all_ethnicities.csv")
save(int_effects_all_ethnicities, file = "results/01_Consultations/tables/int_effects_all_ethnicities.Rdata")

########################### SENSITIVITY ANALYSES ############################### ----

## 11_Certainty of migration status ----

# Reload to include 2020
load(file = "primary_care_cleaned_files/02_Cons_Rdata/cohort_England_2015_2020_annual_conscounts.Rdata")

cohort_England_2015_2020_annual_conscounts <- dplyr::select(cohort_England_2015_2020_annual_conscounts, 
                                                            c(patid, pyears, conscount, migrant_status,
                                                              age_subcohort, imd, ethnicat6, gender,
                                                              studyyear, studyyear_agecat, prac_region,
                                                              migcertainty, cohort_entry))

# Relevel
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, gender <- relevel (gender, ref="Male")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, imd <- relevel (imd, ref="IMD 1"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, prac_region <- relevel (prac_region, ref="London")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, ethnicat6 <- relevel (ethnicat6, ref="White British"))

# IRs by migcertainty + studyyear

## All ages

pyears_allpatients_overall_migcertainty_studyyear <- aggregate(pyears ~ migcertainty + studyyear, cohort_England_2015_2020_annual_conscounts, sum)
pyears_allpatients_overall_migcertainty_studyyear <- pyears_allpatients_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty_studyyear <- aggregate(conscount ~ migcertainty + studyyear, cohort_England_2015_2020_annual_conscounts, sum)
IR_allpatients_overall_migcertainty_studyyear <- inner_join(conscount_allpatients_overall_migcertainty_studyyear,pyears_allpatients_overall_migcertainty_studyyear, by= c("migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migcertainty_studyyear, output_name = 'IR_allpatients_overall_migcertainty_studyyear')
save(IR_allpatients_overall_migcertainty_studyyear, file = 'S:/CALIBER_19_062R/01_gold/01_all_patients/results/01_Consultations/tables/IR_allpatients_overall_migcertainty_studyyear.Rdata')
IR_allpatients_overall_migcertainty_studyyear$age_subcohort <- 'All_ages'

## Age_subcohorts

# IR by migcertainty
pyears_agesubcohorts_overall_migcertainty_studyyear <- aggregate(pyears ~ age_subcohort + migcertainty + studyyear, cohort_England_2015_2020_annual_conscounts, sum)
pyears_agesubcohorts_overall_migcertainty_studyyear <- pyears_agesubcohorts_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty_studyyear <- aggregate(conscount ~ age_subcohort+ migcertainty + studyyear, cohort_England_2015_2020_annual_conscounts, sum)
IR_agesubcohorts_overall_migcertainty_studyyear <- inner_join(conscount_agesubcohorts_overall_migcertainty_studyyear,pyears_agesubcohorts_overall_migcertainty_studyyear,
                                                              by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty_studyyear, output_name = 'IR_agesubcohorts_overall_migcertainty_studyyear')

## Combine

IR_migcertainty_studyyear <- dplyr::full_join(IR_allpatients_overall_migcertainty_studyyear, IR_agesubcohorts_overall_migcertainty_studyyear, by = c('age_subcohort' = 'age_subcohort', 'migcertainty' = 'migcertainty',
                                                                                                                                                     'studyyear' = 'studyyear', 'events' = 'events', 'person_years' = 'person_years', 'ir_ci' = 'ir_ci'))
## Reorder results (can change depending on how we want to present it)
IR_migcertainty_studyyear <- arrange(IR_migcertainty_studyyear, migcertainty, age_subcohort) %>%
  relocate(age_subcohort)

# IR_migcertainty_studyyear <- arrange(IR_migcertainty_studyyear, studyyear) %>%
#   relocate(age_subcohort)


write_csv(IR_migcertainty_studyyear, "results/01_Consultations/full_cohort/IR_migcertainty_studyyear.csv")

# Filter out 2020 for modelling
cohort_England_2015_2020_annual_conscounts <- cohort_England_2015_2020_annual_conscounts %>%
  filter(studyyear != 2020)

## Univariable model

## All patients (i.e. all ages)

# IR
pyears_allpatients_overall_migcertainty <- aggregate(pyears ~ migcertainty, cohort_England_2015_2020_annual_conscounts, sum)
pyears_allpatients_overall_migcertainty <- pyears_allpatients_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty <- aggregate(conscount ~ migcertainty, cohort_England_2015_2020_annual_conscounts, sum)
IR_allpatients_overall_migcertainty <- inner_join(conscount_allpatients_overall_migcertainty,pyears_allpatients_overall_migcertainty, by= c("migcertainty" = "migcertainty"))

function_to_calculate_IR(IR_allpatients_overall_migcertainty, output_name = 'IR_allpatients_overall_migcertainty')
save(IR_allpatients_overall_migcertainty, file = 'S:/CALIBER_19_062R/01_gold/01_all_patients/results/01_Consultations/tables/IR_allpatients_overall_migcertainty.Rdata')
IR_allpatients_overall_migcertainty$age_subcohort <- 'All_ages'

# IRR (univariable glm)
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)), data = cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'glm_migcertainty', cohort_England_2015_2020_annual_conscounts, migcertainty)
round(ci.lin(x,Exp=T),3)

# Save for forestplot
save(glm_migcertainty, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/glm_migcertainty.Rdata")

# Join glm + IR + univariable_mig
univariable_migcertainty <- dplyr::full_join(IR_allpatients_overall_migcertainty, glm_migcertainty_table, by = c("age_subcohort" = "age_subcohort","migcertainty" = "names")) %>%
  relocate(age_subcohort)

## Age_subcohorts (not included in paper)

# IR
pyears_agesubcohorts_overall_migcertainty <- aggregate(pyears ~ age_subcohort + migcertainty, cohort_England_2015_2020_annual_conscounts, sum)
pyears_agesubcohorts_overall_migcertainty <- pyears_agesubcohorts_overall_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_agesubcohorts_overall_migcertainty <- aggregate(conscount ~ age_subcohort+ migcertainty, cohort_England_2015_2020_annual_conscounts, sum)
IR_agesubcohorts_overall_migcertainty <- inner_join(conscount_agesubcohorts_overall_migcertainty,pyears_agesubcohorts_overall_migcertainty,
                                                    by= c("age_subcohort" = 'age_subcohort',"migcertainty" = "migcertainty"))
function_to_calculate_IR(IR_agesubcohorts_overall_migcertainty, output_name = 'IR_agesubcohorts_overall_migcertainty')

# IRR (univariable glm)
### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)),
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)),
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)),
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, migcertainty)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)),
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, migcertainty)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)),
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migcertainty)

## 65 years and over 
x <- glm.nb(conscount ~ as.factor(migcertainty) + offset(log(pyears)),
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_migcertainty_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migcertainty)


# Join age_subcohorts glm+IR and all_ages glm+IR
glm_migcertainty_combined_table <- dplyr::bind_rows(glm_migcertainty_0to15_table,glm_migcertainty_16to24_table,glm_migcertainty_25to34_table,glm_migcertainty_35to49_table,glm_migcertainty_50to64_table,glm_migcertainty_65plus_table)

univariable_migcertainty_agesubcohorts <- left_join(IR_agesubcohorts_overall_migcertainty, glm_migcertainty_combined_table,
                                                    by = c("migcertainty" = "names", "age_subcohort" = "age_subcohort")) 
univariable_migcertainty <- dplyr::bind_rows(univariable_migcertainty, univariable_migcertainty_agesubcohorts)  
univariable_migcertainty$age_subcohort <- as.factor(univariable_migcertainty$age_subcohort)
univariable_migcertainty$age_subcohort <- factor(univariable_migcertainty$age_subcohort,
                                                 levels = c('All_ages', '0-15 years','16-24 years', '25-34 years', '35-49 years','50-64 years', '>=65 years'))
univariable_migcertainty <- arrange(univariable_migcertainty, age_subcohort)


write_csv(univariable_migcertainty, "results/01_Consultations/full_cohort/univariable_migcertainty.csv")

## Multivariable model controlling for year + gender + as.factor(studyyear_agecat) + as.factor(imd) + prac_region

## All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts)
extract_glm_results_allages(x, 'multivariable_migcertainty_allages', cohort_England_2015_2020_annual_conscounts, migcertainty)
round(ci.lin(x,Exp=T),3)

## Age_subcohorts (not included in paper)

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migcertainty)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migcertainty)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migcertainty)

## 65 years and over
x <- glm.nb(conscount ~ as.factor(migcertainty) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts[cohort_England_2015_2020_annual_conscounts$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_multivariable_migcertainty_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migcertainty)

# Save .Rdata for forestplots
save(multivariable_migcertainty_allages, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/multivariable_migcertainty_allages.Rdata")
save(glm_multivariable_migcertainty_0to15, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/glm_multivariable_migcertainty_0to15.Rdata")
save(glm_multivariable_migcertainty_16to24, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/glm_multivariable_migcertainty_16to24.Rdata")
save(glm_multivariable_migcertainty_25to34, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/glm_multivariable_migcertainty_25to34.Rdata")
save(glm_multivariable_migcertainty_35to49, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/glm_multivariable_migcertainty_35to49.Rdata")
save(glm_multivariable_migcertainty_50to64, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/glm_multivariable_migcertainty_50to64.Rdata")
save(glm_multivariable_migcertainty_65plus, file="results/01_Consultations/full_cohort/multivariable_migcertainty_forestplot_inputs/glm_multivariable_migcertainty_65plus.Rdata")

# Join all ages and age_subcohorts glm results
multivariable_migcertainty <- dplyr::bind_rows(multivariable_migcertainty_allages_table, glm_multivariable_migcertainty_0to15_table,glm_multivariable_migcertainty_16to24_table,
                                               glm_multivariable_migcertainty_25to34_table,glm_multivariable_migcertainty_35to49_table, glm_multivariable_migcertainty_50to64_table,glm_multivariable_migcertainty_65plus_table)

write_csv(multivariable_migcertainty, "results/01_Consultations/full_cohort/multivariable_migcertainty.csv")


############################ APPENDIX (not included in paper) ###################### -----
## 12_Crude consultation rates yearly table ------

# Migrant_status
## annual
pyears_allpatients_overall_migstatus_studyyear <- aggregate(pyears ~ migrant_status + studyyear, cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migstatus_studyyear <- pyears_allpatients_overall_migstatus_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus_studyyear <- aggregate(conscount ~ migrant_status + studyyear, cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migstatus_studyyear <- inner_join(conscount_allpatients_overall_migstatus_studyyear,pyears_allpatients_overall_migstatus_studyyear, by= c("migrant_status" = "migrant_status", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migstatus_studyyear, output_name = 'IR_allpatients_overall_migstatus_studyyear')
## overall
pyears_allpatients_overall_migstatus <- aggregate(pyears ~ migrant_status, cohort_England_2015_2020_annual_conscounts, sum) 
pyears_allpatients_overall_migstatus <- pyears_allpatients_overall_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migstatus <- aggregate(conscount ~ migrant_status, cohort_England_2015_2020_annual_conscounts, sum) 
IR_allpatients_overall_migstatus <- inner_join(conscount_allpatients_overall_migstatus,pyears_allpatients_overall_migstatus, by= c("migrant_status" = "migrant_status"))
function_to_calculate_IR(IR_allpatients_overall_migstatus, output_name = 'IR_allpatients_overall_migstatus')


# Migcertainty
pyears_allpatients_overall_migcertainty_studyyear <- aggregate(pyears ~ migcertainty + studyyear, cohort_England_2015_2020_annual_conscounts, sum)
pyears_allpatients_overall_migcertainty_studyyear <- pyears_allpatients_overall_migcertainty_studyyear %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_allpatients_overall_migcertainty_studyyear <- aggregate(conscount ~ migcertainty + studyyear, cohort_England_2015_2020_annual_conscounts, sum)
IR_allpatients_overall_migcertainty_studyyear <- inner_join(conscount_allpatients_overall_migcertainty_studyyear,pyears_allpatients_overall_migcertainty_studyyear, by= c("migcertainty" = "migcertainty", 'studyyear' = 'studyyear'))
function_to_calculate_IR(IR_allpatients_overall_migcertainty_studyyear, output_name = 'IR_allpatients_overall_migcertainty_studyyear')
IR_allpatients_overall_migcertainty_studyyear <- IR_allpatients_overall_migcertainty_studyyear %>%
  filter(migcertainty %in% c('Definite', 'Probable')) %>%
  rename(migrant_status = migcertainty)

# Combine and convert to wide format for table
IRs_pre_pandemic_annual <- bind_rows(IR_allpatients_overall_migstatus_studyyear,
                                     IR_allpatients_overall_migcertainty_studyyear) %>%
  arrange(studyyear) %>%
  dplyr::select(c(studyyear, migrant_status,ir_ci)) %>%
  tidyr::pivot_wider(names_from = migrant_status, values_from = ir_ci) 
## 13_Study year----

## Multivariable IR and IRRs 

## 2015

studyyear_2015 <- cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2015)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2015)
extract_glm_results_allages(x, 'multivariable_2015_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2015[studyyear_2015$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2015[studyyear_2015$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2015[studyyear_2015$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2015[studyyear_2015$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2015[studyyear_2015$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2015[studyyear_2015$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2015_multivariable_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

# Save .Rdata for forestplots
save(multivariable_2015_allages, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/2015_multivariable_allages.Rdata")
save(glm_2015_multivariable_0to15, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2015_multivariable_0to15.Rdata")
save(glm_2015_multivariable_16to24, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2015_multivariable_16to24.Rdata")
save(glm_2015_multivariable_25to34, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2015_multivariable_25to34.Rdata")
save(glm_2015_multivariable_35to49, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2015_multivariable_35to49.Rdata")
save(glm_2015_multivariable_50to64, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2015_multivariable_50to64.Rdata")
save(glm_2015_multivariable_65plus, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2015_multivariable_65plus.Rdata")

# Join all ages and age_subcohorts glm results
multivariable_2015_migrant_status <- bind_rows(multivariable_2015_allages_table, glm_2015_multivariable_0to15_table,glm_2015_multivariable_16to24_table,
                                               glm_2015_multivariable_25to34_table,glm_2015_multivariable_35to49_table, glm_2015_multivariable_50to64_table,glm_2015_multivariable_65plus_table)


write_csv(multivariable_2015_migrant_status, "results/01_Consultations/full_cohort/multivariable_2015_migrant_status.csv")

## 2016

studyyear_2016 <- cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2016)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2016)
extract_glm_results_allages(x, 'multivariable_2016_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2016[studyyear_2016$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2016[studyyear_2016$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2016[studyyear_2016$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2016[studyyear_2016$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2016[studyyear_2016$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2016[studyyear_2016$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2016_multivariable_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

# Save .Rdata for forestplots
save(multivariable_2016_allages, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/2016_multivariable_allages.Rdata")
save(glm_2016_multivariable_0to15, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2016_multivariable_0to15.Rdata")
save(glm_2016_multivariable_16to24, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2016_multivariable_16to24.Rdata")
save(glm_2016_multivariable_25to34, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2016_multivariable_25to34.Rdata")
save(glm_2016_multivariable_35to49, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2016_multivariable_35to49.Rdata")
save(glm_2016_multivariable_50to64, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2016_multivariable_50to64.Rdata")
save(glm_2016_multivariable_65plus, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2016_multivariable_65plus.Rdata")

# Join all ages and age_subcohorts glm results
multivariable_2016_migrant_status <- bind_rows(multivariable_2016_allages_table, glm_2016_multivariable_0to15_table,glm_2016_multivariable_16to24_table,
                                               glm_2016_multivariable_25to34_table,glm_2016_multivariable_35to49_table,glm_2016_multivariable_50to64_table,glm_2016_multivariable_65plus_table)

write_csv(multivariable_2016_migrant_status, "results/01_Consultations/full_cohort/multivariable_2016_migrant_status.csv")

## 2017

studyyear_2017 <- cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2017)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2017)
extract_glm_results_allages(x, 'multivariable_2017_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2017[studyyear_2017$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2017[studyyear_2017$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2017[studyyear_2017$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2017[studyyear_2017$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2017[studyyear_2017$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2017[studyyear_2017$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2017_multivariable_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

# Save .Rdata for forestplots
save(multivariable_2017_allages, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/2017_multivariable_allages.Rdata")
save(glm_2017_multivariable_0to15, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2017_multivariable_0to15.Rdata")
save(glm_2017_multivariable_16to24, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2017_multivariable_16to24.Rdata")
save(glm_2017_multivariable_25to34, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2017_multivariable_25to34.Rdata")
save(glm_2017_multivariable_35to49, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2017_multivariable_35to49.Rdata")
save(glm_2017_multivariable_50to64, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2017_multivariable_50to64.Rdata")
save(glm_2017_multivariable_65plus, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2017_multivariable_65plus.Rdata")

# Join all ages and age_subcohorts glm results
multivariable_2017_migrant_status <- bind_rows(multivariable_2017_allages_table, glm_2017_multivariable_0to15_table,glm_2017_multivariable_16to24_table,
                                               glm_2017_multivariable_25to34_table,glm_2017_multivariable_35to49_table,glm_2017_multivariable_50to64_table,glm_2017_multivariable_65plus_table)

write_csv(multivariable_2017_migrant_status, "results/01_Consultations/full_cohort/multivariable_2017_migrant_status.csv")

## 2018

studyyear_2018 <- cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2018)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2018)
extract_glm_results_allages(x, 'multivariable_2018_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2018[studyyear_2018$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2018[studyyear_2018$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2018[studyyear_2018$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2018[studyyear_2018$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2018[studyyear_2018$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2018[studyyear_2018$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2018_multivariable_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

# Save .Rdata for forestplots
save(multivariable_2018_allages, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/2018_multivariable_allages.Rdata")
save(glm_2018_multivariable_0to15, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2018_multivariable_0to15.Rdata")
save(glm_2018_multivariable_16to24, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2018_multivariable_16to24.Rdata")
save(glm_2018_multivariable_25to34, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2018_multivariable_25to34.Rdata")
save(glm_2018_multivariable_35to49, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2018_multivariable_35to49.Rdata")
save(glm_2018_multivariable_50to64, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2018_multivariable_50to64.Rdata")
save(glm_2018_multivariable_65plus, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2018_multivariable_65plus.Rdata")

# Join all ages and age_subcohorts glm results
multivariable_2018_migrant_status <- bind_rows(multivariable_2018_allages_table, glm_2018_multivariable_0to15_table,glm_2018_multivariable_16to24_table,
                                               glm_2018_multivariable_25to34_table,glm_2018_multivariable_35to49_table,glm_2018_multivariable_50to64_table,glm_2018_multivariable_65plus_table)

write_csv(multivariable_2018_migrant_status, "results/01_Consultations/full_cohort/multivariable_2018_migrant_status.csv")

## 2019

studyyear_2019 <- cohort_England_2015_2020_annual_conscounts %>% filter(studyyear == 2019)

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2019)
extract_glm_results_allages(x, 'multivariable_2019_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)

## Age_subcohorts

### 0-15 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2019[studyyear_2019$age_subcohort=="0-15 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_0to15", "0-15 years", cohort_England_2015_2020_annual_conscounts, variable = migrant_status)

## 16-24 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2019[studyyear_2019$age_subcohort=="16-24 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_16to24", "16-24 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 25-34 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2019[studyyear_2019$age_subcohort=="25-34 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_25to34", "25-34 years", cohort_England_2015_2020_annual_conscounts, studyyear)

## 35-49 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2019[studyyear_2019$age_subcohort=="35-49 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_35to49", "35-49 years", cohort_England_2015_2020_annual_conscounts, studyyear)

# 50-64 years
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2019[studyyear_2019$age_subcohort=="50-64 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_50to64", "50-64 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

## 65 years and over
x <- glm.nb(conscount ~ as.factor(migrant_status) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = studyyear_2019[studyyear_2019$age_subcohort==">=65 years",])
extract_glm_results_agesubcohorts(x, "glm_2019_multivariable_65plus", ">=65 years", cohort_England_2015_2020_annual_conscounts, migrant_status)

# Save .Rdata for forestplots
save(multivariable_2019_allages, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/2019_multivariable_allages.Rdata")
save(glm_2019_multivariable_0to15, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2019_multivariable_0to15.Rdata")
save(glm_2019_multivariable_16to24, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2019_multivariable_16to24.Rdata")
save(glm_2019_multivariable_25to34, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2019_multivariable_25to34.Rdata")
save(glm_2019_multivariable_35to49, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2019_multivariable_35to49.Rdata")
save(glm_2019_multivariable_50to64, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2019_multivariable_50to64.Rdata")
save(glm_2019_multivariable_65plus, file="results/01_Consultations/full_cohort/multivariable_studyyear_forestplot_inputs/glm_2019_multivariable_65plus.Rdata")

# Join all ages and age_subcohorts glm results
multivariable_2019_migrant_status <- bind_rows(multivariable_2019_allages_table, glm_2019_multivariable_0to15_table,glm_2019_multivariable_16to24_table,
                                               glm_2019_multivariable_25to34_table,glm_2019_multivariable_35to49_table,glm_2019_multivariable_50to64_table,glm_2019_multivariable_65plus_table)


write_csv(multivariable_2019_migrant_status, "results/01_Consultations/full_cohort/multivariable_2019_migrant_status.csv")

## 14_Overdispersion and zero counts --------------------------------

# Visualise overdispersion
hist(cohort_England_2015_2020_annual_conscounts$conscount,
     breaks = 1000,
     xlim = c(0,150))

# Proportion of zero counts
zero_counts <- sum(cohort_England_2015_2020_annual_conscounts$conscount == 0)
total_counts <- nrow(cohort_England_2015_2020_annual_conscounts)
zero_counts/total_counts # 31.64%

## 15_Poisson vs NB vs ZIP vs ZINB model fit tests -------------------------------------------------

library(pscl) # Poisson, negative binomial, ZIP and ZINB
library(lmtest) # Likelihood ratio tests

## Multivariable model with all covariates

# Poisson
poisson <- glm(formula = conscount ~ as.factor(migrant_status) + offset(log(pyears)) +
                 as.factor(studyyear) + as.factor(studyyear_agecat) +
                 as.factor(imd) + as.factor(prac_region) + as.factor(gender),
               data = cohort_England_2015_2020_annual_conscounts, family = "poisson")
names <- names(coef(poisson))
estimate<- exp(coef(poisson))
confint.default <- exp(confint.default(poisson))
p <- coef(summary(poisson))[,4]
poisson_outputs <- cbind(names,estimate,confint.default,p)

# Negative binomial
nb <- glm.nb(formula = conscount ~ as.factor(migrant_status) + offset(log(pyears)) +
               as.factor(studyyear) + as.factor(studyyear_agecat) +
               as.factor(imd) + as.factor(prac_region) + as.factor(gender),
             data = cohort_England_2015_2020_annual_conscounts)
names <- names(coef(nb))
estimate<- exp(coef(nb))
confint.default <- exp(confint.default(nb))
p <- coef(summary(nb))[,4]
nb_outputs <- cbind(names,estimate,confint.default,p)

# ZIP
zip <- zeroinfl(conscount ~ as.factor(migrant_status) + offset(log(pyears)) +
                  as.factor(studyyear) + as.factor(studyyear_agecat) +
                  as.factor(imd) + as.factor(prac_region) + as.factor(gender),
                data=cohort_England_2015_2020_annual_conscounts)
names <- names(coef(zip))
estimate<- exp(coef(zip))
confint.default <- exp(confint.default(zip))
p <- coef(summary(zip))[,4] 
zip_outputs <- cbind(names,estimate,confint.default,p)

# ZINB
zinb <- zeroinfl(formula = conscount ~ as.factor(migrant_status) + offset(log(pyears)) +
                   as.factor(studyyear) + as.factor(studyyear_agecat) +
                   as.factor(imd) + as.factor(prac_region) + as.factor(gender),
                 dist = "negbin",
                 data = cohort_England_2015_2020_annual_conscounts)
names <- names(coef(zinb))
estimate<- exp(coef(zinb))
confint.default <- exp(confint.default(zinb))
summary <- summary(zinb)$coef
p <- coef(summary(zinb))[,4] 
zinb_outputs <- cbind(names,estimate,confint.default,p)

# Likelihood ratio test NB vs. Poisson
lrtest(poisson, nb) #p<0.001 NB better

# Vuong test Poisson vs. ZIP
vuong(poisson, zip) # ZIP better p<0.001

# Likelihood ratio test ZINB vs. ZIP
lrtest(zip,zinb) #p<0.001, ZINB better

# Vuong test ZINB vs NB
vuong(nb,zinb) # ZINB better (minimally) - but not in line with our research questions

## 16_Remove files ------------------------------------------------------------

rm(list = ls()[!(ls() %in% c('cohort_England_2015_2020_annual_conscounts',
                             'function_to_calculate_IR', 'extract_glm_results_allages',
                             'extract_glm_results_agesubcohorts'))])
gc()



