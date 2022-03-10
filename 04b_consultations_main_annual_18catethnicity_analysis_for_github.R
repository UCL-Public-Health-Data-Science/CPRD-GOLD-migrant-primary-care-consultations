## 0_Description ---------------------------------------------------------------------------

# Migrants' primary care utilisation before and during the COVID-19 pandemic in England: A cohort study
# Annual pre-pandemic analysis: Migration status and 18-category ethnicity 
# Date started: 30/11/2020
# Author(s): Yamina Boukari / Claire Zhang / Neha Pathak 
# QC (Date): Claire Zhang (05/03/2022)

## 1_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc, epitools, data.table, epiR, MASS, psych, forestplot, dsr, gmodels, Epi, msm)

## 2_Set working directory ------------------------------------------------------------------

setwd("filepath")

## 3_Load datasets and functions ------------------------------------------------------

# Datasets 

load(file = "filepath/cohort_England_2015_2020_annual_conscounts.Rdata")

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

## 4_Select needed variables and prepare data for modelling  -------------------------------------------------

cohort_England_2015_2020_annual_conscounts <- dplyr::select(cohort_England_2015_2020_annual_conscounts, 
                                                            c(patid, pyears, conscount, migrant_status,
                                                              age_subcohort, imd, ethnicat, gender,
                                                              studyyear, studyyear_agecat, prac_region,
                                                              migcertainty, cohort_entry))

# Relevel
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, migrant_status <- relevel (migrant_status, ref="Non-migrant")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, gender <- relevel (gender, ref="Male")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, imd <- relevel (imd, ref="IMD 1"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, prac_region <- relevel (prac_region, ref="London")) 
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, age_subcohort <- relevel (age_subcohort, ref="0-15 years"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, studyyear_agecat <- relevel (studyyear_agecat, ref="0-5 years"))
cohort_England_2015_2020_annual_conscounts <- within(cohort_England_2015_2020_annual_conscounts, ethnicat <- relevel (ethnicat, ref="British"))


## 5_IMPORTANT STEP: Filter out 2020 for modelling ------
cohort_England_2015_2020_annual_conscounts <- cohort_England_2015_2020_annual_conscounts %>%
  filter(studyyear != 2020)

## 6_Ethnicity (ethnicat) - interaction term, additive and multiplicative effects ----

# All patients (i.e. all ages)

# Generate IRR via negative binomial regression
x <- glm.nb(conscount ~ as.factor(migrant_status)*as.factor(ethnicat) + as.factor(studyyear) +
              as.factor(gender) + as.factor(studyyear_agecat) + as.factor(imd) + as.factor(prac_region) + offset(log(pyears)) ,
            data = cohort_England_2015_2020_annual_conscounts)

round(ci.lin(x,Exp=T),3)

extract_glm_results_allages(x, 'multivariable_ethnicat18interaction_allages', cohort_England_2015_2020_annual_conscounts, migrant_status)

# Save .Rdata for forestplots
save(multivariable_ethnicat18interaction_allages, file="filepath/multivariable_ethnicat18interaction_allages.Rdata")

write_csv(multivariable_ethnicat18interaction_allages, "filepath/multivariable_ethnicat18interaction_allages.csv")

# RR (95% CI) for non-migrants in each ethnicity compared to White British non-migrants
# For table

ethnic_groups <- multivariable_ethnicat18interaction_allages$names[3:20]

british_non_migrant <- data.frame(Ethnicity = 'British', RR_nonmigrantsVsWBNM = '1.0')
non_migrants_by_ethnicity <- multivariable_ethnicat18interaction_allages[3:20,] %>%
  mutate(ci = paste(lower, upper, sep ="-"),
         RR_nonmigrantsVsWBNM = paste0(estimate, ' (', ci, ")")) %>%
  mutate(Ethnicity = ethnic_groups) %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, estimate, upper, lower)
non_migrants_by_ethnicity <- bind_rows(british_non_migrant, non_migrants_by_ethnicity)

## 7_Calculate interaction effects - migrants in each ethnicity compared to British non-migrant (BNM) reference group ----

get_RR_and_95_CI <- function(dataframe) {
  output <- dataframe %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$Lower.CI, output$Upper.CI, sep ="-")
  output$RR <- paste0(output$Estimate, ' (', output$ci, ")")
  output_name <- deparse(substitute(dataframe))
  assign(x = output_name, value = output, envir = globalenv())
}

## British migrants ----
british_migrant_vsBNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(british_migrant_vsBNM)
save(british_migrant_vsBNM, file="filepath/british_migrant_vsBNM.Rdata")
british_migrant <- british_migrant_vsBNM %>%
  dplyr::select(RR)

## Irish migrants ----
irish_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Irish' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Irish' = 1), conf=.95))
get_RR_and_95_CI(irish_migrant_vsBNM)
save(irish_migrant_vsBNM, file="filepath/irish_migrant_vsBNM.Rdata")
irish_migrant <- irish_migrant_vsBNM %>%
  dplyr::select(RR)

## Gypsy or Irish traveller (GIT) migrants ----
GIT_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Gypsy or Irish Traveller' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Gypsy or Irish Traveller' = 1), conf=.95))
get_RR_and_95_CI(GIT_migrant_vsBNM)
save(GIT_migrant_vsBNM, file="filepath/GIT_migrant_vsBNM.Rdata")
GIT_migrant <- GIT_migrant_vsBNM %>%
  dplyr::select(RR)

## Other White ----
otherWhite_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Other White' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other White' = 1), conf=.95))
get_RR_and_95_CI(otherWhite_migrant_vsBNM)
save(otherWhite_migrant_vsBNM, file="filepath/otherWhite_migrant_vsBNM.Rdata")
otherWhite_migrant <- otherWhite_migrant_vsBNM %>%
  dplyr::select(RR)

## Mixed White and Black Caribbean ----
mixedWhiteBlackCaribbean_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Mixed White and Black Caribbean' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Mixed White and Black Caribbean' = 1), conf=.95))
get_RR_and_95_CI(mixedWhiteBlackCaribbean_migrant_vsBNM)
save(mixedWhiteBlackCaribbean_migrant_vsBNM, file="filepath/mixedWhiteBlackCaribbean_migrant_vsBNM.Rdata")
mixedWhiteBlackCaribbean_migrant <- mixedWhiteBlackCaribbean_migrant_vsBNM %>%
  dplyr::select(RR)

## Mixed White and Black African ----
mixedWhiteBlackAfrican_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Mixed White and Black African' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Mixed White and Black African' = 1), conf=.95))
get_RR_and_95_CI(mixedWhiteBlackAfrican_migrant_vsBNM)
save(mixedWhiteBlackAfrican_migrant_vsBNM, file="filepath/mixedWhiteBlackAfrican_migrant_vsBNM.Rdata")
mixedWhiteBlackAfrican_migrant <- mixedWhiteBlackAfrican_migrant_vsBNM %>%
  dplyr::select(RR)

## Mixed White and Asian ----
mixedWhiteAsian_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Mixed White and Asian' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Mixed White and Asian' = 1), conf=.95))
get_RR_and_95_CI(mixedWhiteAsian_migrant_vsBNM)
save(mixedWhiteAsian_migrant_vsBNM, file="filepath/mixedWhiteAsian_migrant_vsBNM.Rdata")
mixedWhiteAsian_migrant <- mixedWhiteAsian_migrant_vsBNM %>%
  dplyr::select(RR)

## Other Mixed ----
otherMixed_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Other Mixed' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other Mixed' = 1), conf=.95))
get_RR_and_95_CI(otherMixed_migrant_vsBNM)
save(otherMixed_migrant_vsBNM, file="filepath/otherMixed_migrant_vsBNM.Rdata")
otherMixed_migrant <- otherMixed_migrant_vsBNM %>%
  dplyr::select(RR)



## Indian----
indian_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Indian' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Indian' = 1), conf=.95))
get_RR_and_95_CI(indian_migrant_vsBNM)
save(indian_migrant_vsBNM, file="filepath/indian_migrant_vsBNM.Rdata")
indian_migrant <- indian_migrant_vsBNM %>%
  dplyr::select(RR)

## Pakistani----
pakistani_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Pakistani' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Pakistani' = 1), conf=.95))
get_RR_and_95_CI(pakistani_migrant_vsBNM)
save(pakistani_migrant_vsBNM, file="filepath/pakistani_migrant_vsBNM.Rdata")
pakistani_migrant <- pakistani_migrant_vsBNM %>%
  dplyr::select(RR)

## Bangladeshi----
bangladeshi_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Bangladeshi' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Bangladeshi' = 1), conf=.95))
get_RR_and_95_CI(bangladeshi_migrant_vsBNM)
save(bangladeshi_migrant_vsBNM, file="filepath/bangladeshi_migrant_vsBNM.Rdata")
bangladeshi_migrant <- bangladeshi_migrant_vsBNM %>%
  dplyr::select(RR)

## Chinese----
chinese_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Chinese' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Chinese' = 1), conf=.95))
get_RR_and_95_CI(chinese_migrant_vsBNM)
save(chinese_migrant_vsBNM, file="filepath/chinese_migrant_vsBNM.Rdata")
chinese_migrant <- chinese_migrant_vsBNM %>%
  dplyr::select(RR)

## Other Asian----
otherAsian_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Other Asian' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other Asian' = 1), conf=.95))
get_RR_and_95_CI(otherAsian_migrant_vsBNM)
save(otherAsian_migrant_vsBNM, file="filepath/otherAsian_migrant_vsBNM.Rdata")
otherAsian_migrant <- otherAsian_migrant_vsBNM %>%
  dplyr::select(RR)

## African----
african_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)African' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)African' = 1), conf=.95))
get_RR_and_95_CI(african_migrant_vsBNM)
save(african_migrant_vsBNM, file="filepath/african_migrant_vsBNM.Rdata")
african_migrant <- african_migrant_vsBNM %>%
  dplyr::select(RR)

## Caribbean----
caribbean_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Caribbean' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Caribbean' = 1), conf=.95))
get_RR_and_95_CI(caribbean_migrant_vsBNM)
save(caribbean_migrant_vsBNM, file="filepath/caribbean_migrant_vsBNM.Rdata")
caribbean_migrant <- caribbean_migrant_vsBNM %>%
  dplyr::select(RR)

## Other Black----
otherBlack_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Other Black' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other Black' = 1), conf=.95))
get_RR_and_95_CI(otherBlack_migrant_vsBNM)
save(otherBlack_migrant_vsBNM, file="filepath/otherBlack_migrant_vsBNM.Rdata")
otherBlack_migrant <- otherBlack_migrant_vsBNM %>%
  dplyr::select(RR)

## Arab ----
arab_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Arab' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Arab' = 1), conf=.95))
get_RR_and_95_CI(arab_migrant_vsBNM)
save(arab_migrant_vsBNM, file="filepath/arab_migrant_vsBNM.Rdata")
arab_migrant <- arab_migrant_vsBNM %>%
  dplyr::select(RR)

## Any other ethnic group ----
other_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Any other ethnic group' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Any other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_migrant_vsBNM)
save(other_migrant_vsBNM, file="filepath/other_migrant_vsBNM.Rdata")
other_migrant <- other_migrant_vsBNM %>%
  dplyr::select(RR)

## Unknown ----
unknown_migrant_vsBNM <- exp(estimable(x, c('as.factor(ethnicat)Unknown' = 1, 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_migrant_vsBNM)
save(unknown_migrant_vsBNM, file="filepath/unknown_migrant_vsBNM.Rdata")
unknown_migrant <- unknown_migrant_vsBNM %>%
  dplyr::select(RR)


## 8_Migrants within each ethnicity strata ----

## British migrants ----
british_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1), conf=.95))
get_RR_and_95_CI(british_strata_vsNM)
save(british_strata_vsNM, file="filepath/british_strata_vsNM.Rdata")
british_strata <- british_strata_vsNM %>%
  dplyr::select(RR)

## Irish migrants ----
irish_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Irish' = 1), conf=.95))
get_RR_and_95_CI(irish_strata_vsNM)
save(irish_strata_vsNM, file="filepath/irish_strata_vsNM.Rdata")
irish_strata <- irish_strata_vsNM %>%
  dplyr::select(RR)

## Gypsy or Irish traveller (GIT) migrants ----
GIT_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Gypsy or Irish Traveller' = 1), conf=.95))
get_RR_and_95_CI(GIT_strata_vsNM)
save(GIT_strata_vsNM, file="filepath/GIT_strata_vsNM.Rdata")
GIT_strata <- GIT_strata_vsNM %>%
  dplyr::select(RR)

## Other White ----
otherWhite_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other White' = 1), conf=.95))
get_RR_and_95_CI(otherWhite_strata_vsNM)
save(otherWhite_strata_vsNM, file="filepath/otherWhite_strata_vsNM.Rdata")
otherWhite_strata <- otherWhite_strata_vsNM %>%
  dplyr::select(RR)

## Mixed White and Black Caribbean ----
mixedWhiteBlackCaribbean_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Mixed White and Black Caribbean' = 1), conf=.95))
get_RR_and_95_CI(mixedWhiteBlackCaribbean_strata_vsNM)
save(mixedWhiteBlackCaribbean_strata_vsNM, file="filepath/mixedWhiteBlackCaribbean_strata_vsNM.Rdata")
mixedWhiteBlackCaribbean_strata <- mixedWhiteBlackCaribbean_strata_vsNM %>%
  dplyr::select(RR)

## Mixed White and Black African ----
mixedWhiteBlackAfrican_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Mixed White and Black African' = 1), conf=.95))
get_RR_and_95_CI(mixedWhiteBlackAfrican_strata_vsNM)
save(mixedWhiteBlackAfrican_strata_vsNM, file="filepath/mixedWhiteBlackAfrican_strata_vsNM.Rdata")
mixedWhiteBlackAfrican_strata <- mixedWhiteBlackAfrican_strata_vsNM %>%
  dplyr::select(RR)

## Mixed White and Asian ----
mixedWhiteAsian_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Mixed White and Asian' = 1), conf=.95))
get_RR_and_95_CI(mixedWhiteAsian_strata_vsNM)
save(mixedWhiteAsian_strata_vsNM, file="filepath/mixedWhiteAsian_strata_vsNM.Rdata")
mixedWhiteAsian_strata <- mixedWhiteAsian_strata_vsNM %>%
  dplyr::select(RR)

## Other Mixed ----
otherMixed_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other Mixed' = 1), conf=.95))
get_RR_and_95_CI(otherMixed_strata_vsNM)
save(otherMixed_strata_vsNM, file="filepath/otherMixed_strata_vsNM.Rdata")
otherMixed_strata <- otherMixed_strata_vsNM %>%
  dplyr::select(RR)



## Indian----
indian_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Indian' = 1), conf=.95))
get_RR_and_95_CI(indian_strata_vsNM)
save(indian_strata_vsNM, file="filepath/indian_strata_vsNM.Rdata")
indian_strata <- indian_strata_vsNM %>%
  dplyr::select(RR)

## Pakistani----
pakistani_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Pakistani' = 1), conf=.95))
get_RR_and_95_CI(pakistani_strata_vsNM)
save(pakistani_strata_vsNM, file="filepath/pakistani_strata_vsNM.Rdata")
pakistani_strata <- pakistani_strata_vsNM %>%
  dplyr::select(RR)

## Bangladeshi----
bangladeshi_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Bangladeshi' = 1), conf=.95))
get_RR_and_95_CI(bangladeshi_strata_vsNM)
save(bangladeshi_strata_vsNM, file="filepath/bangladeshi_strata_vsNM.Rdata")
bangladeshi_strata <- bangladeshi_strata_vsNM %>%
  dplyr::select(RR)

## Chinese----
chinese_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Chinese' = 1), conf=.95))
get_RR_and_95_CI(chinese_strata_vsNM)
save(chinese_strata_vsNM, file="filepath/chinese_strata_vsNM.Rdata")
chinese_strata <- chinese_strata_vsNM %>%
  dplyr::select(RR)

## Other Asian----
otherAsian_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other Asian' = 1), conf=.95))
get_RR_and_95_CI(otherAsian_strata_vsNM)
save(otherAsian_strata_vsNM, file="results/01_Consultations/full_cohort/multivariable_ethnicity_forestplot_inputs/otherAsian_strata_vsNM.Rdata")
otherAsian_strata <- otherAsian_strata_vsNM %>%
  dplyr::select(RR)

## African----
african_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)African' = 1), conf=.95))
get_RR_and_95_CI(african_strata_vsNM)
save(african_strata_vsNM, file="filepath/african_strata_vsNM.Rdata")
african_strata <- african_strata_vsNM %>%
  dplyr::select(RR)

## Caribbean----
caribbean_strata_vsNM <- exp(estimable(x, c( 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Caribbean' = 1), conf=.95))
get_RR_and_95_CI(caribbean_strata_vsNM)
save(caribbean_strata_vsNM, file="filepath/caribbean_strata_vsNM.Rdata")
caribbean_strata <- caribbean_strata_vsNM %>%
  dplyr::select(RR)

## Other Black----
otherBlack_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Other Black' = 1), conf=.95))
get_RR_and_95_CI(otherBlack_strata_vsNM)
save(otherBlack_strata_vsNM, file="filepath/otherBlack_strata_vsNM.Rdata")
otherBlack_strata <- otherBlack_strata_vsNM %>%
  dplyr::select(RR)

## Arab ----
arab_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Arab' = 1), conf=.95))
get_RR_and_95_CI(arab_strata_vsNM)
save(arab_strata_vsNM, file="filepath/arab_strata_vsNM.Rdata")
arab_strata <- arab_strata_vsNM %>%
  dplyr::select(RR)

## Any other ethnic group ----
other_strata_vsNM <- exp(estimable(x, c( 'as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Any other ethnic group' = 1), conf=.95))
get_RR_and_95_CI(other_strata_vsNM)
save(other_strata_vsNM, file="filepath/other_strata_vsNM.Rdata")
other_strata <- other_strata_vsNM %>%
  dplyr::select(RR)

## Unknown ----
unknown_strata_vsNM <- exp(estimable(x, c('as.factor(migrant_status)Migrant' = 1, 'as.factor(migrant_status)Migrant:as.factor(ethnicat)Unknown' = 1), conf=.95))
get_RR_and_95_CI(unknown_strata_vsNM)
save(unknown_strata_vsNM, file="filepath/unknown_strata_vsNM.Rdata")
unknown_strata <- unknown_strata_vsNM %>%
  dplyr::select(RR)

## 9_Combine into one dataframe ----

Ethnicity <- data.frame(levels(cohort_England_2015_2020_annual_conscounts$ethnicat))
migrants_by_ethnicity <- bind_rows(british_migrant, irish_migrant, GIT_migrant,
                                   otherWhite_migrant, mixedWhiteBlackCaribbean_migrant,
                                   mixedWhiteBlackAfrican_migrant, mixedWhiteAsian_migrant,
                                   otherMixed_migrant, indian_migrant, pakistani_migrant,
                                   bangladeshi_migrant, chinese_migrant, otherAsian_migrant,
                                   african_migrant, caribbean_migrant, otherBlack_migrant,
                                   arab_migrant,other_migrant, unknown_migrant)
Ethniciy_strata <- bind_rows(british_strata, irish_strata, GIT_strata, otherWhite_strata,
                             mixedWhiteBlackCaribbean_strata, mixedWhiteBlackAfrican_strata,
                             mixedWhiteAsian_strata, otherMixed_strata, indian_strata, pakistani_strata,
                             bangladeshi_strata, chinese_strata, otherAsian_strata, african_strata,
                             caribbean_strata, otherBlack_strata, arab_strata, other_strata, unknown_strata)
ethnicity_table <- bind_cols(Ethnicity, migrants_by_ethnicity, Ethniciy_strata)
row.names(ethnicity_table) <- NULL

ethnicity_table <- ethnicity_table %>%
  rename(
    Ethnicity = levels.cohort_England_2015_2020_annual_conscounts.ethnicat.,
    RR_migrantsVsWBNM = RR...2,
    RR_MvsNM_each_ethnicity = RR...3
  ) %>%
  right_join(non_migrants_by_ethnicity, by = 'Ethnicity') %>%
  dplyr::select(Ethnicity, RR_nonmigrantsVsWBNM, RR_migrantsVsWBNM, RR_MvsNM_each_ethnicity)

# Save dataframe for knitting in RMD file
save(ethnicity_table, file="filepath/ethnicity18cat_interactions_table.Rdata")


## 10_RERI and multiplicative effect -----

# get coefficients and RRs

get_RERI_and_multiplicative_effect <- function(ethnicity){
  exposure1 <- names(coef(x))[2]
  exposure2 <- paste0('as.factor(ethnicat)',ethnicity)
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
  output <- data.frame(ethnicity = ethnicity, 
                       RERI = RERI,
                       RERI_upper_CI = upper_CI,
                       RERI_lower_CI = lower_CI,
                       multiplicative_effect = RR_interaction,
                       mult_lower_CI = beta_interaction_CI[1],
                       mult_upper_CI = beta_interaction_CI[2])
  row.names(output) <- NULL
  assign(x = paste0('interaction_effects_', ethnicity), value = output, envir = globalenv())}

get_RERI_and_multiplicative_effect('Irish')
get_RERI_and_multiplicative_effect('Gypsy or Irish Traveller')
get_RERI_and_multiplicative_effect('Other White')
get_RERI_and_multiplicative_effect('Mixed White and Black Caribbean')
get_RERI_and_multiplicative_effect('Mixed White and Black African')
get_RERI_and_multiplicative_effect('Mixed White and Asian')
get_RERI_and_multiplicative_effect('Other Mixed')
get_RERI_and_multiplicative_effect('Indian')
get_RERI_and_multiplicative_effect('Pakistani')
get_RERI_and_multiplicative_effect('Bangladeshi')
get_RERI_and_multiplicative_effect('Chinese')
get_RERI_and_multiplicative_effect('Other Asian')
get_RERI_and_multiplicative_effect('African')
get_RERI_and_multiplicative_effect('Caribbean')
get_RERI_and_multiplicative_effect('Other Black')
get_RERI_and_multiplicative_effect('Arab')
get_RERI_and_multiplicative_effect('Any other ethnic group')
get_RERI_and_multiplicative_effect('Unknown')

int_effects_all_ethnicat18 <- bind_rows(interaction_effects_Irish, `interaction_effects_Gypsy or Irish Traveller`,
                                         `interaction_effects_Other White`, `interaction_effects_Mixed White and Black Caribbean`,
                                         `interaction_effects_Mixed White and Black African`, `interaction_effects_Mixed White and Asian`,
                                        `interaction_effects_Other Mixed`, interaction_effects_Indian, interaction_effects_Pakistani,
                                        interaction_effects_Bangladeshi, interaction_effects_Chinese, `interaction_effects_Other Asian`,
                                        interaction_effects_African, interaction_effects_Caribbean, `interaction_effects_Other Black`,
                                        interaction_effects_Arab, `interaction_effects_Any other ethnic group`, interaction_effects_Unknown)

save(int_effects_all_ethnicat18, file="filepath/int_effects_all_ethnicat18.Rdata")

