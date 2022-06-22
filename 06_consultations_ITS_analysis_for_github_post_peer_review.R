## 0_Description ---------------------------------------------------------------------------

# Migrants' primary care utilisation before and during the COVID-19 pandemic in England: A cohort study
# 2015-2020 interrupted time series analysis - code for analyses included in the paper
# Date started: 25/03/2021
# Author(s): Yamina Boukari / Claire Zhang 
# QA (date): Claire Zhang (26/02/2022)

## 1_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc, epitools, data.table, epiR, MASS, psych, forestplot, foreign,tsModel,lmtest,
               Epi, multcomp,splines, vcd, here, stringr, patchwork, gmodels, gtsummary, scales, lmtest, cowplot, grid, gridExtra, sandwich, car, msm)


## 2_Set working directory ------------------------------------------------------------------

setwd("filepath")

## 3_Load datasets and functions ------------------------------------------------------

# Dataset 
## Main analysis
load(file = "filepath")
## Secondary analyses
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
## Sensitivity analyses
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")
load(file = "filepath")

# Functions

function_to_calculate_IR <- function(dataset, output_name) {
  output <- pois.exact(dataset$conscount, dataset$pyears, conf.level=0.95)
  output <- output %>% 
    rename(incidence_rate = rate)
  output <- bind_cols(dataset, output) 
  output <- output %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$ir_ci <- paste0(output$incidence_rate, ' (', output$ci, ')')
  output <-  dplyr::select(output, -c(conf.level, x, pt, ci))
  output <- output %>% 
    rename(person_years = pyears)
  assign(x = output_name, value = output, envir = globalenv())
} # per 1 person-year


function_to_calculate_IR_F2F <- function(dataset, output_name) {
  output <- pois.exact(dataset$facetoface, dataset$pyears, conf.level=0.95)
  output <- output %>% 
    rename(incidence_rate = rate)
  output <- bind_cols(dataset, output) 
  output <- output %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$ir_ci <- paste0(output$incidence_rate, ' (', output$ci, ')')
  output <-  dplyr::select(output, -c(conf.level, x, pt, ci))
  output <- output %>% 
    rename(person_years = pyears)
  assign(x = output_name, value = output, envir = globalenv())
}

function_to_calculate_IR_phone <- function(dataset, output_name) {
  output <- pois.exact(dataset$phone, dataset$pyears, conf.level=0.95)
  output <- output %>% 
    rename(incidence_rate = rate)
  output <- bind_cols(dataset, output) 
  output <- output %>% mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$ir_ci <- paste0(output$incidence_rate, ' (',output$ci, ')')
  output <-  dplyr::select(output, -c(conf.level, x, pt, ci))
  output <- output %>% 
    rename(person_years = pyears)
  assign(x = output_name, value = output, envir = globalenv())
}

extract_glm_results_allages <- function (x, output_name, dataset, variable) {
  names <- names(coef(x))
  names <- sub(".*)", "", names) # extracts all text after ')' from names to avoid having to recode later (excluding the intercept, the below lines of code are deal with this)
  if (is.factor(eval(substitute(variable), dataset)) == F) {
    variable <- as.factor(eval(substitute(variable), dataset))
  } # if a variable isn't a factor (which is needed for the below line), e.g. studyyear, it converts it to a factor
  #names[1] <- levels(eval(substitute(variable), dataset))[1] # takes the reference level of a specified factor variable and replaces it with 'Intercept' in names
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
  #output$estimate[1] <- 1.00
  #output$lower[1] <- 1.00
  #output$upper[1] <- 1.00
  #output$ci[1] <- "1.00-1.00"
  output$irr_ci <- paste(output$estimate, output$ci, sep =",")
  output_table <- dplyr::select(output,names, irr_ci )
  output_table <- output_table %>%
    mutate(age_subcohort = 'All_ages')
  assign(x = output_name, value = output, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
  assign(x = paste0(output_name,'_table'), value = output_table, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
}

# 4_MAIN ANALYSIS: England ####################################

# 4a_Summary statistics #########

# Number of events 
number_cons_migstatus <- aggregate(conscount ~migrant_status, England_2015_2020_aggregate_weekly_conscounts, sum)
number_cons_total <- data.frame(migrant_status = 'All', conscount = as.numeric(sum(number_cons_migstatus$conscount)))
number_cons_migstatus <- bind_rows(number_cons_migstatus, number_cons_total)
write.csv(number_cons_migstatus, "filepath" )


# Aggregate IRs pre and during pandemic by migrant_status and mig_certainty

## All consultations
# migrant_status
pyears_pandemic_migstatus <- aggregate(pyears ~ migrant_status + lockdown1, England_2015_2020_aggregate_weekly_conscounts, sum) 
pyears_pandemic_migstatus <- pyears_pandemic_migstatus %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_pandemic_migstatus <- aggregate(conscount ~ migrant_status +lockdown1, England_2015_2020_aggregate_weekly_conscounts, sum) 
IR_pandemic_migstatus <- inner_join(conscount_pandemic_migstatus,pyears_pandemic_migstatus, by= c("migrant_status" = 'migrant_status', 'lockdown1' = 'lockdown1'))
function_to_calculate_IR(IR_pandemic_migstatus, output_name = 'IR_pandemic_migstatus')

# migcertainty
pyears_pandemic_migcertainty <- aggregate(pyears ~ migcertainty + lockdown1, England_2015_2020_aggregate_weekly_conscounts_migcertainty, sum) 
pyears_pandemic_migcertainty <- pyears_pandemic_migcertainty %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_pandemic_migcertainty <- aggregate(conscount ~ migcertainty +lockdown1, England_2015_2020_aggregate_weekly_conscounts_migcertainty, sum) 
IR_pandemic_migcertainty <- inner_join(conscount_pandemic_migcertainty,pyears_pandemic_migcertainty, by= c("migcertainty" = 'migcertainty', 'lockdown1' = 'lockdown1'))
function_to_calculate_IR(IR_pandemic_migcertainty, output_name = 'IR_pandemic_migcertainty')
IR_pandemic_migcertainty <- IR_pandemic_migcertainty %>%
  filter(migcertainty %in% c('Definite', 'Probable')) %>%
  rename(migrant_status = migcertainty)

# Combine and convert to wide format for table
IRs_pandemic <- bind_rows(IR_pandemic_migstatus,
                          IR_pandemic_migcertainty) %>%
  arrange(lockdown1) %>%
  dplyr::select(c(migrant_status, lockdown1, ir_ci)) %>%
  tidyr::pivot_wider(names_from = migrant_status, values_from = ir_ci) %>%
  mutate(type = 'All consultations') %>%
  relocate(type)

## F2F consultations
# migrant_status
pyears_pandemic_migstatus_F2F <- aggregate(pyears ~ migrant_status + lockdown1, England_2015_2020_aggregate_weekly_conscounts, sum) 
pyears_pandemic_migstatus_F2F <- pyears_pandemic_migstatus_F2F %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_pandemic_migstatus_F2F <- aggregate(facetoface ~ migrant_status +lockdown1, England_2015_2020_aggregate_weekly_conscounts, sum) 
IR_pandemic_migstatus_F2F <- inner_join(conscount_pandemic_migstatus_F2F,pyears_pandemic_migstatus_F2F, by= c("migrant_status" = 'migrant_status', 'lockdown1' = 'lockdown1'))
function_to_calculate_IR_F2F(IR_pandemic_migstatus_F2F, output_name = 'IR_pandemic_migstatus_F2F')

# migcertainty
pyears_pandemic_migcertainty_F2F <- aggregate(pyears ~ migcertainty + lockdown1, England_2015_2020_aggregate_weekly_conscounts_migcertainty, sum) 
pyears_pandemic_migcertainty_F2F <- pyears_pandemic_migcertainty_F2F %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_pandemic_migcertainty_F2F <- aggregate(facetoface ~ migcertainty +lockdown1, England_2015_2020_aggregate_weekly_conscounts_migcertainty, sum) 
IR_pandemic_migcertainty_F2F <- inner_join(conscount_pandemic_migcertainty_F2F,pyears_pandemic_migcertainty_F2F, by= c("migcertainty" = 'migcertainty', 'lockdown1' = 'lockdown1'))
function_to_calculate_IR_F2F(IR_pandemic_migcertainty_F2F, output_name = 'IR_pandemic_migcertainty_F2F')
IR_pandemic_migcertainty_F2F <- IR_pandemic_migcertainty_F2F %>%
  filter(migcertainty %in% c('Definite', 'Probable')) %>%
  rename(migrant_status = migcertainty)

# Combine and convert to wide format for table
IRs_pandemic_F2F <- bind_rows(IR_pandemic_migstatus_F2F,
                          IR_pandemic_migcertainty_F2F) %>%
  arrange(lockdown1) %>%
  dplyr::select(c(migrant_status, lockdown1, ir_ci)) %>%
  tidyr::pivot_wider(names_from = migrant_status, values_from = ir_ci) %>%
  mutate(type = 'Face-to-face consultations') %>%
  relocate(type)

## Phone consultations
# migrant_status
pyears_pandemic_migstatus_phone <- aggregate(pyears ~ migrant_status + lockdown1, England_2015_2020_aggregate_weekly_conscounts, sum) 
pyears_pandemic_migstatus_phone <- pyears_pandemic_migstatus_phone %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_pandemic_migstatus_phone <- aggregate(phone ~ migrant_status +lockdown1, England_2015_2020_aggregate_weekly_conscounts, sum) 
IR_pandemic_migstatus_phone <- inner_join(conscount_pandemic_migstatus_phone,pyears_pandemic_migstatus_phone, by= c("migrant_status" = 'migrant_status', 'lockdown1' = 'lockdown1'))
function_to_calculate_IR_phone(IR_pandemic_migstatus_phone, output_name = 'IR_pandemic_migstatus_phone')

# migcertainty
pyears_pandemic_migcertainty_phone <- aggregate(pyears ~ migcertainty + lockdown1, England_2015_2020_aggregate_weekly_conscounts_migcertainty, sum) 
pyears_pandemic_migcertainty_phone <- pyears_pandemic_migcertainty_phone %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_pandemic_migcertainty_phone <- aggregate(phone ~ migcertainty +lockdown1, England_2015_2020_aggregate_weekly_conscounts_migcertainty, sum) 
IR_pandemic_migcertainty_phone <- inner_join(conscount_pandemic_migcertainty_phone,pyears_pandemic_migcertainty_phone, by= c("migcertainty" = 'migcertainty', 'lockdown1' = 'lockdown1'))
function_to_calculate_IR_phone(IR_pandemic_migcertainty_phone, output_name = 'IR_pandemic_migcertainty_phone')
IR_pandemic_migcertainty_phone <- IR_pandemic_migcertainty_phone %>%
  filter(migcertainty %in% c('Definite', 'Probable')) %>%
  rename(migrant_status = migcertainty)

# Combine and convert to wide format for table
IRs_pandemic_phone <- bind_rows(IR_pandemic_migstatus_phone,
                              IR_pandemic_migcertainty_phone) %>%
  arrange(lockdown1) %>%
  dplyr::select(c(migrant_status, lockdown1, ir_ci)) %>%
  tidyr::pivot_wider(names_from = migrant_status, values_from = ir_ci) %>%
  mutate(type = 'Telephone consultations') %>%
  relocate(type)

# Combine all, F2F and phone cons
IRs_pandemic_combined <- bind_rows(IRs_pandemic, IRs_pandemic_F2F, IRs_pandemic_phone)

# Save as data file for tables 
save(IRs_pandemic_combined, file="filepath")

# 4b_Prepare data for interrupted time series ----

## Remove studyweek1 (as it's not a full week)
covid_weekly_conscounts_aggregate_England <- England_2015_2020_aggregate_weekly_conscounts %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
covid_weekly_conscounts_aggregate_England_ap <- covid_weekly_conscounts_aggregate_England
covid_weekly_conscounts_aggregate_England_ap$lockdown1[covid_weekly_conscounts_aggregate_England_ap$date %in% adjustment_period] <- NA

# 4c_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_phone')

# 4d_Segmented regression model: All consultations -----

## Fit model 
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))

round(ci.lin(full_model_nb_all,Exp=T),3)
extract_glm_results_allages(full_model_nb_all, 'all_incladjustmentperiod_nolaggedresiduals', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_nolaggedresiduals, "filepath" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

## Check autocorrelation
pacf(res, lag = 368)
acf(res,lag = 368)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_all <- full_join(covid_weekly_conscounts_aggregate_England_ap_all, lagres_timing, by = c('migrant_status','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr, 'all_incladjustmentperiod_laggedresiduals', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_laggedresiduals, "filepath" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

##  Make predictions
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_all
prediction_data$lagres <- 0
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_England_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_England_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>% 
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs
allcons_RRs <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
row.names(allcons_RRs) <- NULL
write.csv(allcons_RRs, "filepath")

# RERI - additive effect 
#get coefficients and RRs
exposure1 <- names(coef(final_model_nb_all_lr))[2]
exposure2 <- names(coef(final_model_nb_all_lr))[14]
interaction <- paste(exposure1, exposure2, sep = ':')

beta_lockdown <- coef(final_model_nb_all_lr)[exposure1]
beta_migrant_status <- coef(final_model_nb_all_lr)[exposure2]

beta_interaction <- coef(final_model_nb_all_lr)[interaction]
beta_interaction_CI <- exp(confint.default(final_model_nb_all_lr)[interaction, ])
RR_ref <- 1.0
RR_migrant_status <- exp(beta_migrant_status)
RR_lockdown <- exp(beta_lockdown)
RR_interaction <- exp(beta_interaction)
RR_migrant_status_lockdown <- exp(beta_migrant_status + beta_lockdown + beta_interaction)

# Calculate RERI
RERI <- RR_migrant_status_lockdown - RR_migrant_status - RR_lockdown + RR_ref

# Delta standard error
variance_covariance_matrix <-  as.matrix(vcov(final_model_nb_all_lr))
keepers <- c(exposure1, exposure2, interaction)
variance_covariance_matrix <- variance_covariance_matrix[keepers, keepers] # select variables of interest 

SE_RERI <- deltamethod(~ exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 1,
                       mean = c(beta_migrant_status, beta_lockdown, beta_interaction),
                       cov = variance_covariance_matrix)

# Delta 95% CIs
alpha = 1 - 0.95
z = qnorm(1 - alpha/2)
lower_CI <- RERI - z*SE_RERI 
upper_CI <- RERI + z*SE_RERI

# Combine and save RERI
RERI_all <- data.frame(analysis = 'Main analysis', RERI = RERI, CI = paste0(lower_CI, '-', upper_CI))
write.csv(RERI_all, "filepath")

# Predicted rates
# Aggregate IRs pre and during pandemic by migrant_status 

## All consultations
# migrant_status
pyears_pandemic_migstatus_pred <- aggregate(person_years ~ migrant_status + lockdown1, outcome_plot_nb, sum) 
pyears_pandemic_migstatus_pred <- pyears_pandemic_migstatus_pred %>% mutate(across(where(is.numeric), ~ round(.,)))
conscount_pandemic_migstatus_pred <- aggregate(pred ~ migrant_status +lockdown1, outcome_plot_nb, sum) 
IR_pandemic_migstatus_pred <- inner_join(conscount_pandemic_migstatus_pred,pyears_pandemic_migstatus_pred, by= c("migrant_status" = 'migrant_status', 'lockdown1' = 'lockdown1'))
function_to_calculate_IR(IR_pandemic_migstatus_pred, output_name = 'IR_pandemic_migstatus_pred')

output <- pois.exact(IR_pandemic_migstatus_pred$pred, IR_pandemic_migstatus_pred$person_years, conf.level=0.95)
output <- output %>% 
  rename(incidence_rate = rate)
output <- bind_cols(IR_pandemic_migstatus_pred, output) 
output <- output %>% mutate(across(where(is.numeric), ~ round(.,2)))
output$ci <- paste(output$lower, output$upper, sep ="-")
output$ir_ci <- paste0(output$incidence_rate, ' (', output$ci, ')')
output <-  dplyr::select(output, -c(conf.level, x, pt, ci))
write.csv(output, "filepath")



# Format for tables rmarkdown
allcons_RRs <- allcons_RRs %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(allcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, allcons) 
  
# 4e_Segmented regression model: Face-to-face consultations -----

# Fit model and calculate lagged residuals
full_model_nb_F2F <- glm.nb(facetoface ~ offset(log(person_years)) + studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F,Exp=T),3)

extract_glm_results_allages(full_model_nb_F2F, 'F2F_nolaggedresiduals', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_nolaggedresiduals, "filepath" )

lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res, lag = 364)
acf(res,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_F2F <- full_join(covid_weekly_conscounts_aggregate_England_ap_F2F, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_F2F_lr, 'F2F_laggedresiduals', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_laggedresiduals, "filepath" )


# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

prediction_data <- covid_weekly_conscounts_aggregate_England_ap_F2F
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_F2F <- nb_pred_F2F$fit
stbp_nb_F2F <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_F2F <- covid_weekly_conscounts_aggregate_England_ap_F2F %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb_F2F, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn_F2F <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_F2F <- bind_cols(stbp = stbp_nb_F2F, stbp_noLdn= stbp_nb_noLdn_F2F, pred = pred_values_F2F, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_F2F <- bind_cols(df_se_nb_F2F, covid_weekly_conscounts_aggregate_England_ap_F2F) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_F2F <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates_F2F %>%
  mutate(var = rownames(parameter_estimates_F2F)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status_F2F <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_F2F <- post_lockdown_migrant_status_F2F %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status_F2F)
row.names(F2Fcons_RRs) <- NULL
write.csv(F2Fcons_RRs, "filepath" )

# Format for tables rmarkdown
F2Fcons_RRs <- F2Fcons_RRs %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(F2Fcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, F2Fcons) 


# 4f_Segmented regression model: Phone consultations -----

# Fit model and calculate lagged residuals
full_model_nb_phone <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)

extract_glm_results_allages(full_model_nb_phone, 'phone_nolaggedresiduals', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_nolaggedresiduals, "filepath" )


lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res, lag = 364)
acf(res,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box')

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_phone <- full_join(covid_weekly_conscounts_aggregate_England_ap_phone, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_phone_lr, 'phone_laggedresiduals', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_laggedresiduals, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

# Create prediction data set 
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_phone
prediction_data$lagres <- 0 # set to 0 based on discussion with Ali and Amy
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_phone <- nb_pred_phone$fit
stbp_nb_phone <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_phone <- covid_weekly_conscounts_aggregate_England_ap_phone %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_phone_nolockdown <- predict(final_model_nb_phone_lr, newdata = datanew2_nb_phone, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_phone_noLockdown <-nb_pred_phone_nolockdown$fit
stbp_nb_noLdn_phone <- nb_pred_phone_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_phone <- bind_cols(stbp = stbp_nb_phone, stbp_noLdn= stbp_nb_noLdn_phone, pred = pred_values_phone, pred_noLdn = predicted_vals_nb_phone_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_phone <- bind_cols(df_se_nb_phone, covid_weekly_conscounts_aggregate_England_ap_phone) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_phone <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates_phone %>%
  mutate(var = rownames(parameter_estimates_phone)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

## Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status_phone <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_phone <- post_lockdown_migrant_status_phone %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

## Combine RRs and save
phonecons_RRs <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status_phone)
row.names(phonecons_RRs) <- NULL
write.csv(phonecons_RRs, "filepath" )

# Format for tables rmarkdown
phonecons_RRs <- phonecons_RRs %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(phonecons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, phonecons) 

# 4g_Combine plots to make a master plot -----

outcome_plot_nb$constype <- 'All consultations'
outcome_plot_nb_F2F$constype <- 'Face-to-face consultations'
outcome_plot_nb_phone$constype <- 'Telephone consultations'


all_outcome_plot_nb <- bind_rows(outcome_plot_nb, outcome_plot_nb_F2F, outcome_plot_nb_phone)

## Create dataframe for the vertical line lockdown labels (i.e. the adjustment-to-restrictions period in the COVID collateral paper - 8th March to the 28th March)
lockdowns <- data.frame(date = as_date(c('2020-03-08', '2020-03-29')))

combined_England_plot <- ggplot(filter(all_outcome_plot_nb, date <= as.Date('2020-11-30') & date > as.Date('2019-06-29')), 
                                          aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, color = migrant_status)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  geom_line(aes(y = incidence_rate), linetype = 2) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = NULL, date_breaks = '3 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) + 
  facet_wrap(~constype, ncol = 1, scales = 'free_y')

ggsave('filepath')

# 4h_Combine RRs for rmarkdown and save as Rdata file ----
RRs_England <- bind_cols(allcons_RRs, F2Fcons_RRs, phonecons_RRs) %>%
  dplyr::select(`var...1`, allcons, F2Fcons, phonecons) %>%
  rename('Variable' = `var...1`)
save(RRs_England, file = 'filepath')

# 5_SECONDARY ANALYSIS: age_subcohort ----

# 5a_Prepare data for interrupted time series -----

## Remove studyweek1 
covid_weekly_conscounts_aggregate_England_agesubcohort <- England_2015_2020_aggregate_weekly_conscounts_agesubcohort %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020
covid_weekly_conscounts_aggregate_England_agesubcohort_ap <- covid_weekly_conscounts_aggregate_England_agesubcohort
covid_weekly_conscounts_aggregate_England_agesubcohort_ap$lockdown1[covid_weekly_conscounts_aggregate_England_agesubcohort_ap$date %in% adjustment_period] <- NA

# 5b_Calculate IRs ----
# all data 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_agesubcohort, 'covid_weekly_conscounts_England_agesubcohort_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_agesubcohort, 'covid_weekly_conscounts_England_agesubcohort_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_agesubcohort, 'covid_weekly_conscounts_England_agesubcohort_phone')
# truncated data with lockdown1 set to NA during the adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone')

# 5c_Segmented regression model: All consultations ----

# 5d_0-15 years -----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all, age_subcohort == '0-15 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
round(ci.lin(full_model_nb_all,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Save df for combination plot
outcome_plot_nb_0to15years <- outcome_plot_nb
outcome_plot_nb_0to15years$age_subgroup <- '0-15 years'

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_0_15years <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
rownames(allcons_RRs_0_15years) <- NULL
write.csv(allcons_RRs_0_15years, "filepath" )

# Format for tables markdown
allcons_RRs_0_15years <- allcons_RRs_0_15years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group0to15 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group0to15) 

# 5e_16-24 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all, age_subcohort == '16-24 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
round(ci.lin(full_model_nb_all,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_all_lr)
round(ci.lin(final_model_nb_all_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Set up prediction data
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Save df for combination plot
outcome_plot_nb_16to24years <- outcome_plot_nb
outcome_plot_nb_16to24years$age_subgroup <- '16-24 years'

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_16_24years <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
rownames(allcons_RRs_16_24years) <- NULL
write.csv(allcons_RRs_16_24years, "filepath" )

# Format for tables markdown
allcons_RRs_16_24years <- allcons_RRs_16_24years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group16to24 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group16to24) 

# 5f_25-34 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all, age_subcohort == '25-34 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
round(ci.lin(full_model_nb_all,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_all_lr)
round(ci.lin(final_model_nb_all_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Save df for combination plot
outcome_plot_nb_25to34years <- outcome_plot_nb
outcome_plot_nb_25to34years$age_subgroup <- '25-34 years'

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_25_34years <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
rownames(allcons_RRs_25_34years) <- NULL
write.csv(allcons_RRs_25_34years, "filepath")

# Format for tables markdown
allcons_RRs_25_34years <- allcons_RRs_25_34years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group25to34 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group25to34) 

# 5g_35-49 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all, age_subcohort == '35-49 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
round(ci.lin(full_model_nb_all,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_all_lr)
round(ci.lin(final_model_nb_all_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Save df for combination plot
outcome_plot_nb_35to49years <- outcome_plot_nb
outcome_plot_nb_35to49years$age_subgroup <- '35-49 years'

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_35_49years <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
rownames(allcons_RRs_35_49years) <- NULL
write.csv(allcons_RRs_35_49years, "filepath")

# Format for tables markdown
allcons_RRs_35_49years <- allcons_RRs_35_49years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group35to49 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group35to49) 

# 5h_50-64 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all, age_subcohort == '50-64 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
round(ci.lin(full_model_nb_all,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_all_lr)
round(ci.lin(final_model_nb_all_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Save df for combination plot
outcome_plot_nb_50to64years <- outcome_plot_nb
outcome_plot_nb_50to64years$age_subgroup <- '50-64 years'

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_50_64years <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
rownames(allcons_RRs_50_64years) <- NULL
write.csv(allcons_RRs_50_64years, "filepath")

# Format for tables markdown
allcons_RRs_50_64years <- allcons_RRs_50_64years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group50to64 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group50to64) 

# 5i_65+ years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all, age_subcohort == '>=65 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
round(ci.lin(full_model_nb_all,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_all_lr)
round(ci.lin(final_model_nb_all_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

# Set up prediction dataset 
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Save df for combination plot
outcome_plot_nb_65plus <- outcome_plot_nb
outcome_plot_nb_65plus$age_subgroup <- '65 years and older'

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_65plusyears <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
rownames(allcons_RRs_65plusyears) <- NULL
write.csv(allcons_RRs_65plusyears, "filepath" )

# Format for tables markdown
allcons_RRs_65plusyears <- allcons_RRs_65plusyears %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group65plus = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group65plus) 

# 5j_Combine plots to make a master plot -------

all_outcome_plot_nb <- bind_rows(outcome_plot_nb_0to15years, outcome_plot_nb_16to24years,
                                 outcome_plot_nb_25to34years, outcome_plot_nb_35to49years,
                                 outcome_plot_nb_50to64years, outcome_plot_nb_65plus)

combined_age_subgroup_plot <- ggplot(filter(all_outcome_plot_nb, date > as.Date('2019-06-30') & date < as.Date('2020-11-29')), 
                                     aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, colour = migrant_status)) +
  geom_line()+
  geom_line(aes(y = incidence_rate), linetype = 'dotted') +
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  labs(y = 'Rate (per person-year)\n', x = NULL, title = NULL) + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = ' ', date_breaks = '6 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b')+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_wrap(~age_subgroup)
combined_age_subgroup_plot

ggsave('filepath')

# 5k_Combine RRs for rmarkdown and save as Rdata file ----
RRs_agesubcohort <- bind_cols(allcons_RRs_0_15years, allcons_RRs_16_24years, allcons_RRs_25_34years,
                              allcons_RRs_35_49years, allcons_RRs_50_64years, allcons_RRs_65plusyears) %>%
  dplyr::select(`var...1`, group0to15, group16to24, group25to34, group35to49, group50to64, group65plus) %>%
  rename('Variable' = `var...1`)
save(RRs_agesubcohort, file = 'filepath')

# 6_SECONDARY ANALYSIS: London -----


# 6a_Prepare data for interrupted time series ----

## Remove studyweek1 
covid_weekly_conscounts_aggregate_London <- London_2015_2020_aggregate_weekly_conscounts %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
covid_weekly_conscounts_aggregate_London_ap <- covid_weekly_conscounts_aggregate_London
covid_weekly_conscounts_aggregate_London_ap$lockdown1[covid_weekly_conscounts_aggregate_London_ap$date %in% adjustment_period] <- NA

# 6b_Calculate IRs ----
# Without accounting for an adjustment period
function_to_calculate_IR(covid_weekly_conscounts_aggregate_London, 'covid_weekly_conscounts_London_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_London, 'covid_weekly_conscounts_London_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_London, 'covid_weekly_conscounts_London_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_London_ap, 'covid_weekly_conscounts_aggregate_London_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_London_ap, 'covid_weekly_conscounts_aggregate_London_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_London_ap, 'covid_weekly_conscounts_aggregate_London_ap_phone')

lockdowns <- data.frame(date = as_date(c('2020-03-08', '2020-03-29')))

# 6c_Segmented regression model: All consultations -----

# Fit model and calculate lagged residuals
full_model_nb_all_london <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                                   data = filter(covid_weekly_conscounts_aggregate_London_ap_all, !is.na(lockdown1)))
summary(full_model_nb_all_london)
summary(full_model_nb_all_london$dispersion)
round(ci.lin(full_model_nb_all_london,Exp=T),3)

extract_glm_results_allages(full_model_nb_all_london, 'all_incladjustmentperiod_nolaggedresiduals_london', covid_weekly_conscounts_aggregate_London_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_nolaggedresiduals, "filepath" )

lagresiduals_m <- lag(residuals(full_model_nb_all_london)) %>% as.numeric()
res <- residuals(full_model_nb_all_london, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all_london, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_London_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_London_ap_all <- full_join(covid_weekly_conscounts_aggregate_London_ap_all, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr_london <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                         lagres,
                                       data =  filter(covid_weekly_conscounts_aggregate_London_ap_all, !is.na(lockdown1)))
summary(final_model_nb_all_lr_london)
round(ci.lin(final_model_nb_all_lr_london,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr_london, 'all_incladjustmentperiod_laggedresiduals_london', covid_weekly_conscounts_aggregate_London_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_laggedresiduals_london, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr_london, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr_london, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

prediction_data <- covid_weekly_conscounts_aggregate_London_ap_all
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr_london, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_London_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr_london, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions  
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_London_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr_london))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr_london))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr_london, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_london <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
row.names(allcons_RRs_london) <- NULL
write.csv(allcons_RRs_london, "filepath" )

# Format for tables rmarkdown
allcons_RRs_london <- allcons_RRs_london %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(allcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, allcons) 

# 6d_Segmented regression model: Face-to-face consultations ----------------------------------

# Fit model and calculate lagged residuals
full_model_nb_F2F_london <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                                   data = filter(covid_weekly_conscounts_aggregate_London_ap_F2F, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F_london,Exp=T),3)

extract_glm_results_allages(full_model_nb_F2F_london, 'F2F_nolaggedresiduals_london', covid_weekly_conscounts_aggregate_London_ap_F2F, 'lockdown1')
write.csv(F2F_nolaggedresiduals_london, "filepath" )

lagresiduals_m <- lag(residuals(full_model_nb_F2F_london)) %>% as.numeric()
res <- residuals(full_model_nb_F2F_london, type = 'deviance')

pacf(res, lag = 364)
acf(res,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F_london, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_London_ap_F2F, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_London_ap_F2F <- full_join(covid_weekly_conscounts_aggregate_London_ap_F2F, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr_london <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                         lagres,
                                       data =  filter(covid_weekly_conscounts_aggregate_London_ap_F2F, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr_london,Exp=T),3)
extract_glm_results_allages(final_model_nb_F2F_lr_london, 'F2F_laggedresiduals', covid_weekly_conscounts_aggregate_London_ap_F2F, 'lockdown1')
write.csv(F2F_laggedresiduals, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr_london, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr_london, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

prediction_data <- covid_weekly_conscounts_aggregate_London_ap_F2F
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr_london, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_F2F <- nb_pred_F2F$fit
stbp_nb_F2F <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_F2F <- covid_weekly_conscounts_aggregate_London_ap_F2F %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr_london, newdata = datanew2_nb_F2F, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn_F2F <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_F2F <- bind_cols(stbp = stbp_nb_F2F, stbp_noLdn= stbp_nb_noLdn_F2F, pred = pred_values_F2F, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_F2F <- bind_cols(df_se_nb_F2F, covid_weekly_conscounts_aggregate_London_ap_F2F) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates_F2F <- as.data.frame(ci.exp(final_model_nb_F2F_lr_london))
effect_of_lockdown_cons_F2F <- parameter_estimates_F2F %>%
  mutate(var = rownames(parameter_estimates_F2F)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr_london))
post_lockdown_migrant_status_F2F <- exp(estimable(final_model_nb_F2F_lr_london, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_F2F <- post_lockdown_migrant_status_F2F %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_london <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status_F2F)
row.names(F2Fcons_RRs_london) <- NULL
write.csv(F2Fcons_RRs_london, "filepath" )

# Format for tables rmarkdown
F2Fcons_RRs_london <- F2Fcons_RRs_london %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(F2Fcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, F2Fcons) 

# 6e_Segmented regression model: Phone consultations -----------------------------

# Fit model and calculate lagged residuals
full_model_nb_phone_london <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                                     data = filter(covid_weekly_conscounts_aggregate_London_ap_phone, !is.na(lockdown1)))
round(ci.lin(full_model_nb_phone_london,Exp=T),3)

extract_glm_results_allages(full_model_nb_phone_london, 'phone_nolaggedresiduals_london', covid_weekly_conscounts_aggregate_London_ap_phone, 'lockdown1')
write.csv(phone_nolaggedresiduals_london, "filepath" )

lagresiduals_1 <- lag(residuals(full_model_nb_phone_london),1) %>% as.numeric()
lagresiduals_2 <- lag(residuals(full_model_nb_phone_london),2) %>% as.numeric()
lagresiduals_3 <- lag(residuals(full_model_nb_phone_london),3) %>% as.numeric()
lagresiduals_4 <- lag(residuals(full_model_nb_phone_london),4) %>% as.numeric()


res <- residuals(full_model_nb_phone_london)

pacf(res) # first four spikes are above the significance line - so 4th order lagged residuals added 
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone_london, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_London_ap_phone, !is.na(lockdown1)),
                           'lagres1' = lagresiduals_1,
                           'lagres2' = lagresiduals_2,
                           'lagres3' = lagresiduals_3,
                           'lagres4' = lagresiduals_4,) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres1',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres2',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres3',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres4',~replace(., is.na(.), 0))

lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres1, lagres2, lagres3, lagres4, migrant_status))
covid_weekly_conscounts_aggregate_London_ap_phone <- full_join(covid_weekly_conscounts_aggregate_London_ap_phone, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr_london <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                           lagres1 +lagres2 + lagres3 + lagres4,
                                         data =  filter(covid_weekly_conscounts_aggregate_London_ap_phone, !is.na(lockdown1)))
round(ci.lin(final_model_nb_phone_lr_london,Exp=T),3)
extract_glm_results_allages(final_model_nb_phone_lr_london, 'phone_laggedresiduals_london', covid_weekly_conscounts_aggregate_London_ap_phone, 'lockdown1')
write.csv(phone_laggedresiduals_london, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr_london, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr_london, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Create prediction data set 
prediction_data <- covid_weekly_conscounts_aggregate_London_ap_phone
prediction_data$lagres1 <- 0
prediction_data$lagres2 <- 0
prediction_data$lagres3 <- 0
prediction_data$lagres4 <- 0

#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr_london, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_phone <- nb_pred_phone$fit
stbp_nb_phone <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_phone <- covid_weekly_conscounts_aggregate_London_ap_phone %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_phone_nolockdown <- predict(final_model_nb_phone_lr_london, newdata = datanew2_nb_phone, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_phone_noLockdown <-nb_pred_phone_nolockdown$fit
stbp_nb_noLdn_phone <- nb_pred_phone_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_phone <- bind_cols(stbp = stbp_nb_phone, stbp_noLdn= stbp_nb_noLdn_phone, pred = pred_values_phone, pred_noLdn = predicted_vals_nb_phone_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_phone <- bind_cols(df_se_nb_phone, covid_weekly_conscounts_aggregate_London_ap_phone) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates_phone <- as.data.frame(ci.exp(final_model_nb_phone_lr_london))
effect_of_lockdown_cons_phone <- parameter_estimates_phone %>%
  mutate(var = rownames(parameter_estimates_phone)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr_london))
post_lockdown_migrant_status_phone <- exp(estimable(final_model_nb_phone_lr_london, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_phone <- post_lockdown_migrant_status_phone %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
phonecons_RRs_london <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status_phone)
row.names(phonecons_RRs_london) <- NULL
write.csv(phonecons_RRs_london, "filepath" )

# Format for tables rmarkdown
phonecons_RRs_london <- phonecons_RRs_london %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(phonecons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, phonecons) 

# 6f_Combine plots to make a master plot -----

outcome_plot_nb$constype <- 'All consultations'
outcome_plot_nb_F2F$constype <- 'Face-to-face consultations'
outcome_plot_nb_phone$constype <- 'Telephone consultations'


all_outcome_plot_nb <- bind_rows(outcome_plot_nb, outcome_plot_nb_F2F, outcome_plot_nb_phone)

combined_London_plot <- ggplot(filter(all_outcome_plot_nb, date <= as.Date('2020-11-30') & date > as.Date('2019-06-29')), 
                               aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, color = migrant_status)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  geom_line(aes(y = incidence_rate), linetype = 2) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = NULL, date_breaks = '3 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) + 
  facet_wrap(~constype, ncol = 1, scales = 'free_y')


combined_London_plot

ggsave('filepath')

# 6g_Combine RRs for rmarkdown and save as Rdata file ----
RRs_London <- bind_cols(allcons_RRs_london, F2Fcons_RRs_london, phonecons_RRs_london) %>%
  dplyr::select(`var...1`, allcons, F2Fcons, phonecons) %>%
  rename('Variable' = `var...1`)
save(RRs_London, file = 'filepath')


# 7_SECONDARY ANALYSIS: Migrant_status and ethnicity (6 categories) ------

# 7a_Prepare data for interrupted time series ----

## Remove studyweek1 
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap$lockdown1[England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap$date %in% adjustment_period] <- NA

# 7b_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, 'England_weekly_conscounts_migstatus_ethnicity_all')
function_to_calculate_IR_F2F(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, 'England_weekly_conscounts_migstatus_ethnicity_F2F')
function_to_calculate_IR_phone(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, 'England_weekly_conscounts_migstatus_ethnicity_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap, 'England_weekly_conscounts_migstatus_ethnicity_ap_all')
function_to_calculate_IR_F2F(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap, 'England_weekly_conscounts_migstatus_ethnicity_ap_F2F')
function_to_calculate_IR_phone(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap, 'England_weekly_conscounts_migstatus_ethnicity_ap_phone')

# 7c_Plot IRs ----

## Create dataframe for the vertical line lockdown labels (i.e. the adjustment-to-restrictions period in the COVID collateral paper - 8th March to the 28th March)
lockdowns <- data.frame(date = as_date(c('2020-03-08', '2020-03-29')))

## All consultations
all_cons_IR_England_migstatus_ethnicity <- ggplot() + 
  geom_line(data = filter(England_weekly_conscounts_migstatus_ethnicity_all, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = migrant_status))+
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'All consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '2 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~ethnicat6)
all_cons_IR_England_migstatus_ethnicity
ggsave('filepath')

## Face-to-face consultations 
F2F_cons_IR_England_migstatus_ethnicity <-  ggplot() + 
  geom_line(data = filter(England_weekly_conscounts_migstatus_ethnicity_F2F, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = migrant_status))+
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Face-to-face consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '2 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~ethnicat6)
F2F_cons_IR_England_migstatus_ethnicity
ggsave('filepath')

## Phone consultations 
phone_cons_IR_England_migstatus_ethnicity <- ggplot() + 
  geom_line(data = filter(England_weekly_conscounts_migstatus_ethnicity_phone, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = migrant_status))+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Telephone consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '2 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_wrap(~ethnicat6)
phone_cons_IR_England_migstatus_ethnicity
ggsave('filepath')

# 7d_Segmented regression model: All consultations -----

## Fit model 
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years))+ studyweek + as.factor(studymonth) + migrant_status*ethnicat6*lockdown1,
                            data = filter(England_weekly_conscounts_migstatus_ethnicity_ap_all, !is.na(lockdown1)))
round(ci.lin(full_model_nb_all,Exp=T),3)
extract_glm_results_allages(full_model_nb_all, 'migstatus_ethnicity_all_incladjustmentperiod_nolaggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_all, 'lockdown1')
write.csv(migstatus_ethnicity_all_incladjustmentperiod_nolaggedresiduals, "results/01_Consultations/its_analysis/migstatus_ethnicity_all_incladjustmentperiod_nolaggedresiduals.csv" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

## Check autocorrelation
pacf(res, lag = 368)
acf(res,lag = 368)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(England_weekly_conscounts_migstatus_ethnicity_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status, ethnicat6))
England_weekly_conscounts_migstatus_ethnicity_ap_all <- full_join(England_weekly_conscounts_migstatus_ethnicity_ap_all, lagres_timing, by = c('migrant_status', 'ethnicat6','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + studyweek+ as.factor(studymonth)+ migrant_status*ethnicat6*lockdown1 + 
                                  lagres,
                                data =  filter(England_weekly_conscounts_migstatus_ethnicity_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr, 'migstatus_ethnicity_all_incladjustmentperiod_laggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_all, 'lockdown1')
write.csv(migstatus_ethnicity_all_incladjustmentperiod_laggedresiduals, "filepath" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'response')
pacf(res2)
acf(res2)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

##  Make predictions
prediction_data <- England_weekly_conscounts_migstatus_ethnicity_ap_all
prediction_data$lagres <- 0
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- England_weekly_conscounts_migstatus_ethnicity_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, England_weekly_conscounts_migstatus_ethnicity_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)


# 7e_Plot the results and save -----

all_cons_migstatus_ethnicity_nb <- ggplot(filter(outcome_plot_nb, date <= as.Date('2020-11-30') & date > as.Date('2019-06-29')), 
                                            aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, color = migrant_status)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  geom_line(aes(y = incidence_rate), linetype = 2) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = '', date_breaks = '3 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) + 
  facet_wrap(~ethnicat6, ncol = 4)


all_cons_migstatus_ethnicity_nb
ggsave('filepath')

# 7f_Get rate ratios ------

## Get RR migrants of different ethnicities versus the reference group (white British non-migrants)
names(coef(final_model_nb_all_lr))

### White British migrants vs WBNM
post_lockdown_WBvsWBNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1), conf=.95))
post_lockdown_WBvsWBNM <- post_lockdown_WBvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1')
### White non-British migrants vs WBNM
post_lockdown_WNBvsWBNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'ethnicat6White non-British:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6White non-British:lockdown1' = 1), conf=.95))
post_lockdown_WNBvsWBNM <- post_lockdown_WNBvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6White non-British:lockdown1 + migrant_statusMigrant:ethnicat6White non-British:lockdown1')
### Mixed migrants vs WBNM
post_lockdown_MixedvsWBNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                              'ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1,
                                                              'migrant_statusMigrant:ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1), conf=.95))
post_lockdown_MixedvsWBNM <- post_lockdown_MixedvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Mixed:lockdown1 + migrant_statusMigrant:ethnicat6Mixed:lockdown1')
### Asian migrants vs WBNM
post_lockdown_AsianvsWBNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'ethnicat6Asian/Asian British:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Asian/Asian British:lockdown1' = 1), conf=.95))
post_lockdown_AsianvsWBNM <- post_lockdown_AsianvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Asian:lockdown1 + migrant_statusMigrant:ethnicat6Asian:lockdown1')
### Black migrants vs WBNM
post_lockdown_BlackvsWBNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                              'ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1,
                                                              'migrant_statusMigrant:ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1), conf=.95))
post_lockdown_BlackvsWBNM <- post_lockdown_BlackvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Black:lockdown1 + migrant_statusMigrant:ethnicat6Black:lockdown1')
### Other migrants vs WBNM
post_lockdown_OthervsWBNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'ethnicat6Other ethnic group:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Other ethnic group:lockdown1' = 1), conf=.95))
post_lockdown_OthervsWBNM <- post_lockdown_OthervsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Other:lockdown1 + migrant_statusMigrant:ethnicat6Other:lockdown1')
### Unknown migrants vs WBNM
post_lockdown_UnknownvsWBNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                'ethnicat6Unknown:lockdown1' = 1,
                                                                'migrant_statusMigrant:ethnicat6Unknown:lockdown1' = 1), conf=.95))
post_lockdown_UnknownvsWBNM <- post_lockdown_UnknownvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Unknown:lockdown1 + migrant_statusMigrant:ethnicat6Unknown:lockdown1')

## Combine RRs and save
allcons_migstatus_ethnicity_RRs_vsWBNMs <- bind_rows(post_lockdown_WBvsWBNM, post_lockdown_WNBvsWBNM,post_lockdown_MixedvsWBNM,
                                             post_lockdown_AsianvsWBNM, post_lockdown_BlackvsWBNM, post_lockdown_OthervsWBNM,
                                             post_lockdown_UnknownvsWBNM)
row.names(allcons_migstatus_ethnicity_RRs_vsWBNMs) <- NULL
write.csv(allcons_migstatus_ethnicity_RRs_vsWBNMs, "filepath" )

## Get RR for migrants vs non-migrants of different ethnicities 
names(coef(final_model_nb_all_lr))

### White British migrants vs non-migrants
post_lockdown_WBvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1), conf=.95))
post_lockdown_WBvsNM <- post_lockdown_WBvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1')
### White non-British migrants vs non-migrants
post_lockdown_WNBvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6White non-British:lockdown1' = 1), conf=.95))
post_lockdown_WNBvsNM <- post_lockdown_WNBvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
          'upper' = Upper.CI,
          'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6White non-British:lockdown1')
### Mixed migrants vs non-migrants
post_lockdown_MixedvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1), conf=.95))
post_lockdown_MixedvsNM <- post_lockdown_MixedvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
          'upper' = Upper.CI,
          'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Mixed:lockdown1')
### Asian migrants vs non-migrants
post_lockdown_AsianvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Asian/Asian British:lockdown1' = 1), conf=.95))
post_lockdown_AsianvsNM <- post_lockdown_AsianvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Asian:lockdown1')
### Black migrants vs non-migrants
post_lockdown_BlackvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1), conf=.95))
post_lockdown_BlackvsNM <- post_lockdown_BlackvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Black:lockdown1')
### Other migrants vs non-migrants
post_lockdown_OthervsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Other ethnic group:lockdown1' = 1), conf=.95))
post_lockdown_OthervsNM <- post_lockdown_OthervsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Other:lockdown1')
### Unknown migrants vs non-migrants
post_lockdown_UnknownvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicat6Unknown:lockdown1' = 1), conf=.95))
post_lockdown_UnknownvsNM <- post_lockdown_UnknownvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Unknown:lockdown1')

# Combine RRs and save
allcons_migstatus_ethnicity_RRs_vsNMs <- bind_rows(post_lockdown_WBvsNM, post_lockdown_WNBvsNM,post_lockdown_MixedvsNM,
                                             post_lockdown_AsianvsNM, post_lockdown_BlackvsNM, post_lockdown_OthervsNM,
                                             post_lockdown_UnknownvsNM)
row.names(allcons_migstatus_ethnicity_RRs_vsNMs) <- NULL
write.csv(allcons_migstatus_ethnicity_RRs_vsNMs, "filepath")

# Format for tables rmarkdown and save

RRs_England_migstatus_ethnicity <- allcons_migstatus_ethnicity_RRs_vsNMs %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(RR95CI = paste0(estimate, ' (', lower, '-', upper)) %>%
  dplyr::select(var, RR95CI)
save(RRs_England_migstatus_ethnicity, file = 'filepath')

# 8_SECONDARY ANALYSIS: Migrant_status and ethnicity (18 categories) ----
# 8a_Prepare data for interrupted time series ----

## Remove studyweek1 
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity18cat %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap$lockdown1[England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap$date %in% adjustment_period] <- NA

# 8b_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, 'England_weekly_conscounts_migstatus_ethnicity_all')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap, 'England_weekly_conscounts_migstatus_ethnicity_ap_all')


# 8c_Segmented regression model: All consultations -----

## Fit model 
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years))+ studyweek + as.factor(studymonth) + migrant_status*ethnicat*lockdown1,
                            data = filter(England_weekly_conscounts_migstatus_ethnicity_ap_all, !is.na(lockdown1)))
round(ci.lin(full_model_nb_all,Exp=T),3)
#extract_glm_results_allages(full_model_nb_all, 'migstatus_ethnicity_all_incladjustmentperiod_nolaggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_all, 'lockdown1')
#write.csv(migstatus_ethnicity_all_incladjustmentperiod_nolaggedresiduals, "filepath")

res <- residuals(full_model_nb_all, type = 'deviance')

## Check autocorrelation
pacf(res)
acf(res)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

## Calculate lagged residuals 

lagresiduals_1 <- lag(residuals(full_model_nb_all),1) %>% as.numeric()
lagresiduals_2 <- lag(residuals(full_model_nb_all),2) %>% as.numeric()
lagresiduals_3 <- lag(residuals(full_model_nb_all),3) %>% as.numeric()
lagresiduals_4 <- lag(residuals(full_model_nb_all),4) %>% as.numeric()

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(England_weekly_conscounts_migstatus_ethnicity_ap_all, !is.na(lockdown1)),
                           'lagres1' = lagresiduals_1,
                           'lagres2' = lagresiduals_2,
                           'lagres3' = lagresiduals_3,
                           'lagres4' = lagresiduals_4,) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres1',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres2',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres3',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres4',~replace(., is.na(.), 0))

lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres1, lagres2, lagres3, lagres4, migrant_status, ethnicat))
England_weekly_conscounts_migstatus_ethnicity_ap_all <- full_join(England_weekly_conscounts_migstatus_ethnicity_ap_all, lagres_timing, by = c('migrant_status','studyweek', 'ethnicat'))

## Re-run final model (including lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + studyweek+ as.factor(studymonth)+ migrant_status*ethnicat*lockdown1 + 
                                  lagres1 + lagres2 + lagres3 + lagres4,
                                data =  filter(England_weekly_conscounts_migstatus_ethnicity_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
#extract_glm_results_allages(final_model_nb_all_lr, 'migstatus_ethnicity_all_incladjustmentperiod_laggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_all, 'lockdown1')
#write.csv(migstatus_ethnicity_all_incladjustmentperiod_laggedresiduals, "filepath")

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'response')
pacf(res2)
acf(res2)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

##  Make predictions
prediction_data <- England_weekly_conscounts_migstatus_ethnicity_ap_all
prediction_data$lagres1 <- 0
prediction_data$lagres2 <- 0
prediction_data$lagres3 <- 0
prediction_data$lagres4 <- 0
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- England_weekly_conscounts_migstatus_ethnicity_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, England_weekly_conscounts_migstatus_ethnicity_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)


# 8d_Plot the results and save----

library(stringr)
var_width = 20
outcome_plot_nb <- mutate(outcome_plot_nb, tidy_ethnicat_name = str_wrap(ethnicat, width = var_width))
outcome_plot_nb$tidy_ethnicat6_name <- factor(outcome_plot_nb$tidy_ethnicat6_name,levels = c("White British", "White non-British", "Mixed/Multiple\nethnic groups", 
                                                                                             "Asian/Asian British", "Black/African/\nCaribbean/Black\nBritish", "Other ethnic group", 'Unknown'))

all_cons_migstatus_ethnicity_nb <- ggplot(filter(outcome_plot_nb, date <= as.Date('2020-11-30') & date > as.Date('2019-06-29')), 
                                          aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, color = migrant_status)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  geom_line(aes(y = incidence_rate), linetype = 2) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '3 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) + 
  facet_wrap(~ethnicat, scales = 'free_y')


all_cons_migstatus_ethnicity_nb
ggsave('filepath')

# 8e_Get rate ratios -----

## Get RR for migrants vs non-migrants of the same ethnicity 
names(coef(final_model_nb_all_lr))

### British migrants vs non-migrants
post_lockdown_BvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1), conf=.95))
post_lockdown_BvsNM <- post_lockdown_BvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1')

### Irish migrants vs non-migrants
post_lockdown_IrishvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicatIrish:lockdown1' = 1), conf=.95))
post_lockdown_IrishvsNM <- post_lockdown_IrishvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatIrish:lockdown1')

### Gypsy or Irish Traveller migrants vs non-migrants
post_lockdown_GITvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                'migrant_statusMigrant:ethnicatGypsy or Irish Traveller:lockdown1' = 1), conf=.95))
post_lockdown_GITvsNM <- post_lockdown_GITvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatGypsy or Irish Traveller:lockdown1')

### Other White migrants vs non-migrants
post_lockdown_otherWhitevsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                       'migrant_statusMigrant:ethnicatOther White:lockdown1' = 1), conf=.95))
post_lockdown_otherWhitevsNM <- post_lockdown_otherWhitevsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatOther White:lockdown1')

### Mixed White and Black Caribbean migrants vs non-migrants
post_lockdown_mixedWhiteBlackCaribbeanvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                                     'migrant_statusMigrant:ethnicatMixed White and Black Caribbean:lockdown1' = 1), conf=.95))
post_lockdown_mixedWhiteBlackCaribbeanvsNM <- post_lockdown_mixedWhiteBlackCaribbeanvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatMixed White and Black Caribbean:lockdown1')

### Mixed White and Black African migrants vs non-migrants
post_lockdown_mixedWhiteBlackAfricanvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                                   'migrant_statusMigrant:ethnicatMixed White and Black African:lockdown1' = 1), conf=.95))
post_lockdown_mixedWhiteBlackAfricanvsNM <- post_lockdown_mixedWhiteBlackAfricanvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatMixed White and Black African:lockdown1')

### Mixed White and Asian migrants vs non-migrants
post_lockdown_mixedWhiteAsianvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                            'migrant_statusMigrant:ethnicatMixed White and Asian:lockdown1' = 1), conf=.95))
post_lockdown_mixedWhiteAsianvsNM <- post_lockdown_mixedWhiteAsianvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatMixed White and Asian:lockdown1')

### Other Mixed migrants vs non-migrants
post_lockdown_otherMixedvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                       'migrant_statusMigrant:ethnicatOther Mixed:lockdown1' = 1), conf=.95))
post_lockdown_otherMixedvsNM <- post_lockdown_otherMixedvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatOther Mixed:lockdown1')

### Indian migrants vs non-migrants
post_lockdown_indianvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                   'migrant_statusMigrant:ethnicatIndian:lockdown1' = 1), conf=.95))
post_lockdown_indianvsNM <- post_lockdown_indianvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatIndian:lockdown1')

### Pakistani migrants vs non-migrants
post_lockdown_pakistanivsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicatPakistani:lockdown1' = 1), conf=.95))
post_lockdown_pakistanivsNM <- post_lockdown_pakistanivsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatPakistani:lockdown1')

### Bangladesh migrants vs non-migrants
post_lockdown_bangladeshivsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                        'migrant_statusMigrant:ethnicatBangladeshi:lockdown1' = 1), conf=.95))
post_lockdown_bangladeshivsNM <- post_lockdown_bangladeshivsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatBangladeshi:lockdown1')

### Chinese migrants vs non-migrants
post_lockdown_chinesevsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicatChinese:lockdown1' = 1), conf=.95))
post_lockdown_chinesevsNM<- post_lockdown_chinesevsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatChinese:lockdown1')

### Other Asian migrants vs non-migrants
post_lockdown_otherAsianvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                       'migrant_statusMigrant:ethnicatOther Asian:lockdown1' = 1), conf=.95))
post_lockdown_otherAsianvsNM <- post_lockdown_otherAsianvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatOther Asian:lockdown1')

### African migrants vs non-migrants
post_lockdown_AfricanvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicatAfrican:lockdown1' = 1), conf=.95))
post_lockdown_AfricanvsNM <- post_lockdown_AfricanvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatAfrican:lockdown1')

### Caribbean migrants vs non-migrants
post_lockdown_CaribbeanvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicatCaribbean:lockdown1' = 1), conf=.95))
post_lockdown_CaribbeanvsNM <- post_lockdown_CaribbeanvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatCaribbean:lockdown1')

### Other black migrants vs non-migrants
post_lockdown_otherBlackvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                       'migrant_statusMigrant:ethnicatOther Black:lockdown1' = 1), conf=.95))
post_lockdown_otherBlackvsNM <- post_lockdown_otherBlackvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatOther Black:lockdown1')

### Arab migrants vs non-migrants
post_lockdown_ArabvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                 'migrant_statusMigrant:ethnicatArab:lockdown1' = 1), conf=.95))
post_lockdown_ArabvsNM <- post_lockdown_ArabvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatArab:lockdown1')

### Any other ethnic group migrants vs non-migrants
post_lockdown_othervsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicatAny other ethnic group:lockdown1' = 1), conf=.95))
post_lockdown_othervsNM <- post_lockdown_othervsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatAny other ethnic group:lockdown1')

### Unknown migrants vs non-migrants
post_lockdown_unknownvsNM <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicatUnknown:lockdown1' = 1), conf=.95))
post_lockdown_unknownvsNM <- post_lockdown_unknownvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicatUnknown:lockdown1')

# Combine RRs and save
allcons_migstatus_ethnicity18cat_RRs_vsNMs <- bind_rows(post_lockdown_BvsNM, post_lockdown_IrishvsNM,
                                                        post_lockdown_GITvsNM, post_lockdown_otherWhitevsNM, 
                                                        post_lockdown_mixedWhiteBlackCaribbeanvsNM, post_lockdown_mixedWhiteBlackAfricanvsNM,
                                                        post_lockdown_mixedWhiteAsianvsNM, post_lockdown_otherMixedvsNM,
                                                        post_lockdown_indianvsNM, post_lockdown_pakistanivsNM,
                                                        post_lockdown_bangladeshivsNM, post_lockdown_chinesevsNM,
                                                        post_lockdown_otherAsianvsNM, post_lockdown_AfricanvsNM,
                                                        post_lockdown_CaribbeanvsNM, post_lockdown_otherBlackvsNM,
                                                        post_lockdown_ArabvsNM, post_lockdown_othervsNM,
                                                        post_lockdown_unknownvsNM)

row.names(allcons_migstatus_ethnicity18cat_RRs_vsNMs) <- NULL
write.csv(allcons_migstatus_ethnicity18cat_RRs_vsNMs, "filepath")

# Format for tables rmarkdown and save

RRs_England_migstatus_ethnicity18cat <- allcons_migstatus_ethnicity18cat_RRs_vsNMs %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(RR95CI = paste0(estimate, ' (', lower, '-', upper,')')) %>%
  dplyr::select(var, RR95CI)
save(RRs_England_migstatus_ethnicity18cat, file = 'filepath')

# 9_SENSITIVITY ANALYSIS: Migcertainty ------

# 9a_Prepare data for interrupted time series ----

## Remove studyweek1 and renumber studyweeks 
covid_weekly_conscounts_aggregate_England_migcertainty <- England_2015_2020_aggregate_weekly_conscounts_migcertainty %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
covid_weekly_conscounts_aggregate_England_migcertainty_ap <- covid_weekly_conscounts_aggregate_England_migcertainty
covid_weekly_conscounts_aggregate_England_migcertainty_ap$lockdown1[covid_weekly_conscounts_aggregate_England_migcertainty_ap$date %in% adjustment_period] <- NA

# 9b_Calculate IRs ----
# all data 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_migcertainty, 'covid_weekly_conscounts_migcertainty_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_migcertainty, 'covid_weekly_conscounts_migcertainty_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_migcertainty, 'covid_weekly_conscounts_migcertainty_phone')
# truncated data with lockdown1 set to NA during the adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_migcertainty_ap, 'covid_weekly_conscounts_aggregate_England_migcertainty_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_migcertainty_ap, 'covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_migcertainty_ap, 'covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone')

# 9c_Plot IRs  ----

lockdowns <- data.frame(date = as_date(c('2020-03-08', '2020-03-29')))

# All consultations
all_cons_IR_England_migcertainty <- ggplot() + 
  geom_line(data = filter(covid_weekly_conscounts_migcertainty_all, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = migcertainty))+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'All consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))

all_cons_IR_England_migcertainty
#ggsave('filepath')

# Face-to-face consultations 
F2F_cons_IR_England_migcertainty <- ggplot() + 
  geom_line(data = filter(covid_weekly_conscounts_migcertainty_F2F, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = migcertainty))+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Face-to-face consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))

F2F_cons_IR_England_migcertainty
#ggsave('filepath')

# Phone consultations 
phone_cons_IR_England_migcertainty <- ggplot() + 
  geom_line(data = filter(covid_weekly_conscounts_migcertainty_phone, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = migcertainty))+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Telephone consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  ylim(0,1)

phone_cons_IR_England_migcertainty
#ggsave('filepath')

# 9d_Segmented regression model: All consultations -------

# Fit model and calculate lagged residuals
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_all, !is.na(lockdown1)))
round(ci.lin(full_model_nb_all,Exp=T),3)

extract_glm_results_allages(full_model_nb_all, 'all_incladjustmentperiod_nolaggedresiduals_migcertainty', covid_weekly_conscounts_aggregate_England_migcertainty_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_nolaggedresiduals_migcertainty, "filepath" )

lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migcertainty))
covid_weekly_conscounts_aggregate_England_migcertainty_ap_all <- full_join(covid_weekly_conscounts_aggregate_England_migcertainty_ap_all, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr, 'all_incladjustmentperiod_laggedresiduals_migcertainty', covid_weekly_conscounts_aggregate_England_migcertainty_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_laggedresiduals_migcertainty, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

prediction_data <- covid_weekly_conscounts_aggregate_England_migcertainty_ap_all
prediction_data$lagres <- 0
#  Make predictions
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_England_migcertainty_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_England_migcertainty_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all_migcertainty <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migcertaintyDefinite', 'lockdown1:migcertaintyDefinite', 
                    'migcertaintyProbable', 'lockdown1:migcertaintyProbable')) 


# Get RR for the effect of migcertainty on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
# Definite migrants
post_lockdown_migcertaintyDefinite <- exp(estimable(final_model_nb_all_lr, c('migcertaintyDefinite' = 1, 'lockdown1:migcertaintyDefinite' = 1), conf=.95))
post_lockdown_migcertaintyDefinite <- post_lockdown_migcertaintyDefinite %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migcertaintyDefinite+lockdown1:migcertaintyDefinite')
# Probable migrants
post_lockdown_migcertaintyProbable <- exp(estimable(final_model_nb_all_lr, c('migcertaintyProbable' = 1, 'lockdown1:migcertaintyProbable' = 1), conf=.95))
post_lockdown_migcertaintyProbable <- post_lockdown_migcertaintyProbable %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migcertaintyProbable+lockdown1:migcertaintyProbable')

# Combine RRs and save
allcons_RRs_migcertainty <- bind_rows(effect_of_lockdown_cons_all_migcertainty,
                                      post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(allcons_RRs_migcertainty) = NULL
write.csv(allcons_RRs_migcertainty, "filepath")

# Format for tables rmarkdown
allcons_RRs_migcertainty <- allcons_RRs_migcertainty %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(allcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, allcons) 

# 9e_Segmented regression model: F2F consultations -----

# Fit model and calculate lagged residuals
full_model_nb_F2F <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F,Exp=T),3)

extract_glm_results_allages(full_model_nb_F2F, 'F2F_incladjustmentperiod_nolaggedresiduals_migcertainty', covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F, 'lockdown1')
write.csv(F2F_incladjustmentperiod_nolaggedresiduals_migcertainty, "filepath" )

res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res) # lags 1, 3, 4 and 5 are significant 
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Create fourth-order lagged residuals and merge lag residuals with dataframe

lagresiduals_1 <- lag(residuals(full_model_nb_F2F),1) %>% as.numeric()
lagresiduals_3 <- lag(residuals(full_model_nb_F2F),3) %>% as.numeric()
lagresiduals_4 <- lag(residuals(full_model_nb_F2F),4) %>% as.numeric()
lagresiduals_5 <- lag(residuals(full_model_nb_F2F),5) %>% as.numeric()

lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F, !is.na(lockdown1)),
                           'lagres1' = lagresiduals_1,
                           'lagres3' = lagresiduals_3,
                           'lagres4' = lagresiduals_4,
                           'lagres5' = lagresiduals_5,) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres1',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres3',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres4',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres5',~replace(., is.na(.), 0))
  

lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres1, lagres3, lagres4, lagres5, migcertainty))
covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F <- full_join(covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
                                  lagres1 + lagres3 + lagres4 + lagres5,
                                data =  filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_F2F_lr, 'F2F_incladjustmentperiod_laggedresiduals_migcertainty', covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F, 'lockdown1')
write.csv(F2F_incladjustmentperiod_laggedresiduals_migcertainty, "filepath")

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2) # 
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')

prediction_data <- covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F
prediction_data$lagres1 <- 0
prediction_data$lagres3 <- 0
prediction_data$lagres4 <- 0
prediction_data$lagres5 <- 0

#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb_F2F <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_England_migcertainty_ap_F2F) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

F2F_cons_nb_migcertainty <- ggplot(filter(outcome_plot_nb_F2F, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), 
                                   aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migcertainty, colour = migcertainty)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  geom_line(aes(y = incidence_rate), linetype = 2) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Face-to-face consultations') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) +
  ylim(0,7)

F2F_cons_nb_migcertainty
#ggsave('filepath')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F_migcertainty <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migcertaintyDefinite', 'lockdown1:migcertaintyDefinite', 
                    'migcertaintyProbable', 'lockdown1:migcertaintyProbable')) 

# Get RR for the effect of migcertainty on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
# Definite migrants
post_lockdown_migcertaintyDefinite <- exp(estimable(final_model_nb_F2F_lr, c('migcertaintyDefinite' = 1, 'lockdown1:migcertaintyDefinite' = 1), conf=.95))
post_lockdown_migcertaintyDefinite <- post_lockdown_migcertaintyDefinite %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migcertaintyDefinite+lockdown1:migcertaintyDefinite')
# Probable migrants
post_lockdown_migcertaintyProbable <- exp(estimable(final_model_nb_F2F_lr, c('migcertaintyProbable' = 1, 'lockdown1:migcertaintyProbable' = 1), conf=.95))
post_lockdown_migcertaintyProbable <- post_lockdown_migcertaintyProbable %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migcertaintyProbable+lockdown1:migcertaintyProbable')

# Combine RRs and save
F2Fcons_RRs_migcertainty <- bind_rows(effect_of_lockdown_cons_F2F_migcertainty,
                                      post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(F2Fcons_RRs_migcertainty) = NULL
write.csv(F2Fcons_RRs_migcertainty, "filepath" )

# Format for tables rmarkdown
F2Fcons_RRs_migcertainty <- F2Fcons_RRs_migcertainty %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(F2Fcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, F2Fcons) 

# 9f_Segmented regression model: Phone consultations -----

# Fit model and calculate lagged residuals
full_model_nb_phone <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone, !is.na(lockdown1)))
round(ci.lin(full_model_nb_phone,Exp=T),3)

extract_glm_results_allages(full_model_nb_phone, 'phone_incladjustmentperiod_nolaggedresiduals_migcertainty', covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone, 'lockdown1')
write.csv(phone_incladjustmentperiod_nolaggedresiduals_migcertainty, "filepath" )

res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res) # lags 1-5 and 7 are significant 
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Create lagged residuals and merge lag residuals with dataframe

lagresiduals_1 <- lag(residuals(full_model_nb_phone),1) %>% as.numeric()
lagresiduals_2 <- lag(residuals(full_model_nb_phone),2) %>% as.numeric()
lagresiduals_3 <- lag(residuals(full_model_nb_phone),3) %>% as.numeric()
lagresiduals_4 <- lag(residuals(full_model_nb_phone),4) %>% as.numeric()
lagresiduals_5 <- lag(residuals(full_model_nb_phone),5) %>% as.numeric()
lagresiduals_7 <- lag(residuals(full_model_nb_phone),7) %>% as.numeric()

lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone, !is.na(lockdown1)),
                           'lagres1' = lagresiduals_1,
                           'lagres2' = lagresiduals_2,
                           'lagres3' = lagresiduals_3,
                           'lagres4' = lagresiduals_4,
                           'lagres5' = lagresiduals_5,
                           'lagres7' = lagresiduals_7,) 

lagres_timing <- lagres_timing %>%
  mutate_at('lagres1',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres2',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres3',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres4',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres5',~replace(., is.na(.), 0)) %>%
  mutate_at('lagres7',~replace(., is.na(.), 0)) 


lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres1, lagres2, lagres3, lagres4, lagres5, lagres7, migcertainty))
covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone <- full_join(covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
                                  lagres1 + lagres2 + lagres3 + lagres4 + lagres5 + lagres7,
                                data =  filter(covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone, !is.na(lockdown1)))
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_phone_lr, 'phone_incladjustmentperiod_laggedresiduals_migcertainty', covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone, 'lockdown1')
write.csv(phone_incladjustmentperiod_laggedresiduals_migcertainty, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')

# Set up prediction data frame
prediction_data <- covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone
prediction_data$lagres1 <- 0
prediction_data$lagres2 <- 0
prediction_data$lagres3 <- 0
prediction_data$lagres4 <- 0
prediction_data$lagres5 <- 0
prediction_data$lagres7 <- 0

#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_phone_nolockdown <- predict(final_model_nb_phone_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_phone_noLockdown <-nb_pred_phone_nolockdown$fit
stbp_nb_noLdn <- nb_pred_phone_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_phone_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn),
  )

# Combine data set and predictions
outcome_plot_nb_phone <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_England_migcertainty_ap_phone) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

## Plot results

phone_cons_nb_migcertainty <- ggplot(filter(outcome_plot_nb_phone, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), 
                                   aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migcertainty, colour = migcertainty)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  geom_line(aes(y = incidence_rate), linetype = 2) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Telephone consultations') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) +
  ylim(0,1)

phone_cons_nb_migcertainty
#ggsave('filepath')

## Get rate ratios 

### Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone_migcertainty <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migcertaintyDefinite', 'lockdown1:migcertaintyDefinite', 
                    'migcertaintyProbable', 'lockdown1:migcertaintyProbable')) 

### Get RR for the effect of migcertainty on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
#### Definite migrants
post_lockdown_migcertaintyDefinite <- exp(estimable(final_model_nb_phone_lr, c('migcertaintyDefinite' = 1, 'lockdown1:migcertaintyDefinite' = 1), conf=.95))
post_lockdown_migcertaintyDefinite <- post_lockdown_migcertaintyDefinite %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migcertaintyDefinite+lockdown1:migcertaintyDefinite')
#### Probable migrants
post_lockdown_migcertaintyProbable <- exp(estimable(final_model_nb_phone_lr, c('migcertaintyProbable' = 1, 'lockdown1:migcertaintyProbable' = 1), conf=.95))
post_lockdown_migcertaintyProbable <- post_lockdown_migcertaintyProbable %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migcertaintyProbable+lockdown1:migcertaintyProbable')

### Combine RRs and save
phonecons_RRs_migcertainty <- bind_rows(effect_of_lockdown_cons_phone_migcertainty,
                                      post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(phonecons_RRs_migcertainty) = NULL
write.csv(phonecons_RRs_migcertainty, "filepath" )

# Format for tables markdown
phonecons_RRs_migcertainty <- phonecons_RRs_migcertainty %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(phonecons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, phonecons) 

# 9g_Combine plots to make a master plot -----

outcome_plot_nb$constype <- 'All consultations'
outcome_plot_nb_F2F$constype <- 'Face-to-face consultations'
outcome_plot_nb_phone$constype <- 'Telephone consultations'


all_outcome_plot_nb <- bind_rows(outcome_plot_nb, outcome_plot_nb_F2F, outcome_plot_nb_phone)

combined_migcertainty_plot <- ggplot(filter(all_outcome_plot_nb, date > as.Date('2019-06-30') & date < as.Date('2020-11-29')), 
                               aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migcertainty, colour = migcertainty)) +
  geom_line()+
  geom_line(aes(y = incidence_rate), linetype = 2) +
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  labs(y = 'Rate (per person-year)', x = NULL, title = NULL) + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = NULL, date_breaks = '3 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_wrap(~constype, ncol = 1, scales = 'free_y')

combined_migcertainty_plot

ggsave('filepath')

# 9h_Combine RRs for rmarkdown and save as Rdata file ----
RRs_migcertainty <- bind_cols(allcons_RRs_migcertainty, F2Fcons_RRs_migcertainty, phonecons_RRs_migcertainty) %>%
  dplyr::select(`var...1`, allcons, F2Fcons, phonecons) %>%
  rename('Variable' = `var...1`)
save(RRs_migcertainty, file = 'filepath')


# 10_SENSITIVITY ANALYSIS 1: Cohort matched on age_data_start, year_data_start, prac_region, IMD and gender -------

# 10a_Prepare data for interrupted time series ----

## Remove studyweek1 (as it's not a full week)
covid_weekly_conscounts_aggregate_England <- England_2015_2020_aggregate_weekly_conscounts_ds %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
covid_weekly_conscounts_aggregate_England_ap <- covid_weekly_conscounts_aggregate_England
covid_weekly_conscounts_aggregate_England_ap$lockdown1[covid_weekly_conscounts_aggregate_England_ap$date %in% adjustment_period] <- NA

# 10b_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_phone')

# 10c_Segmented regression model: All consultations -----

## Fit model 
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))

round(ci.lin(full_model_nb_all,Exp=T),3)
extract_glm_results_allages(full_model_nb_all, 'all_incladjustmentperiod_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_nolaggedresiduals_ds, "filepath" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

## Check autocorrelation
pacf(res)
acf(res)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box')

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_all <- full_join(covid_weekly_conscounts_aggregate_England_ap_all, lagres_timing, by = c('migrant_status','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + studyweek+ as.factor(studymonth)+lockdown1+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr, 'all_incladjustmentperiod_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_laggedresiduals_ds, "filepath" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2)
acf(res2)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

##  Make predictions
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_all
prediction_data$lagres <- 0
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_England_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_England_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_ds <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
row.names(allcons_RRs_ds) <- NULL
write.csv(allcons_RRs_ds, "filepath")

# Format for tables rmarkdown
allcons_RRs_ds <- allcons_RRs_ds %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(allcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, allcons) 

# 10d_Segmented regression model: Face-to-face consultations -----

# Fit model and calculate lagged residuals
full_model_nb_F2F <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F,Exp=T),3)

extract_glm_results_allages(full_model_nb_F2F, 'F2F_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_nolaggedresiduals_ds, "filepath" )

lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_F2F <- full_join(covid_weekly_conscounts_aggregate_England_ap_F2F, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_F2F_lr, 'F2F_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_laggedresiduals_ds, "filepath" )


# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

# Set up prediction dataset
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_F2F
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_F2F <- nb_pred_F2F$fit
stbp_nb_F2F <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_F2F <- covid_weekly_conscounts_aggregate_England_ap_F2F %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb_F2F, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn_F2F <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_F2F <- bind_cols(stbp = stbp_nb_F2F, stbp_noLdn= stbp_nb_noLdn_F2F, pred = pred_values_F2F, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_F2F <- bind_cols(df_se_nb_F2F, covid_weekly_conscounts_aggregate_England_ap_F2F) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_F2F <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates_F2F %>%
  mutate(var = rownames(parameter_estimates_F2F)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status_F2F <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_F2F <- post_lockdown_migrant_status_F2F %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_ds <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status_F2F)
row.names(F2Fcons_RRs_ds) <- NULL
write.csv(F2Fcons_RRs_ds, "filepath" )

# Format for tables rmarkdown
F2Fcons_RRs_ds <- F2Fcons_RRs_ds %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(F2Fcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, F2Fcons) 

# 10e_Segmented regression model: Phone consultations -----

# Fit model and calculate lagged residuals
full_model_nb_phone <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)

extract_glm_results_allages(full_model_nb_phone, 'phone_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_nolaggedresiduals_ds, "filepath" )


lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_phone <- full_join(covid_weekly_conscounts_aggregate_England_ap_phone, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_phone_lr, 'phone_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_laggedresiduals_ds, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Create prediction data set 
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_phone
prediction_data$lagres <- 0 # set to 0 based on discussion with Ali and Amy
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_phone <- nb_pred_phone$fit
stbp_nb_phone <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_phone <- covid_weekly_conscounts_aggregate_England_ap_phone %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_phone_nolockdown <- predict(final_model_nb_phone_lr, newdata = datanew2_nb_phone, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_phone_noLockdown <-nb_pred_phone_nolockdown$fit
stbp_nb_noLdn_phone <- nb_pred_phone_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_phone <- bind_cols(stbp = stbp_nb_phone, stbp_noLdn= stbp_nb_noLdn_phone, pred = pred_values_phone, pred_noLdn = predicted_vals_nb_phone_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_phone <- bind_cols(df_se_nb_phone, covid_weekly_conscounts_aggregate_England_ap_phone) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_phone <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates_phone %>%
  mutate(var = rownames(parameter_estimates_phone)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

## Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status_phone <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_phone <- post_lockdown_migrant_status_phone %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

## Combine RRs and save
phonecons_RRs_ds <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status_phone)
row.names(phonecons_RRs_ds) <- NULL
write.csv(phonecons_RRs_ds, "filepath" )

# Format for tables rmarkdown
phonecons_RRs_ds <- phonecons_RRs_ds %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(phonecons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, phonecons) 

# 10f_Combine plots to make a master plot -----

outcome_plot_nb$constype <- 'All consultations'
outcome_plot_nb_F2F$constype <- 'Face-to-face consultations'
outcome_plot_nb_phone$constype <- 'Telephone consultations'


all_outcome_plot_nb <- bind_rows(outcome_plot_nb, outcome_plot_nb_F2F, outcome_plot_nb_phone)

combined_England_plot_ds <- ggplot(filter(all_outcome_plot_nb, date > as.Date('2019-06-30') & date < as.Date('2020-11-29')), 
                                   aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, colour = migrant_status)) +
  geom_line()+
  geom_line(aes(y = incidence_rate), linetype = 2) +
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  labs(y = 'Rate (per person-year)', x = NULL, title = NULL) + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = NULL, date_breaks = '3 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_wrap(~constype, ncol = 1, scales = 'free_y')
combined_England_plot_ds

ggsave('filepath')

# 10g_Combine RRs for rmarkdown and save as Rdata file ----
RRs_England_SA <- bind_cols(allcons_RRs_ds, F2Fcons_RRs_ds, phonecons_RRs_ds) %>%
  dplyr::select(`var...1`, allcons, F2Fcons, phonecons) %>%
  rename('Variable' = `var...1`)
save(RRs_England_SA, file = 'filepath')

# 11_SENSITIVITY ANALYSIS 2: Cohort matched on age at study start and gender -------

# 11a_Prepare data for interrupted time series ----

## Remove studyweek1 (as it's not a full week)
covid_weekly_conscounts_aggregate_England <- England_2015_2020_aggregate_weekly_conscounts_SA2 %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
covid_weekly_conscounts_aggregate_England_ap <- covid_weekly_conscounts_aggregate_England
covid_weekly_conscounts_aggregate_England_ap$lockdown1[covid_weekly_conscounts_aggregate_England_ap$date %in% adjustment_period] <- NA

# 11b_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_phone')

# 11c_Segmented regression model: All consultations -----

## Fit model 
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))

round(ci.lin(full_model_nb_all,Exp=T),3)
extract_glm_results_allages(full_model_nb_all, 'all_incladjustmentperiod_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_nolaggedresiduals_ds, "filepath" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

## Check autocorrelation
pacf(res)
acf(res)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_all <- full_join(covid_weekly_conscounts_aggregate_England_ap_all, lagres_timing, by = c('migrant_status','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + studyweek+ as.factor(studymonth)+lockdown1+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr, 'all_incladjustmentperiod_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_laggedresiduals_ds, "filepath" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2)
acf(res2)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

##  Make predictions
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_all
prediction_data$lagres <- 0
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_England_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_England_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_SA2 <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
row.names(allcons_RRs_SA2) <- NULL
write.csv(allcons_RRs_SA2, "filepath")

# Format for tables rmarkdown
allcons_RRs_SA2 <- allcons_RRs_SA2 %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(allcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, allcons) 

# 11d_Segmented regression model: Face-to-face consultations -----

# Fit model and calculate lagged residuals
full_model_nb_F2F <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F,Exp=T),3)

extract_glm_results_allages(full_model_nb_F2F, 'F2F_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_nolaggedresiduals_ds, "filepath" )

lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_F2F <- full_join(covid_weekly_conscounts_aggregate_England_ap_F2F, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_F2F_lr, 'F2F_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_laggedresiduals_ds, "filepath" )


# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

# Set up prediction dataset
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_F2F
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_F2F <- nb_pred_F2F$fit
stbp_nb_F2F <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_F2F <- covid_weekly_conscounts_aggregate_England_ap_F2F %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb_F2F, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn_F2F <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_F2F <- bind_cols(stbp = stbp_nb_F2F, stbp_noLdn= stbp_nb_noLdn_F2F, pred = pred_values_F2F, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_F2F <- bind_cols(df_se_nb_F2F, covid_weekly_conscounts_aggregate_England_ap_F2F) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_F2F <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates_F2F %>%
  mutate(var = rownames(parameter_estimates_F2F)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status_F2F <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_F2F <- post_lockdown_migrant_status_F2F %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_SA2 <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status_F2F)
row.names(F2Fcons_RRs_SA2) <- NULL
write.csv(F2Fcons_RRs_SA2, "filepath" )

# Format for tables rmarkdown
F2Fcons_RRs_SA2 <- F2Fcons_RRs_SA2 %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(F2Fcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, F2Fcons) 

# 11e_Segmented regression model: Phone consultations -----

# Fit model and calculate lagged residuals
full_model_nb_phone <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)

extract_glm_results_allages(full_model_nb_phone, 'phone_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_nolaggedresiduals_ds, "filepath" )


lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_phone <- full_join(covid_weekly_conscounts_aggregate_England_ap_phone, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_phone_lr, 'phone_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_laggedresiduals_ds, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

# Create prediction data set 
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_phone
prediction_data$lagres <- 0 # set to 0 based on discussion with Ali and Amy
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_phone <- nb_pred_phone$fit
stbp_nb_phone <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_phone <- covid_weekly_conscounts_aggregate_England_ap_phone %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_phone_nolockdown <- predict(final_model_nb_phone_lr, newdata = datanew2_nb_phone, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_phone_noLockdown <-nb_pred_phone_nolockdown$fit
stbp_nb_noLdn_phone <- nb_pred_phone_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_phone <- bind_cols(stbp = stbp_nb_phone, stbp_noLdn= stbp_nb_noLdn_phone, pred = pred_values_phone, pred_noLdn = predicted_vals_nb_phone_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_phone <- bind_cols(df_se_nb_phone, covid_weekly_conscounts_aggregate_England_ap_phone) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_phone <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates_phone %>%
  mutate(var = rownames(parameter_estimates_phone)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

## Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status_phone <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_phone <- post_lockdown_migrant_status_phone %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

## Combine RRs and save
phonecons_RRs_SA2 <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status_phone)
row.names(phonecons_RRs_SA2) <- NULL
write.csv(phonecons_RRs_ds, "filepath" )

# Format for tables rmarkdown
phonecons_RRs_SA2 <- phonecons_RRs_SA2 %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(phonecons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, phonecons) 

# 11f_Combine plots to make a master plot -----

outcome_plot_nb$constype <- 'All consultations'
outcome_plot_nb_F2F$constype <- 'Face-to-face consultations'
outcome_plot_nb_phone$constype <- 'Telephone consultations'


all_outcome_plot_nb <- bind_rows(outcome_plot_nb, outcome_plot_nb_F2F, outcome_plot_nb_phone)

combined_England_plot_SA2 <- ggplot(filter(all_outcome_plot_nb, date > as.Date('2019-06-30') & date < as.Date('2020-11-29')), 
                                   aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, colour = migrant_status)) +
  geom_line()+
  geom_line(aes(y = incidence_rate), linetype = 2) +
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  labs(y = 'Rate (per person-year)', x = NULL, title = NULL) + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = NULL, date_breaks = '3 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_wrap(~constype, ncol = 1, scales = 'free_y')
combined_England_plot_SA2

ggsave('filepath')

# 11g_Combine RRs for rmarkdown and save as Rdata file ----
RRs_England_SA2 <- bind_cols(allcons_RRs_SA2, F2Fcons_RRs_SA2, phonecons_RRs_SA2) %>%
  dplyr::select(`var...1`, allcons, F2Fcons, phonecons) %>%
  rename('Variable' = `var...1`)
save(RRs_England_SA2, file = 'filepath')

# 11_SENSITIVITY ANALYSIS 3: Cohort matched on age at study start, prac_region and gender -------

# 11a_Prepare data for interrupted time series ----

## Remove studyweek1 (as it's not a full week)
covid_weekly_conscounts_aggregate_England <- England_2015_2020_aggregate_weekly_conscounts_SA3 %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
covid_weekly_conscounts_aggregate_England_ap <- covid_weekly_conscounts_aggregate_England
covid_weekly_conscounts_aggregate_England_ap$lockdown1[covid_weekly_conscounts_aggregate_England_ap$date %in% adjustment_period] <- NA

# 11b_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_phone')

# 11c_Segmented regression model: All consultations -----

## Fit model 
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))

round(ci.lin(full_model_nb_all,Exp=T),3)
extract_glm_results_allages(full_model_nb_all, 'all_incladjustmentperiod_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_nolaggedresiduals_ds, "filepath" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

## Check autocorrelation
pacf(res)
acf(res)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_all <- full_join(covid_weekly_conscounts_aggregate_England_ap_all, lagres_timing, by = c('migrant_status','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + studyweek+ as.factor(studymonth)+lockdown1+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr, 'all_incladjustmentperiod_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_all, 'lockdown1')
write.csv(all_incladjustmentperiod_laggedresiduals_ds, "filename" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'deviance')
pacf(res2)
acf(res2)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

##  Make predictions
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_all
prediction_data$lagres <- 0
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- covid_weekly_conscounts_aggregate_England_ap_all %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_all_nolockdown <- predict(final_model_nb_all_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_all_noLockdown <-nb_pred_all_nolockdown$fit
stbp_nb_noLdn <- nb_pred_all_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_all_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, covid_weekly_conscounts_aggregate_England_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_all_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
allcons_RRs_SA3 <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_migrant_status)
row.names(allcons_RRs_SA3) <- NULL
write.csv(allcons_RRs_SA3, "filepath")

# Format for tables rmarkdown
allcons_RRs_SA3 <- allcons_RRs_SA3 %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(allcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, allcons) 

# 11d_Segmented regression model: Face-to-face consultations -----

# Fit model and calculate lagged residuals
full_model_nb_F2F <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F,Exp=T),3)

extract_glm_results_allages(full_model_nb_F2F, 'F2F_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_nolaggedresiduals_ds, "filepath" )

lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_F2F <- full_join(covid_weekly_conscounts_aggregate_England_ap_F2F, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(covid_weekly_conscounts_aggregate_England_ap_F2F, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_F2F_lr, 'F2F_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_F2F, 'lockdown1')
write.csv(F2F_laggedresiduals_ds, "filepath" )


# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Set up prediction dataset
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_F2F
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_F2F <- nb_pred_F2F$fit
stbp_nb_F2F <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_F2F <- covid_weekly_conscounts_aggregate_England_ap_F2F %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb_F2F, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn_F2F <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_F2F <- bind_cols(stbp = stbp_nb_F2F, stbp_noLdn= stbp_nb_noLdn_F2F, pred = pred_values_F2F, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_F2F <- bind_cols(df_se_nb_F2F, covid_weekly_conscounts_aggregate_England_ap_F2F) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_F2F <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates_F2F %>%
  mutate(var = rownames(parameter_estimates_F2F)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status_F2F <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_F2F <- post_lockdown_migrant_status_F2F %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_SA3 <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status_F2F)
row.names(F2Fcons_RRs_SA3) <- NULL
write.csv(F2Fcons_RRs_SA3, "filepath" )

# Format for tables rmarkdown
F2Fcons_RRs_SA3 <- F2Fcons_RRs_SA3 %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(F2Fcons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, F2Fcons) 

# 11e_Segmented regression model: Phone consultations -----

# Fit model and calculate lagged residuals
full_model_nb_phone <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)

extract_glm_results_allages(full_model_nb_phone, 'phone_nolaggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_nolaggedresiduals_ds, "filepath" )


lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res)
acf(res)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') 

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
covid_weekly_conscounts_aggregate_England_ap_phone <- full_join(covid_weekly_conscounts_aggregate_England_ap_phone, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(covid_weekly_conscounts_aggregate_England_ap_phone, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_phone_lr, 'phone_laggedresiduals_ds', covid_weekly_conscounts_aggregate_England_ap_phone, 'lockdown1')
write.csv(phone_laggedresiduals_ds, "filepath" )

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2)
acf(res2)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box')  

# Create prediction data set 
prediction_data <- covid_weekly_conscounts_aggregate_England_ap_phone
prediction_data$lagres <- 0 # set to 0 based on discussion with Ali and Amy
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values_phone <- nb_pred_phone$fit
stbp_nb_phone <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb_phone <- covid_weekly_conscounts_aggregate_England_ap_phone %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_phone_nolockdown <- predict(final_model_nb_phone_lr, newdata = datanew2_nb_phone, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_phone_noLockdown <-nb_pred_phone_nolockdown$fit
stbp_nb_noLdn_phone <- nb_pred_phone_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb_phone <- bind_cols(stbp = stbp_nb_phone, stbp_noLdn= stbp_nb_noLdn_phone, pred = pred_values_phone, pred_noLdn = predicted_vals_nb_phone_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

# Combine data set and predictions
outcome_plot_nb_phone <- bind_cols(df_se_nb_phone, covid_weekly_conscounts_aggregate_England_ap_phone) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates_phone <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates_phone %>%
  mutate(var = rownames(parameter_estimates_phone)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

## Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status_phone <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status_phone <- post_lockdown_migrant_status_phone %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

## Combine RRs and save
phonecons_RRs_SA3 <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status_phone)
row.names(phonecons_RRs_SA3) <- NULL
write.csv(phonecons_RRs_SA3, "filepath" )

# Format for tables rmarkdown
phonecons_RRs_SA3 <- phonecons_RRs_SA3 %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(phonecons = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, phonecons) 

# 11f_Combine plots to make a master plot -----

outcome_plot_nb$constype <- 'All consultations'
outcome_plot_nb_F2F$constype <- 'Face-to-face consultations'
outcome_plot_nb_phone$constype <- 'Telephone consultations'


all_outcome_plot_nb <- bind_rows(outcome_plot_nb, outcome_plot_nb_F2F, outcome_plot_nb_phone)

combined_England_plot_SA3 <- ggplot(filter(all_outcome_plot_nb, date > as.Date('2019-06-30') & date < as.Date('2020-11-29')), 
                                    aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, colour = migrant_status)) +
  geom_line()+
  geom_line(aes(y = incidence_rate), linetype = 2) +
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  labs(y = 'Rate (per person-year)', x = NULL, title = NULL) + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = NULL, date_breaks = '3 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_wrap(~constype, ncol = 1, scales = 'free_y')
combined_England_plot_SA3

ggsave('filepath')

# 11g_Combine RRs for rmarkdown and save as Rdata file ----
RRs_England_SA3 <- bind_cols(allcons_RRs_SA3, F2Fcons_RRs_SA3, phonecons_RRs_SA3) %>%
  dplyr::select(`var...1`, allcons, F2Fcons, phonecons) %>%
  rename('Variable' = `var...1`)
save(RRs_England_SA3, file = 'filepath')




