## 0_Description ---------------------------------------------------------------------------

# Migrants' primary care utilisation before and during the COVID-19 pandemic in England: An interrupted time series
# 2015-2020 interrupted time series analysis - code for analyses not included in paper and model choice
# Date started: 25/03/2021
# Author(s): Yamina Boukari / Claire Zhang 

## 0_Load packages -------------------------------------------------------------------------

## 0_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc, epitools, data.table, epiR, MASS, psych, forestplot, foreign,tsModel,lmtest,
               Epi, multcomp,splines, vcd, here, stringr, patchwork, gmodels, gtsummary, scales, lmtest, cowplot, grid, gridExtra, sandwich)


## 0_Set working directory ------------------------------------------------------------------

setwd("filepath")

## 0_Load datasets and functions ------------------------------------------------------

# Datasets
load(file = "filepath/England_2015_2020_aggregate_weekly_conscounts.Rdata")
load(file = "filepath/England_2015_2020_aggregate_weekly_conscounts_agesubcohort.Rdata")
load(file = "filepath/London_2015_2020_aggregate_weekly_conscounts.Rdata")
load(file = "filepath/England_2015_2020_aggregate_weekly_conscounts_ethnicity.Rdata")
load(file = "filepath/England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity.Rdata")
load(file = "filepath/England_2015_2020_aggregate_weekly_conscounts_ds.Rdata")
load(file = "filepath/England_2015_2020_aggregate_weekly_conscounts_migcertainty.Rdata")
load(file = "filepath/England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort.Rdata")

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

# 1_Migrant_status and 6-category ethnicity ------

# 1a_Prepare data for interrupted time series ----

## Remove studyweek1 
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap <- England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity
England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap$lockdown1[England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap$date %in% adjustment_period] <- NA

# 1b_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR_F2F(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, 'England_weekly_conscounts_migstatus_ethnicity_F2F')
function_to_calculate_IR_phone(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity, 'England_weekly_conscounts_migstatus_ethnicity_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR_F2F(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap, 'England_weekly_conscounts_migstatus_ethnicity_ap_F2F')
function_to_calculate_IR_phone(England_2015_2020_aggregate_weekly_conscounts_migstatus_ethnicity_ap, 'England_weekly_conscounts_migstatus_ethnicity_ap_phone')

# 1c_Segmented regression model: Face-to-face consultations -----

## Fit model 
full_model_nb_F2F <- glm.nb(facetoface ~ offset(log(person_years))+ studyweek + as.factor(studymonth) + migrant_status*ethnicat6*lockdown1,
                            data = filter(England_weekly_conscounts_migstatus_ethnicity_ap_F2F, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F,Exp=T),3)
extract_glm_results_allages(full_model_nb_F2F, 'migstatus_ethnicity_F2F_incladjustmentperiod_nolaggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_F2F, 'lockdown1')
write.csv(migstatus_ethnicity_F2F_incladjustmentperiod_nolaggedresiduals, "filepath/migstatus_ethnicity_F2F_incladjustmentperiod_nolaggedresiduals.csv" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

## Check autocorrelation
pacf(res, lag = 368)
acf(res,lag = 368)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05, so evidence of autocorrelation 

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(England_weekly_conscounts_migstatus_ethnicity_ap_F2F, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status, ethnicat6))
England_weekly_conscounts_migstatus_ethnicity_ap_F2F <- full_join(England_weekly_conscounts_migstatus_ethnicity_ap_F2F, lagres_timing, by = c('migrant_status', 'ethnicat6','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_F2F_lr <- glm.nb(facetoface ~ offset(log(person_years)) + studyweek+ as.factor(studymonth)+ migrant_status*ethnicat6*lockdown1 + 
                                  lagres,
                                data =  filter(England_weekly_conscounts_migstatus_ethnicity_ap_F2F, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_F2F_lr, 'migstatus_ethnicity_F2F_incladjustmentperiod_laggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_F2F, 'lockdown1')
write.csv(migstatus_ethnicity_F2F_incladjustmentperiod_laggedresiduals, "filepath/migstatus_ethnicity_F2F_incladjustmentperiod_laggedresiduals.csv" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'response')
pacf(res2, lag = 312)
acf(res2,lag = 312)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p>0.05

##  Make predictions
prediction_data <- England_weekly_conscounts_migstatus_ethnicity_ap_F2F
prediction_data$lagres <- 0
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- England_weekly_conscounts_migstatus_ethnicity_ap_F2F %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, England_weekly_conscounts_migstatus_ethnicity_ap_F2F) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)


# 1d_Plot the results and save ----

## Create dataframe for the vertical line lockdown labels (i.e. the adjustment-to-restrictions period in the COVID collateral paper - 8th March to the 28th March)
lockdowns <- data.frame(date = as_date(c('2020-03-08', '2020-03-29')))

var_width = 20
outcome_plot_nb <- mutate(outcome_plot_nb, tidy_ethnicat6_name = str_wrap(ethnicat6, width = var_width))
outcome_plot_nb$tidy_ethnicat6_name <- factor(outcome_plot_nb$tidy_ethnicat6_name,levels = c("White British", "White non-British", "Mixed/Multiple\nethnic groups", 
                                                                                             "Asian/Asian British", "Black/African/\nCaribbean/Black\nBritish", "Other ethnic group", 'Unknown'))

F2F_cons_migstatus_ethnicity_nb <- outcome_plot_nb %>%
  dplyr::select(c(migrant_status, date, final_pred, final_low, final_upp, incidence_rate, tidy_ethnicat6_name)) %>%
  tidyr::pivot_longer(cols = c('final_pred', 'incidence_rate')) %>%
  filter(date > as.Date('2019-06-30') & date < as.Date('2020-11-29')) %>%
  ggplot(mapping = aes(x = date, y = value, ymin = final_low, ymax = final_upp, fill = migrant_status, colour = migrant_status, linetype = name)) +
  geom_line(aes(colour = migrant_status, linetype = name))+
  geom_ribbon(alpha = 0.3, colour = 0)+ 
  labs(y = 'Rate (per person-year)', x = NULL, title = NULL) + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = NULL, date_breaks = '6 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b')+
  theme_bw() +
  scale_linetype_manual(values = c(1, 2), labels = c('Predicted', 'Observed')) +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  facet_wrap(~tidy_ethnicat6_name, dir = 'h', nrow = 2)

F2F_cons_migstatus_ethnicity_nb
ggsave('filepath/F2F_cons_migstatus_ethnicity_nb.png')

# 1e_Get rate ratios -----

## Get RR migrants of different ethnicities versus the reference group (white British non-migrants)
names(coef(final_model_nb_F2F_lr))

### White British migrants vs WBNM
post_lockdown_WBvsWBNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1), conf=.95))
post_lockdown_WBvsWBNM <- post_lockdown_WBvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1')
### White non-British migrants vs WBNM
post_lockdown_WNBvsWBNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'ethnicat6White non-British:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6White non-British:lockdown1' = 1), conf=.95))
post_lockdown_WNBvsWBNM <- post_lockdown_WNBvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6White non-British:lockdown1 + migrant_statusMigrant:ethnicat6White non-British:lockdown1')
### Mixed migrants vs WBNM
post_lockdown_MixedvsWBNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1), conf=.95))
post_lockdown_MixedvsWBNM <- post_lockdown_MixedvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Mixed:lockdown1 + migrant_statusMigrant:ethnicat6Mixed:lockdown1')
### Asian migrants vs WBNM
post_lockdown_AsianvsWBNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'ethnicat6Asian/Asian British:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Asian/Asian British:lockdown1' = 1), conf=.95))
post_lockdown_AsianvsWBNM <- post_lockdown_AsianvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Asian:lockdown1 + migrant_statusMigrant:ethnicat6Asian:lockdown1')
### Black migrants vs WBNM
post_lockdown_BlackvsWBNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1), conf=.95))
post_lockdown_BlackvsWBNM <- post_lockdown_BlackvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Black:lockdown1 + migrant_statusMigrant:ethnicat6Black:lockdown1')
### Other migrants vs WBNM
post_lockdown_OthervsWBNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'ethnicat6Other ethnic group:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Other ethnic group:lockdown1' = 1), conf=.95))
post_lockdown_OthervsWBNM <- post_lockdown_OthervsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Other:lockdown1 + migrant_statusMigrant:ethnicat6Other:lockdown1')
### Unknown migrants vs WBNM
post_lockdown_UnknownvsWBNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'ethnicat6Unknown:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicat6Unknown:lockdown1' = 1), conf=.95))
post_lockdown_UnknownvsWBNM <- post_lockdown_UnknownvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Unknown:lockdown1 + migrant_statusMigrant:ethnicat6Unknown:lockdown1')

## Combine RRs and save
F2Fcons_migstatus_ethnicity_RRs_vsWBNMs <- bind_rows(post_lockdown_WBvsWBNM, post_lockdown_WNBvsWBNM,post_lockdown_MixedvsWBNM,
                                                     post_lockdown_AsianvsWBNM, post_lockdown_BlackvsWBNM, post_lockdown_OthervsWBNM,
                                                     post_lockdown_UnknownvsWBNM)
row.names(F2Fcons_migstatus_ethnicity_RRs_vsWBNMs) <- NULL
write.csv(F2Fcons_migstatus_ethnicity_RRs_vsWBNMs, "filepath/F2Fcons_migstatus_ethnicity_RRs_vsWBNMs.csv" )

## Prepare data for forest plot

names <- c('White British', 'White non-British', 'Mixed/Multiple ethnic groups', 'Asian/Asian British', 'Black/African/Caribbean/Black British',
           'Other ethnic group', 'Unknown')
F2Fcons_migstatus_ethnicity_RRs_vsWBNMs$names <- names
F2Fcons_migstatus_ethnicity_RRs_vsWBNMs$lower <-  round(F2Fcons_migstatus_ethnicity_RRs_vsWBNMs$lower, digits = 2)
F2Fcons_migstatus_ethnicity_RRs_vsWBNMs$upper <-  round(F2Fcons_migstatus_ethnicity_RRs_vsWBNMs$upper, digits = 2)
F2Fcons_migstatus_ethnicity_RRs_vsWBNMs$estimate <-  round(F2Fcons_migstatus_ethnicity_RRs_vsWBNMs$estimate, digits = 2)
F2Fcons_migstatus_ethnicity_RRs_vsWBNMs <- F2Fcons_migstatus_ethnicity_RRs_vsWBNMs %>%
  mutate(ci = paste0(lower, '-', upper)) %>%
  mutate(outcome = c('Migrants vs. White British non-migrants', ' ', ' ', ' ', ' ', ' ', ' '),
         p = ' ') %>%
  add_row(.before = 1)

## Get RR for migrants vs non-migrants of different ethnicities 
names(coef(final_model_nb_F2F_lr))

### White British migrants vs non-migrants
post_lockdown_WBvsNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1), conf=.95))
post_lockdown_WBvsNM <- post_lockdown_WBvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1')
### White non-British migrants vs non-migrants
post_lockdown_WNBvsNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                'migrant_statusMigrant:ethnicat6White non-British:lockdown1' = 1), conf=.95))
post_lockdown_WNBvsNM <- post_lockdown_WNBvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6White non-British:lockdown1')
### Mixed migrants vs non-migrants
post_lockdown_MixedvsNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1), conf=.95))
post_lockdown_MixedvsNM <- post_lockdown_MixedvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Mixed:lockdown1')
### Asian migrants vs non-migrants
post_lockdown_AsianvsNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6Asian/Asian British:lockdown1' = 1), conf=.95))
post_lockdown_AsianvsNM <- post_lockdown_AsianvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Asian:lockdown1')
### Black migrants vs non-migrants
post_lockdown_BlackvsNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1), conf=.95))
post_lockdown_BlackvsNM <- post_lockdown_BlackvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Black:lockdown1')
### Other migrants vs non-migrants
post_lockdown_OthervsNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6Other ethnic group:lockdown1' = 1), conf=.95))
post_lockdown_OthervsNM <- post_lockdown_OthervsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Other:lockdown1')
### Unknown migrants vs non-migrants
post_lockdown_UnknownvsNM <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Unknown:lockdown1' = 1), conf=.95))
post_lockdown_UnknownvsNM <- post_lockdown_UnknownvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Unknown:lockdown1')

# Combine RRs and save
F2Fcons_migstatus_ethnicity_RRs_vsNMs <- bind_rows(post_lockdown_WBvsNM, post_lockdown_WNBvsNM,post_lockdown_MixedvsNM,
                                                   post_lockdown_AsianvsNM, post_lockdown_BlackvsNM, post_lockdown_OthervsNM,
                                                   post_lockdown_UnknownvsNM)
row.names(F2Fcons_migstatus_ethnicity_RRs_vsNMs) <- NULL
write.csv(F2Fcons_migstatus_ethnicity_RRs_vsNMs, "filepath/F2Fcons_migstatus_ethnicity_RRs_vsNMs.csv")


# 1f_Segmented regression model: Phone consultations -----

## Fit model 
full_model_nb_phone <- glm.nb(phone ~ offset(log(person_years))+ studyweek + as.factor(studymonth) + migrant_status*ethnicat6*lockdown1,
                              data = filter(England_weekly_conscounts_migstatus_ethnicity_ap_phone, !is.na(lockdown1)))
round(ci.lin(full_model_nb_phone,Exp=T),3)
extract_glm_results_allages(full_model_nb_phone, 'migstatus_ethnicity_phone_incladjustmentperiod_nolaggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_phone, 'lockdown1')
write.csv(migstatus_ethnicity_phone_incladjustmentperiod_nolaggedresiduals, "filepath/migstatus_ethnicity_phone_incladjustmentperiod_nolaggedresiduals.csv" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

## Check autocorrelation
pacf(res, lag = 368)
acf(res,lag = 368)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05, so evidence of autocorrelation 

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(England_weekly_conscounts_migstatus_ethnicity_ap_phone, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status, ethnicat6))
England_weekly_conscounts_migstatus_ethnicity_ap_phone <- full_join(England_weekly_conscounts_migstatus_ethnicity_ap_phone, lagres_timing, by = c('migrant_status', 'ethnicat6','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_phone_lr <- glm.nb(phone ~ offset(log(person_years)) + studyweek+ as.factor(studymonth)+ migrant_status*ethnicat6*lockdown1 + 
                                    lagres,
                                  data =  filter(England_weekly_conscounts_migstatus_ethnicity_ap_phone, !is.na(lockdown1)))
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_phone_lr, 'migstatus_ethnicity_phone_incladjustmentperiod_laggedresiduals', England_weekly_conscounts_migstatus_ethnicity_ap_phone, 'lockdown1')
write.csv(migstatus_ethnicity_phone_incladjustmentperiod_laggedresiduals, "filepath/migstatus_ethnicity_phone_incladjustmentperiod_laggedresiduals.csv" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'response')
pacf(res2, lag = 312)
acf(res2,lag = 312)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05 - still evidence of autocorrelation 

##  Make predictions
prediction_data <- England_weekly_conscounts_migstatus_ethnicity_ap_phone
prediction_data$lagres <- 0
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- England_weekly_conscounts_migstatus_ethnicity_ap_phone %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_phone_nolockdown <- predict(final_model_nb_phone_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_phone_noLockdown <-nb_pred_phone_nolockdown$fit
stbp_nb_noLdn <- nb_pred_phone_nolockdown$se.fit

## Combine predictions 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_phone_noLockdown) %>%
  mutate(
    #CIs
    upp = pred + (1.96*stbp),
    low = pred - (1.96*stbp),
    upp_noLdn = pred_noLdn + (1.96*stbp_noLdn),
    low_noLdn = pred_noLdn - (1.96*stbp_noLdn)
  )

## Combine data set and predictions
outcome_plot_nb <- bind_cols(df_se_nb, England_weekly_conscounts_migstatus_ethnicity_ap_phone) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)


# 1g_Plot the results and save -----

phone_cons_migstatus_ethnicity_nb <- ggplot(filter(outcome_plot_nb, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), 
                                            aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = migrant_status, color = migrant_status)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  geom_line(aes(y = incidence_rate), linetype = 2) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Consultations') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '2 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank())+
  ylim(0,6) + 
  facet_wrap(~ethnicat6)

phone_cons_migstatus_ethnicity_nb
ggsave('filepath/phone_cons_migstatus_ethnicity_nb.png')

# 1h_Get rate ratios -----

## Get RR migrants of different ethnicities versus the reference group (white British non-migrants)
names(coef(final_model_nb_phone_lr))

### White British migrants vs WBNM
post_lockdown_WBvsWBNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1), conf=.95))
post_lockdown_WBvsWBNM <- post_lockdown_WBvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1')
### White non-British migrants vs WBNM
post_lockdown_WNBvsWBNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'ethnicat6White non-British:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6White non-British:lockdown1' = 1), conf=.95))
post_lockdown_WNBvsWBNM <- post_lockdown_WNBvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6White non-British:lockdown1 + migrant_statusMigrant:ethnicat6White non-British:lockdown1')
### Mixed migrants vs WBNM
post_lockdown_MixedvsWBNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1), conf=.95))
post_lockdown_MixedvsWBNM <- post_lockdown_MixedvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Mixed:lockdown1 + migrant_statusMigrant:ethnicat6Mixed:lockdown1')
### Asian migrants vs WBNM
post_lockdown_AsianvsWBNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'ethnicat6Asian/Asian British:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicat6Asian/Asian British:lockdown1' = 1), conf=.95))
post_lockdown_AsianvsWBNM <- post_lockdown_AsianvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Asian:lockdown1 + migrant_statusMigrant:ethnicat6Asian:lockdown1')
### Black migrants vs WBNM
post_lockdown_BlackvsWBNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1), conf=.95))
post_lockdown_BlackvsWBNM <- post_lockdown_BlackvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Black:lockdown1 + migrant_statusMigrant:ethnicat6Black:lockdown1')
### Other migrants vs WBNM
post_lockdown_OthervsWBNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'ethnicat6Other ethnic group:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicat6Other ethnic group:lockdown1' = 1), conf=.95))
post_lockdown_OthervsWBNM <- post_lockdown_OthervsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Other:lockdown1 + migrant_statusMigrant:ethnicat6Other:lockdown1')
### Unknown migrants vs WBNM
post_lockdown_UnknownvsWBNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                        'ethnicat6Unknown:lockdown1' = 1,
                                                                        'migrant_statusMigrant:ethnicat6Unknown:lockdown1' = 1), conf=.95))
post_lockdown_UnknownvsWBNM <- post_lockdown_UnknownvsWBNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + ethnicat6Unknown:lockdown1 + migrant_statusMigrant:ethnicat6Unknown:lockdown1')

## Combine RRs and save
phonecons_migstatus_ethnicity_RRs_vsWBNMs <- bind_rows(post_lockdown_WBvsWBNM, post_lockdown_WNBvsWBNM,post_lockdown_MixedvsWBNM,
                                                       post_lockdown_AsianvsWBNM, post_lockdown_BlackvsWBNM, post_lockdown_OthervsWBNM,
                                                       post_lockdown_UnknownvsWBNM)
row.names(phonecons_migstatus_ethnicity_RRs_vsWBNMs) <- NULL
write.csv(phonecons_migstatus_ethnicity_RRs_vsWBNMs, "filepath/phonecons_migstatus_ethnicity_RRs_vsWBNMs.csv" )

## Prepare data for forest plot

names <- c('White British', 'White non-British', 'Mixed/Multiple ethnic groups', 'Asian/Asian British', 'Black/African/Caribbean/Black British',
           'Other ethnic group', 'Unknown')
phonecons_migstatus_ethnicity_RRs_vsWBNMs$names <- names
phonecons_migstatus_ethnicity_RRs_vsWBNMs$lower <-  round(phonecons_migstatus_ethnicity_RRs_vsWBNMs$lower, digits = 2)
phonecons_migstatus_ethnicity_RRs_vsWBNMs$upper <-  round(phonecons_migstatus_ethnicity_RRs_vsWBNMs$upper, digits = 2)
phonecons_migstatus_ethnicity_RRs_vsWBNMs$estimate <-  round(phonecons_migstatus_ethnicity_RRs_vsWBNMs$estimate, digits = 2)
phonecons_migstatus_ethnicity_RRs_vsWBNMs <- phonecons_migstatus_ethnicity_RRs_vsWBNMs %>%
  mutate(ci = paste0(lower, '-', upper)) %>%
  mutate(outcome = c('Migrants vs. White British non-migrants', ' ', ' ', ' ', ' ', ' ', ' '),
         p = ' ') %>%
  add_row(.before = 1)

## Get RR for migrants vs non-migrants of different ethnicities 
names(coef(final_model_nb_phone_lr))

### White British migrants vs non-migrants
post_lockdown_WBvsNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1), conf=.95))
post_lockdown_WBvsNM <- post_lockdown_WBvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1')
### White non-British migrants vs non-migrants
post_lockdown_WNBvsNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                  'migrant_statusMigrant:ethnicat6White non-British:lockdown1' = 1), conf=.95))
post_lockdown_WNBvsNM <- post_lockdown_WNBvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6White non-British:lockdown1')
### Mixed migrants vs non-migrants
post_lockdown_MixedvsNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Mixed/Multiple ethnic groups:lockdown1' = 1), conf=.95))
post_lockdown_MixedvsNM <- post_lockdown_MixedvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Mixed:lockdown1')
### Asian migrants vs non-migrants
post_lockdown_AsianvsNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Asian/Asian British:lockdown1' = 1), conf=.95))
post_lockdown_AsianvsNM <- post_lockdown_AsianvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Asian:lockdown1')
### Black migrants vs non-migrants
post_lockdown_BlackvsNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Black/African/Caribbean/Black British:lockdown1' = 1), conf=.95))
post_lockdown_BlackvsNM <- post_lockdown_BlackvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Black:lockdown1')
### Other migrants vs non-migrants
post_lockdown_OthervsNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                    'migrant_statusMigrant:ethnicat6Other ethnic group:lockdown1' = 1), conf=.95))
post_lockdown_OthervsNM <- post_lockdown_OthervsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Other:lockdown1')
### Unknown migrants vs non-migrants
post_lockdown_UnknownvsNM <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant:lockdown1' = 1,
                                                                      'migrant_statusMigrant:ethnicat6Unknown:lockdown1' = 1), conf=.95))
post_lockdown_UnknownvsNM <- post_lockdown_UnknownvsNM %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('lower' = Lower.CI,
         'upper' = Upper.CI,
         'estimate' = Estimate) %>%
  mutate(var = 'migrant_statusMigrant:lockdown1 + migrant_statusMigrant:ethnicat6Unknown:lockdown1')

# Combine RRs and save
phonecons_migstatus_ethnicity_RRs_vsNMs <- bind_rows(post_lockdown_WBvsNM, post_lockdown_WNBvsNM,post_lockdown_MixedvsNM,
                                                     post_lockdown_AsianvsNM, post_lockdown_BlackvsNM, post_lockdown_OthervsNM,
                                                     post_lockdown_UnknownvsNM)
row.names(phonecons_migstatus_ethnicity_RRs_vsNMs) <- NULL
write.csv(phonecons_migstatus_ethnicity_RRs_vsNMs, "filepath/phonecons_migstatus_ethnicity_RRs_vsNMs.csv")

# Prepare data for forest plot

names <- c('White British', 'White non-British', 'Mixed/Multiple ethnic groups', 'Asian/Asian British', 'Black/African/Caribbean/Black British',
           'Other ethnic group', 'Unknown')
phonecons_migstatus_ethnicity_RRs_vsNMs$names <- names
phonecons_migstatus_ethnicity_RRs_vsNMs$lower <-  round(phonecons_migstatus_ethnicity_RRs_vsNMs$lower, digits = 2)
phonecons_migstatus_ethnicity_RRs_vsNMs$upper <-  round(phonecons_migstatus_ethnicity_RRs_vsNMs$upper, digits = 2)
phonecons_migstatus_ethnicity_RRs_vsNMs$estimate <-  round(phonecons_migstatus_ethnicity_RRs_vsNMs$estimate, digits = 2)
phonecons_migstatus_ethnicity_RRs_vsNMs <- phonecons_migstatus_ethnicity_RRs_vsNMs %>%
  mutate(ci = paste0(lower, '-', upper)) %>%
  mutate(outcome = c('Migrants vs. non-migrants of the same ethnicity', ' ', ' ', ' ', ' ', ' ', ' '),
         p = ' ') %>%
  add_row(.before = 1)

# 1i_Combine forest plot data ----

forest_plot_migstatus_ethnicity <- bind_rows(phonecons_migstatus_ethnicity_RRs_vsWBNMs, phonecons_migstatus_ethnicity_RRs_vsNMs) 

# Make forest plot using forest plot
tabletext <- cbind(c(" ", forest_plot_migstatus_ethnicity$outcome),
                   c(" ", forest_plot_migstatus_ethnicity$names),
                   c("RR", forest_plot_migstatus_ethnicity$estimate),
                   c("95% CI", forest_plot_migstatus_ethnicity$ci))

forest_plot_migstatus_ethnicity <- forest_plot_migstatus_ethnicity %>%
  add_row(.before = 1) # Needs to come after creating the tabletext object so that there are equal number of rows in the labels and means

dev.new()
png(file = "S:/CALIBER_19_062R/01_gold/01_all_patients/results/01_Consultations/forest_plots/forest_plot_migstatus_ethnicity_phone.png", width = 6000, height = 3000)
forest_plot <- forestplot(tabletext, mean = forest_plot_migstatus_ethnicity$estimate, lower=forest_plot_migstatus_ethnicity$lower, upper=forest_plot_migstatus_ethnicity$upper,
                          graph.pos = (ncol(tabletext)-1),
                          clip = c(-2,2),
                          zero = 1, 
                          col = fpColors(box="mediumblue", lines = "black", zero="black", hrz_lines = "black"),
                          xlab = "Rate Ratios",
                          #is.summary = c(rep(FALSE,3), TRUE, rep(FALSE, 26), TRUE, rep(FALSE, 12), TRUE, rep(FALSE, 15)),
                          txt_gp = fpTxtGp(ticks=gpar(cex=4), xlab=gpar(cex=4), cex = 4, summary = gpar(fontface = 'bold')),
                          xticks = c(0.5,1,1.5),
                          ci.vertices = TRUE,
                          ci.vertices.height = 0.2,
                          boxsize = 0.5,
                          title = "Effect of migration and ethnicity on phone consultation rates pre- and during the pandemic",
                          graphwidth = unit(700, 'mm'),
                          align = c('l', 'l', 'l'),
                          alim = c(0.5, 1.5))
dev.off()


# 2_SENSITIVITY ANALYSIS: age_subcohort ----

# 2a_Prepare data for interrupted time series -----

## Remove studyweek1 
covid_weekly_conscounts_aggregate_England_agesubcohort <- England_2015_2020_aggregate_weekly_conscounts_agesubcohort %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020
covid_weekly_conscounts_aggregate_England_agesubcohort_ap <- covid_weekly_conscounts_aggregate_England_agesubcohort
covid_weekly_conscounts_aggregate_England_agesubcohort_ap$lockdown1[covid_weekly_conscounts_aggregate_England_agesubcohort_ap$date %in% adjustment_period] <- NA

# 2b_Calculate IRs ----
# all data 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_agesubcohort, 'covid_weekly_conscounts_England_agesubcohort_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_agesubcohort, 'covid_weekly_conscounts_England_agesubcohort_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_agesubcohort, 'covid_weekly_conscounts_England_agesubcohort_phone')
# truncated data with lockdown1 set to NA during the adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_agesubcohort_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone')

# 2c_Segmented regression model: Face-to-face consultations (not included in paper) ----

# 2d_0-15 years -----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F, age_subcohort == '0-15 years')
full_model_nb_F2F <- glm.nb(facetoface ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
round(ci.lin(full_model_nb_F2F,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p < 0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
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


# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_0_15years <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status)
rownames(F2Fcons_RRs_0_15years) <- NULL
write.csv(F2Fcons_RRs_0_15years, "filepath/F2Fcons_RRs_0_15years.csv" )

# 2e_16-24 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F, age_subcohort == '16-24 years')
full_model_nb_F2F <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_F2F)
round(ci.lin(full_model_nb_F2F,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_F2F_lr)
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p > 0.05 - accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
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

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_16_24years <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status)
rownames(F2Fcons_RRs_16_24years) <- NULL
write.csv(F2Fcons_RRs_16_24years, "filepath/F2Fcons_RRs_16_24years.csv" )

# 2f_25-34 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F, age_subcohort == '25-34 years')
full_model_nb_F2F <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_F2F)
round(ci.lin(full_model_nb_F2F,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p < 0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_F2F_lr)
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p > 0.05 -- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
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

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_25_34years <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status)
rownames(F2Fcons_RRs_25_34years) <- NULL
write.csv(F2Fcons_RRs_25_34years, "filepath/F2Fcons_RRs_25_34years.csv" )

# 2g_35-49 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F, age_subcohort == '35-49 years')
full_model_nb_F2F <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_F2F)
round(ci.lin(full_model_nb_F2F,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_F2F_lr)
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
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

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_35_49years <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status)
rownames(F2Fcons_RRs_35_49years) <- NULL
write.csv(F2Fcons_RRs_35_49years, "filepath/F2Fcons_RRs_35_49years.csv" )

# 2h_50-64 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F, age_subcohort == '50-64 years')
full_model_nb_F2F <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_F2F)
round(ci.lin(full_model_nb_F2F,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_F2F_lr)
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
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

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_50_64years <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status)
rownames(F2Fcons_RRs_50_64years) <- NULL
write.csv(F2Fcons_RRs_50_64years, "filepath/F2Fcons_RRs_50_64years.csv" )

# 2i_65+ years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_F2F, age_subcohort == '>=65 years')
full_model_nb_F2F <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_F2F)
round(ci.lin(full_model_nb_F2F,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_F2F)) %>% as.numeric()
res <- residuals(full_model_nb_F2F, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_F2F, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p <0.005

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_F2F_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                  lagres,
                                data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_F2F_lr)
round(ci.lin(final_model_nb_F2F_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_F2F_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_F2F_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p > 0.05 -- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_F2F<- predict(final_model_nb_F2F_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_F2F$fit
stbp_nb <- nb_pred_F2F$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
  mutate(lockdown1 = 0,
         lagres = 0)
nb_pred_F2F_nolockdown <- predict(final_model_nb_F2F_lr, newdata = datanew2_nb, type = 'response', se.fit = TRUE, interval = 'confidence', population = 1)
predicted_vals_nb_F2F_noLockdown <-nb_pred_F2F_nolockdown$fit
stbp_nb_noLdn <- nb_pred_F2F_nolockdown$se.fit

# Combine predictions and convert from log odds to percentage reporting 
df_se_nb <- bind_cols(stbp = stbp_nb, stbp_noLdn= stbp_nb_noLdn, pred = pred_values, pred_noLdn = predicted_vals_nb_F2F_noLockdown) %>%
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

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_F2F_lr))
effect_of_lockdown_cons_F2F <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_F2F_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_F2F_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
F2Fcons_RRs_65plusyears <- bind_rows(effect_of_lockdown_cons_F2F, post_lockdown_migrant_status)
rownames(F2Fcons_RRs_65plusyears) <- NULL
write.csv(F2Fcons_RRs_65plusyears, "filepath/F2Fcons_RRs_65plusyears.csv")

# 2j_Segmented regression model: Phone consultations (not included in paper) ----

# 2k_0-15 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone, age_subcohort == '0-15 years')
full_model_nb_phone <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(dataset, !is.na(lockdown1)))
round(ci.lin(full_model_nb_phone,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p < 0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(dataset, !is.na(lockdown1)))
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
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
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

## Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
phonecons_RRs_0_15years <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status)
rownames(phonecons_RRs_0_15years) <- NULL
write.csv(phonecons_RRs_0_15years, "filepath/phonecons_RRs_0_15years.csv" )

# 2l_16-24 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone, age_subcohort == '16-24 years')
full_model_nb_phone <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p > 0.05 - accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
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
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
phonecons_RRs_16_24years <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status)
rownames(phonecons_RRs_16_24years) <- NULL
write.csv(phonecons_RRs_16_24years, "filepath/phonecons_RRs_16_24years.csv" )

# 2m_25-34 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone, age_subcohort == '25-34 years')
full_model_nb_phone <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p < 0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p > 0.05 -- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
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
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
phonecons_RRs_25_34years <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status)
rownames(phonecons_RRs_25_34years) <- NULL
write.csv(phonecons_RRs_25_34years, "filepath/phonecons_RRs_25_34years.csv" )

# 2n_35-49 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone, age_subcohort == '35-49 years')
full_model_nb_phone <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
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
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
phonecons_RRs_35_49years <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status)
rownames(phonecons_RRs_35_49years) <- NULL
write.csv(phonecons_RRs_35_49years, "filepath/phonecons_RRs_35_49years.csv" )

# 2o_50-64 years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone, age_subcohort == '50-64 years')
full_model_nb_phone <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
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
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
phonecons_RRs_50_64years <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status)
rownames(phonecons_RRs_50_64years) <- NULL
write.csv(phonecons_RRs_50_64years, "filepath/phonecons_RRs_50_64years.csv" )

# 2p_65+ years ----

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_agesubcohort_ap_phone, age_subcohort == '>=65 years')
full_model_nb_phone <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migrant_status*lockdown1,
                              data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_phone)
round(ci.lin(full_model_nb_phone,Exp=T),3)
lagresiduals_m <- lag(residuals(full_model_nb_phone)) %>% as.numeric()
res <- residuals(full_model_nb_phone, type = 'deviance')

pacf(res, lag = 368)
acf(res,lag = 368)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_phone, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p <0.005

# Merge lag residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(dataset, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, migrant_status))
dataset <- full_join(dataset, lagres_timing, by = c('migrant_status','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_phone_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migrant_status + lockdown1:migrant_status+ 
                                    lagres,
                                  data =  filter(dataset, !is.na(lockdown1)))
summary(final_model_nb_phone_lr)
round(ci.lin(final_model_nb_phone_lr,Exp=T),3)

# Check for autocorrelation
res2 <- residuals(final_model_nb_phone_lr, type = 'deviance')
pacf(res2, lag = 364)
acf(res2,lag = 364)
# Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_phone_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p > 0.05 -- accept null hypothesis that residuals are from a white noise series 

# Set up prediction dataset
prediction_data <- dataset
prediction_data$lagres <- 0
#  Make predictions
nb_pred_phone<- predict(final_model_nb_phone_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_phone$fit
stbp_nb <- nb_pred_phone$se.fit # the standard errors of the predicted means 

# Predict counterfactual situation
datanew2_nb <- dataset %>%
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
outcome_plot_nb <- bind_cols(df_se_nb, dataset) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% # to convert predictions to rate 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_phone_lr))
effect_of_lockdown_cons_phone <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates)) %>%
  filter(var %in% c('lockdown1', 'migrant_statusMigrant', 'lockdown1:migrant_statusMigrant')) 

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_phone_lr))
post_lockdown_migrant_status <- exp(estimable(final_model_nb_phone_lr, c('migrant_statusMigrant' = 1, 'lockdown1:migrant_statusMigrant' = 1), conf=.95))
post_lockdown_migrant_status <- post_lockdown_migrant_status %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'migrant_status+lockdown1:migrant_status')

# Combine RRs and save
phonecons_RRs_65plusyears <- bind_rows(effect_of_lockdown_cons_phone, post_lockdown_migrant_status)
rownames(phonecons_RRs_65plusyears) <- NULL
write.csv(phonecons_RRs_65plusyears, "filepath/phonecons_RRs_65plusyears.csv")

# 3_Check weekly counts of ethnicat6 factor levels to see if sufficient power to include ethnicity interaction -----
load(file = "filepath/covid_alldata_weekly_records_2015_2020.Rdata")
test <- head(covid_alldata_weekly_records_2015_2020)

size_ethnicat6groups_weeklydata <- covid_alldata_weekly_records_2015_2020 %>%
  group_by(studyweek, migrant_status, ethnicat6) %>%
  summarise(count = n())

min_ethnicat6_count <- size_ethnicat6groups_weeklydata  %>%
  group_by(migrant_status, ethnicat6) %>%
  slice_min(count) %>%
  mutate(min_or_max = 'min')

max_ethnicat6_count <- size_ethnicat6groups_weeklydata  %>%
  group_by(migrant_status, ethnicat6) %>%
  slice_max(count) %>%
  mutate(min_or_max = 'max')

min_max_ethnicat6_weekly_counts <- bind_rows(min_ethnicat6_count, max_ethnicat6_count) %>%
  arrange(migrant_status)
write.csv(min_max_ethnicat6_weekly_counts, "filepath/Ethnicat6_counts/min_max_ethnicat6_weekly_counts.csv" )

weekly_ethnicat6_counts <- ggplot(filter(size_ethnicat6groups_weeklydata, studyweek <310), mapping = aes(x = studyweek, y = count, color = ethnicat6)) +
  geom_line()+
  facet_wrap(~migrant_status)
weekly_ethnicat6_counts
ggsave('filepath/weekly_ethnicat6_counts.png')

# 4_Ethnicity-only model------

# 4a_Prepare data for interrupted time series ----

## Remove studyweek1 
England_2015_2020_aggregate_weekly_conscounts_ethnicity <- England_2015_2020_aggregate_weekly_conscounts_ethnicity %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
England_2015_2020_aggregate_weekly_conscounts_ethnicity_ap <- England_2015_2020_aggregate_weekly_conscounts_ethnicity
England_2015_2020_aggregate_weekly_conscounts_ethnicity_ap$lockdown1[England_2015_2020_aggregate_weekly_conscounts_ethnicity_ap$date %in% adjustment_period] <- NA

# 4b_Calculate IRs -----

## Without accounting for an adjustment period 
function_to_calculate_IR(England_2015_2020_aggregate_weekly_conscounts_ethnicity, 'England_weekly_conscounts_ethnicity_all')
function_to_calculate_IR_F2F(England_2015_2020_aggregate_weekly_conscounts_ethnicity, 'England_weekly_conscounts_ethnicity_F2F')
function_to_calculate_IR_phone(England_2015_2020_aggregate_weekly_conscounts_ethnicity, 'England_weekly_conscounts_ethnicity_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(England_2015_2020_aggregate_weekly_conscounts_ethnicity_ap, 'England_weekly_conscounts_ethnicity_ap_all')
function_to_calculate_IR_F2F(England_2015_2020_aggregate_weekly_conscounts_ethnicity_ap, 'England_weekly_conscounts_ethnicity_ap_F2F')
function_to_calculate_IR_phone(England_2015_2020_aggregate_weekly_conscounts_ethnicity_ap, 'England_weekly_conscounts_ethnicity_ap_phone')

# 4c_Plot IRs ----

## Create dataframe for the vertical line lockdown labels (i.e. the adjustment-to-restrictions period in the COVID collateral paper - 8th March to the 28th March)
lockdowns <- data.frame(date = as_date(c('2020-03-08', '2020-03-29')))

## All consultations
all_cons_IR_England_ethnicity <- ggplot() + 
  geom_line(data = filter(England_weekly_conscounts_ethnicity_all, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = ethnicat6))+
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'All consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))
all_cons_IR_England_ethnicity
ggsave('filepath/all_cons_IR_England_ethnicity.png')

## Face-to-face consulations 
F2F_cons_IR_England_ethnicity <-  ggplot() + 
  geom_line(data = filter(England_weekly_conscounts_ethnicity_F2F, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = ethnicat6))+
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Face-to-face consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))
F2F_cons_IR_England_ethnicity
ggsave('results/01_Consultations/its_analysis/F2F_cons_IR_England_ethnicity.png')

## Phone consultations 
phone_cons_IR_England_ethnicity <- ggplot() + 
  geom_line(data = filter(England_weekly_conscounts_ethnicity_phone, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), aes(x = date, y = incidence_rate, colour = ethnicat6))+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Telephone consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1))
phone_cons_IR_England_ethnicity
ggsave('filepath/phone_cons_IR_England_ethnicity.png')

# 4d_Segmented regression model: All consultations -----

## Fit model 
full_model_nb_all <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + ethnicat6*lockdown1,
                            data = filter(England_weekly_conscounts_ethnicity_ap_all, !is.na(lockdown1)))
round(ci.lin(full_model_nb_all,Exp=T),3)
extract_glm_results_allages(full_model_nb_all, 'ethnicity_all_incladjustmentperiod_nolaggedresiduals', England_weekly_conscounts_ethnicity_ap_all, 'lockdown1')
write.csv(ethnicity_all_incladjustmentperiod_nolaggedresiduals, "filepath/ethnicity_all_incladjustmentperiod_nolaggedresiduals.csv" )

## Calculate lagged residuals
lagresiduals_m <- lag(residuals(full_model_nb_all)) %>% as.numeric()
res <- residuals(full_model_nb_all, type = 'deviance')

## Check autocorrelation
pacf(res, lag = 368)
acf(res,lag = 368)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(full_model_nb_all, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p<0.05, so evidence of autocorrelation 

## Merge lagged residuals with dataframe
lagres_timing <- bind_cols('studyweek' = filter(England_weekly_conscounts_ethnicity_ap_all, !is.na(lockdown1)),
                           'lagres' = lagresiduals_m) 
lagres_timing <- lagres_timing %>%
  mutate_at('lagres',~replace(., is.na(.), 0))
lagres_timing <- lagres_timing %>%
  dplyr::select(c(studyweek, lagres, ethnicat6))
England_weekly_conscounts_ethnicity_ap_all <- full_join(England_weekly_conscounts_ethnicity_ap_all, lagres_timing, by = c('ethnicat6','studyweek'))

## Re-run final model (including lagged residuals)
final_model_nb_all_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ lockdown1*ethnicat6+ 
                                  lagres,
                                data =  filter(England_weekly_conscounts_ethnicity_ap_all, !is.na(lockdown1)))
round(ci.lin(final_model_nb_all_lr,Exp=T),3)
extract_glm_results_allages(final_model_nb_all_lr, 'ethnicity_all_incladjustmentperiod_laggedresiduals', England_weekly_conscounts_ethnicity_ap_all, 'lockdown1')
write.csv(ethnicity_all_incladjustmentperiod_laggedresiduals, "filepath/ethnicity_all_incladjustmentperiod_laggedresiduals.csv" )

## Re-check for autocorrelation
res2 <- residuals(final_model_nb_all_lr, type = 'response')
pacf(res2, lag = 312)
acf(res2,lag = 312)
### Formal test for autocorrelation epirhandbook.com/time-series-and-outbreak-detection.html
response_res <- residuals(final_model_nb_all_lr, type = 'response')
Box.test(response_res, type = 'Ljung-Box') # p >0.05 accept null hypothesis that residuals are from a white noise series

##  Make predictions
prediction_data <- England_weekly_conscounts_ethnicity_ap_all
prediction_data$lagres <- 0
nb_pred_all<- predict(final_model_nb_all_lr, newdata = prediction_data, type = 'response', population = 1, interval = 'confidence', se.fit = TRUE)
pred_values <- nb_pred_all$fit
stbp_nb <- nb_pred_all$se.fit # the standard errors of the predicted means 

## Predict counterfactual situation
datanew2_nb <- England_weekly_conscounts_ethnicity_ap_all %>%
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
outcome_plot_nb <- bind_cols(df_se_nb, England_weekly_conscounts_ethnicity_ap_all) %>%
  mutate(final_pred = pred/person_years) %>% # to convert predictions to rate 
  mutate(final_upp = upp/person_years) %>%
  mutate(final_low = low/person_years)%>%
  mutate(final_pred_noLdn = pred_noLdn/person_years) %>% 
  mutate(final_upp_noLdn = upp_noLdn/person_years) %>%
  mutate(final_low_noLdn = low_noLdn/person_years)


## Plot the results and save

all_cons_ethnicity_nb <- ggplot(filter(outcome_plot_nb, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), 
                                aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp, fill = ethnicat6, color = ethnicat6)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Consultations') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '1 month', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank())

all_cons_ethnicity_nb
ggsave('filepath/all_cons_ethnicity_nb.png')

all_cons_ethnicity_nb_grid <- ggplot(filter(outcome_plot_nb, date <= as.Date('2020-06-28') & date > as.Date('2019-12-29')), 
                                     aes(x = date, y = final_pred, ymin = final_low, ymax = final_upp)) +
  geom_line()+
  geom_ribbon(alpha = 0.4, colour = 0)+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Consultations') + 
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black', lwd = 1, linetype = 'dotted') + 
  scale_x_date(name = 'Date', date_breaks = '2 months', expand = c(0.05,0), minor_breaks = NULL, date_labels =  '%b')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        legend.title = element_blank()) +
  facet_wrap(~ethnicat6)

all_cons_ethnicity_nb_grid
ggsave('filepath/all_cons_ethnicity_nb_grid.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_all <- parameter_estimates[grep('ethnicat6|lockdown', rownames(parameter_estimates)), ]
effect_of_lockdown_cons_all$var <- rownames(effect_of_lockdown_cons_all)

# Get RR for the effect of migrant_status on outcome after lockdown (migrant_status+interaction of migrant_status:lockdown1)
names(coef(final_model_nb_all_lr))
# White non-British
post_lockdown_WNB <- exp(estimable(final_model_nb_all_lr, c('ethnicat6White non-British' = 1, 'lockdown1:ethnicat6White non-British' = 1), conf=.95))
post_lockdown_WNB <- post_lockdown_WNB %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'ethnicat6White non-British + lockdown1:ethnicat6White non-British')
# Mixed
post_lockdown_Mixed <- exp(estimable(final_model_nb_all_lr, c('ethnicat6Mixed/Multiple ethnic groups' = 1, 'lockdown1:ethnicat6Mixed/Multiple ethnic groups' = 1), conf=.95))
post_lockdown_Mixed <- post_lockdown_Mixed %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'ethnicat6Mixed + lockdown1:ethnicat6Mixed')
# Asian
post_lockdown_Asian <- exp(estimable(final_model_nb_all_lr, c('ethnicat6Asian/Asian British' = 1, 'lockdown1:ethnicat6Asian/Asian British' = 1), conf=.95))
post_lockdown_Asian <- post_lockdown_Asian %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'ethnicat6Asian + lockdown1:ethnicat6Asian')
# Black
post_lockdown_Black <- exp(estimable(final_model_nb_all_lr, c('ethnicat6Black/African/Caribbean/Black British' = 1, 'lockdown1:ethnicat6Black/African/Caribbean/Black British' = 1), conf=.95))
post_lockdown_Black <- post_lockdown_Black %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'ethnicat6Black + lockdown1:ethnicat6Black')
# Other
post_lockdown_Other <- exp(estimable(final_model_nb_all_lr, c('ethnicat6Other ethnic group' = 1, 'lockdown1:ethnicat6Other ethnic group' = 1), conf=.95))
post_lockdown_Other <- post_lockdown_Other %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'ethnicat6Other + lockdown1:ethnicat6Other')
# Unknown
post_lockdown_Unknown <- exp(estimable(final_model_nb_all_lr, c('ethnicat6Unknown' = 1, 'lockdown1:ethnicat6Unknown' = 1), conf=.95))
post_lockdown_Unknown <- post_lockdown_Unknown %>%
  dplyr::select(Estimate, Lower.CI, Upper.CI) %>%
  rename('2.5%' = Lower.CI,
         '97.5%' = Upper.CI,
         'exp(Est.)' = Estimate) %>%
  mutate(var = 'ethnicat6Unknown + lockdown1:ethnicat6Unknown')

# Combine RRs and save
allcons_ethnicity_RRs <- bind_rows(effect_of_lockdown_cons_all, post_lockdown_WNB, post_lockdown_Mixed, post_lockdown_Asian, post_lockdown_Black, post_lockdown_Other, post_lockdown_Unknown)
row.names(allcons_ethnicity_RRs) <- NULL
write.csv(allcons_ethnicity_RRs, "filepath/allcons_ethnicity_RRs.csv" )

# LRT to check overall significance of ethnicat6 coefficient
## No ethnicat6 in model
### With lagged residuals
final_model_nb_all_lr_noethnicat6 <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ 
                                              lagres,
                                            data =  filter(England_weekly_conscounts_ethnicity_ap_all, !is.na(lockdown1)))
lrtest(final_model_nb_all_lr_noethnicat6, final_model_nb_all_lr) # p<0.001

### Without lagged residuals
full_model_nb_all_noethnicat6 <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth),
                                        data =  filter(England_weekly_conscounts_ethnicity_ap_all, !is.na(lockdown1)))
lrtest(full_model_nb_all_noethnicat6, full_model_nb_all) # p<0.001

# 5_SENSITIVITY ANALYSIS: migcertainty age_subcohort ----

# Remove studyweek1 (DELETE when we have the new file)
England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

# Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 2nd to 29th March 2020

adjustment_period <- c(as_date('2020-03-02'), as_date('2020-03-09'), as_date('2020-03-16'), as_date('2020-03-23'))
covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap <- England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort
covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap$lockdown1[covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap$date %in% adjustment_period] <- NA

# Calculate IRs
# all data 
function_to_calculate_IR(England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort, 'covid_weekly_conscounts_England_migcertainty_agesubcohort_all')
function_to_calculate_IR_F2F(England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort, 'covid_weekly_conscounts_England_migcertainty_agesubcohort_F2F')
function_to_calculate_IR_phone(England_2015_2020_aggregate_weekly_conscounts_migcertainty_agesubcohort, 'covid_weekly_conscounts_England_migcertainty_agesubcohort_phone')
# truncated data with lockdown1 set to NA during the adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap, 'covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_phone')

# 5a_Plot of pop sizes for each age subcohort ------
# All consultations

popsize_total_weekly <- covid_weekly_conscounts_England_migcertainty_agesubcohort_all %>%
  group_by(migcertainty, studyweek) %>%
  summarise(popsize_total = sum(popsize))

test <- covid_weekly_conscounts_England_migcertainty_agesubcohort_all %>%
  left_join(popsize_total_weekly, by = c('migcertainty', 'studyweek')) %>%
  mutate(percentage_popsize = popsize/popsize_total*100)

lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))

migcertainty_pop_size_agesubcohort <- ggplot(filter(test, migcertainty %in% c('Definite', 'Probable')), 
                                             aes(x = date, y = percentage_popsize)) +geom_point(aes(colour = age_subcohort))+
  labs(y = 'Population size', x = 'Year', 
       title = 'Population size broken down into age subcohorts') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap(vars(migcertainty))
migcertainty_pop_size_agesubcohort
ggsave('filepath/migcertainty_pop_size_agesubcohort.png')

# Plot IRs 

#### All data
# All consultations
all_cons_IR_England_migcertainty_agesubcohort <- ggplot(covid_weekly_conscounts_England_migcertainty_agesubcohort_all, aes(x = date, y = incidence_rate)) +geom_point(aes(colour = migcertainty))+ 
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'All consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(rows = vars(age_subcohort))
all_cons_IR_England_migcertainty_agesubcohort
ggsave('filepath/all_cons_IR_England_migcertainty_agesubcohort.png')

# Face-to-face consulations 
F2F_cons_IR_England_migcertainty_agesubcohort <- ggplot(covid_weekly_conscounts_England_migcertainty_agesubcohort_F2F, aes(x = date, y = incidence_rate)) +geom_point(aes(col = migcertainty))+
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Face-to-face consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(rows = vars(age_subcohort))
F2F_cons_IR_England_migcertainty_agesubcohort
ggsave('filepath/F2F_cons_IR_England_migcertainty_agesubcohort.png')

# Phone consultations 
phone_cons_IR_England_migcertainty_agesubcohort <- ggplot(covid_weekly_conscounts_England_migcertainty_agesubcohort_phone, aes(x = date, y = incidence_rate)) +geom_point(aes(colour = migcertainty))+
  labs(y = 'Consultation rate (per person-year)', x = 'Year', 
       title = 'Telephone consultations') +
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  scale_colour_discrete(name = 'Legend')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(rows = vars(age_subcohort))
phone_cons_IR_England_migcertainty_agesubcohort
ggsave('filepath/phone_cons_IR_England_migcertainty_agesubcohort.png')

# 5b_All consultations ----

test <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1 + age_subcohort,
               data = filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all, !is.na(lockdown1)))
summary(test)
round(ci.lin(test,Exp=T),3)

# 0-15 years

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all, age_subcohort == '0-15 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p = 0.6291-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


all_cons_nb_migcertainty_0_15years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'All consultations: 0-15 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
all_cons_nb_migcertainty_0_15years
ggsave('filepath/all_cons_nb_migcertainty_0_15years.png')

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
allcons_RRs_migcertainty_0to15years <- bind_rows(effect_of_lockdown_cons_all_migcertainty,
                                                 post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(allcons_RRs_migcertainty_0to15years) = NULL
write.csv(allcons_RRs_migcertainty_0to15years, "filepath/allcons_RRs_migcertainty_0to15years.csv" )

# Format for tables markdown
allcons_RRs_migcertainty_0to15years <- allcons_RRs_migcertainty_0to15years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group0to15 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group0to15) 

# 16-24 years

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all, age_subcohort == '16-24 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


all_cons_nb_migcertainty_16_24years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'All consultations: 16-24 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
all_cons_nb_migcertainty_16_24years
ggsave('filepath/all_cons_nb_migcertainty_16_24years.png')

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
allcons_RRs_migcertainty_16to24years <- bind_rows(effect_of_lockdown_cons_all_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(allcons_RRs_migcertainty_16to24years) = NULL
write.csv(allcons_RRs_migcertainty_16to24years, "filepath/allcons_RRs_migcertainty_16to24years.csv")

# Format for tables markdown
allcons_RRs_migcertainty_16to24years <- allcons_RRs_migcertainty_16to24years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group16to24 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group16to24) 

# 25-34 years

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all, age_subcohort == '25-34 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


all_cons_nb_migcertainty_25_34years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'All consultations: 25-34 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
all_cons_nb_migcertainty_25_34years
ggsave('filepath/all_cons_nb_migcertainty_25_34years.png')

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
allcons_RRs_migcertainty_25to34years <- bind_rows(effect_of_lockdown_cons_all_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(allcons_RRs_migcertainty_25to34years) = NULL
write.csv(allcons_RRs_migcertainty_25to34years, "filepath/allcons_RRs_migcertainty_25to34years.csv")

# Format for tables markdown
allcons_RRs_migcertainty_25to34years <- allcons_RRs_migcertainty_25to34years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group25to34 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group25to34) 

# 35-49 years

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all, age_subcohort == '35-49 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


all_cons_nb_migcertainty_35_49years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'All consultations: 35-49 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
all_cons_nb_migcertainty_35_49years
ggsave('filepath/all_cons_nb_migcertainty_35_49years.png')

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
allcons_RRs_migcertainty_35to49years <- bind_rows(effect_of_lockdown_cons_all_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(allcons_RRs_migcertainty_35to49years) = NULL
write.csv(allcons_RRs_migcertainty_35to49years, "filepath/allcons_RRs_migcertainty_35to49years.csv")

# Format for tables markdown
allcons_RRs_migcertainty_35to49years <- allcons_RRs_migcertainty_35to49years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group35to49 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group35to49) 

# 50-64 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all, age_subcohort == '50-64 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


all_cons_nb_migcertainty_50_64years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'All consultations: 50-64 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
all_cons_nb_migcertainty_50_64years
ggsave('filepath/all_cons_nb_migcertainty_50_64years.png')

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
allcons_RRs_migcertainty_50to64years <- bind_rows(effect_of_lockdown_cons_all_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(allcons_RRs_migcertainty_50to64years) = NULL
write.csv(allcons_RRs_migcertainty_50to64years, "filepath/allcons_RRs_migcertainty_50to64years.csv")

# Format for tables markdown
allcons_RRs_migcertainty_50to64years <- allcons_RRs_migcertainty_50to64years %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group50to64 = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group50to64) 

# 65+ years

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_all, age_subcohort == '>=65 years')
full_model_nb_all <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
                            data = filter(dataset, !is.na(lockdown1)))
summary(full_model_nb_all)
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(conscount ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p >0.05-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


all_cons_nb_migcertainty_65plusyears <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'All consultations: 65 years and above', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
all_cons_nb_migcertainty_65plusyears
ggsave('filepath/all_cons_nb_migcertainty_65plusyears.png')

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
allcons_RRs_migcertainty_65plusyears <- bind_rows(effect_of_lockdown_cons_all_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(allcons_RRs_migcertainty_65plusyears) = NULL
write.csv(allcons_RRs_migcertainty_65plusyears, "filepath/allcons_RRs_migcertainty_65plusyears.csv")

# Format for tables markdown
allcons_RRs_migcertainty_65plusyears <- allcons_RRs_migcertainty_65plusyears %>%
  mutate(across(where(is.numeric), ~ round(.,2))) %>%
  mutate(group65plus = paste0(`exp(Est.)`, ' (', `2.5%`, '-', `97.5%`, ')')) %>%
  dplyr::select(var, group65plus) 

# 5c_Combine RRs for rmarkdown and save as Rdata file ----
RRs_agesubcohort_migcertainty <- bind_cols(allcons_RRs_migcertainty_0to15years, allcons_RRs_migcertainty_16to24years, allcons_RRs_migcertainty_25to34years,
                                           allcons_RRs_migcertainty_35to49years, allcons_RRs_migcertainty_50to64years, allcons_RRs_migcertainty_65plusyears) %>%
  dplyr::select(`var...1`, group0to15, group16to24, group25to34, group35to49, group50to64, group65plus) %>%
  rename('Variable' = `var...1`)
save(RRs_agesubcohort_migcertainty, file = 'filepath/RRs_agesubcohort_migcertainty.Rdata')


# 5d_Face-to-face consultations (not included in paper) ----

# 0-14 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_F2F, age_subcohort == '0-14 years')
full_model_nb_all <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p = 0.6291-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


F2F_cons_nb_migcertainty_0_14years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Face-to-face consultations: 0-14 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
F2F_cons_nb_migcertainty_0_14years
ggsave('filepath/F2F_cons_nb_migcertainty_0_14years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_F2F_migcertainty <- parameter_estimates %>%
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
F2Fcons_RRs_migcertainty_0to14years <- bind_rows(effect_of_lockdown_cons_F2F_migcertainty,
                                                 post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(F2Fcons_RRs_migcertainty_0to14years) = NULL
write.csv(F2Fcons_RRs_migcertainty_0to14years, "filepath/F2Fcons_RRs_migcertainty_0to14years.csv" )

# 15-24 years

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_F2F, age_subcohort == '15-24 years')
full_model_nb_all <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p = 0.6291-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


F2F_cons_nb_migcertainty_15_24years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Face-to-face consultations: 15-24 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
F2F_cons_nb_migcertainty_15_24years
ggsave('results/01_Consultations/its_analysis/F2F_cons_nb_migcertainty_15_24years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_F2F_migcertainty <- parameter_estimates %>%
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
F2Fcons_RRs_migcertainty_15to24years <- bind_rows(effect_of_lockdown_cons_F2F_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(F2Fcons_RRs_migcertainty_15to24years) = NULL
write.csv(F2Fcons_RRs_migcertainty_15to24years, "filepath/F2Fcons_RRs_migcertainty_15to24years.csv")

# 25-44 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_F2F, age_subcohort == '25-44 years')
full_model_nb_all <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p = 0.6291-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


F2F_cons_nb_migcertainty_25_44years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Face-to-face consultations: 25-44 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
F2F_cons_nb_migcertainty_25_44years
ggsave('filepath/F2F_cons_nb_migcertainty_25_44years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_F2F_migcertainty <- parameter_estimates %>%
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
F2Fcons_RRs_migcertainty_25to44years <- bind_rows(effect_of_lockdown_cons_F2F_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(F2Fcons_RRs_migcertainty_25to44years) = NULL
write.csv(F2Fcons_RRs_migcertainty_25to44years, "filepath/F2Fcons_RRs_migcertainty_25to44years.csv")

# 45-64 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_F2F, age_subcohort == '45-64 years')
full_model_nb_all <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p = 0.6291-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


F2F_cons_nb_migcertainty_45_64years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Face-to-face consultations: 45-64 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
F2F_cons_nb_migcertainty_45_64years
ggsave('filepath/F2F_cons_nb_migcertainty_45_64years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_F2F_migcertainty <- parameter_estimates %>%
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
F2Fcons_RRs_migcertainty_45to64years <- bind_rows(effect_of_lockdown_cons_F2F_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(F2Fcons_RRs_migcertainty_45to64years) = NULL
write.csv(F2Fcons_RRs_migcertainty_45to64years, "filepath/F2Fcons_RRs_migcertainty_45to64years.csv")

# 65+ years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_F2F, age_subcohort == '>=65 years')
full_model_nb_all <- glm.nb(events ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(events ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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
Box.test(response_res, type = 'Ljung-Box') # p = 0.6291-- accept null hypothesis that residuals are from a white noise series 

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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


F2F_cons_nb_migcertainty_65plusyears <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Face-to-face consultations: 65 years and above', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
F2F_cons_nb_migcertainty_65plusyears
ggsave('filepath/F2F_cons_nb_migcertainty_65plusyears.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_F2F_migcertainty <- parameter_estimates %>%
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
F2Fcons_RRs_migcertainty_65plusyears <- bind_rows(effect_of_lockdown_cons_F2F_migcertainty,
                                                  post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(F2Fcons_RRs_migcertainty_65plusyears) = NULL
write.csv(F2Fcons_RRs_migcertainty_65plusyears, "filepath/F2Fcons_RRs_migcertainty_65plusyears.csv")

# 5e_Phone consultations (not included in paper) ----

# 0-14 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_phone, age_subcohort == '0-14 years')
full_model_nb_all <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


phone_cons_nb_migcertainty_0to14years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Phone consultations: 0-14 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
phone_cons_nb_migcertainty_0to14years
ggsave('filepath/phone_cons_nb_migcertainty_0to14years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_phone_migcertainty <- parameter_estimates %>%
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
phonecons_RRs_migcertainty_0to14years <- bind_rows(effect_of_lockdown_cons_phone_migcertainty,
                                                   post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(phonecons_RRs_migcertainty_0to14years) = NULL
write.csv(phonecons_RRs_migcertainty_0to14years, "filepath/phonecons_RRs_migcertainty_0to14years.csv")

# 15-24 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_phone, age_subcohort == '15-24 years')
full_model_nb_all <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


phone_cons_nb_migcertainty_15to24years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Phone consultations: 15-24 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
phone_cons_nb_migcertainty_15to24years
ggsave('filepath/phone_cons_nb_migcertainty_15to24years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_phone_migcertainty <- parameter_estimates %>%
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
phonecons_RRs_migcertainty_15to24years <- bind_rows(effect_of_lockdown_cons_phone_migcertainty,
                                                    post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(phonecons_RRs_migcertainty_15to24years) = NULL
write.csv(phonecons_RRs_migcertainty_15to24years, "filepath/phonecons_RRs_migcertainty_15to24years.csv")

# 25-44 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_phone, age_subcohort == '25-44 years')
full_model_nb_all <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


phone_cons_nb_migcertainty_25to44years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Phone consultations: 25-44 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
phone_cons_nb_migcertainty_25to44years
ggsave('filepath/phone_cons_nb_migcertainty_25to44years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_phone_migcertainty <- parameter_estimates %>%
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
phonecons_RRs_migcertainty_25to44years <- bind_rows(effect_of_lockdown_cons_phone_migcertainty,
                                                    post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(phonecons_RRs_migcertainty_25to44years) = NULL
write.csv(phonecons_RRs_migcertainty_25to44years, "filepath/phonecons_RRs_migcertainty_25to44years.csv")

# 45-64 years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_phone, age_subcohort == '45-64 years')
full_model_nb_all <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


phone_cons_nb_migcertainty_45to64years <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Phone consultations: 45-64 years', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
phone_cons_nb_migcertainty_45to64years
ggsave('filepath/phone_cons_nb_migcertainty_45to64years.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_phone_migcertainty <- parameter_estimates %>%
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
phonecons_RRs_migcertainty_45to64years <- bind_rows(effect_of_lockdown_cons_phone_migcertainty,
                                                    post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(phonecons_RRs_migcertainty_45to64years) = NULL
write.csv(phonecons_RRs_migcertainty_45to64years, "filepath/phonecons_RRs_migcertainty_45to64years.csv")

# 65+ years 

# Fit model and calculate lagged residuals
dataset <- filter(covid_weekly_conscounts_aggregate_England_migcertainty_agesubcohort_ap_phone, age_subcohort == '>=65 years')
full_model_nb_all <- glm.nb(phone ~ offset(log(person_years)) + lockdown1+ studyweek + as.factor(studymonth) + migcertainty*lockdown1,
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
  dplyr::select(c(studyweek, lagres, migcertainty))
dataset <- full_join(dataset, lagres_timing, by = c('migcertainty','studyweek'))

# Run final model (incl. lagged residuals)
final_model_nb_all_lr <- glm.nb(phone ~ offset(log(person_years)) + lockdown1 + studyweek+ as.factor(studymonth)+ migcertainty + lockdown1:migcertainty+ 
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


lockdowns <- data.frame(date = as_date(c('2020-03-30', '2020-06-28')), event = c('Lockdown 1 start', 'Lockdown 1 end'))
#lockdowns_1_2 <- data.frame(date = as_date(c('2020-03-30', '2020-06-28', '2020-11-05')), event = c('Lockdown 1 start', 'Lockdown 1 end', 'Lockdown 2 start'))


phone_cons_nb_migcertainty_65plusyears <- ggplot(outcome_plot_nb, mapping = aes(x = date, y = incidence_rate)) +
  geom_point(col = 'gray60')+
  # the probability if there was no lockdown 
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), mapping = aes(y = final_pred, colour = 'Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), mapping = aes(y = final_pred, colour = 'Migrants'))+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Non-migrant'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Non-migrant'), lty = 0) +
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Definite'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Definite'), lty = 0)+
  geom_ribbon(data = filter(outcome_plot_nb, migcertainty == 'Probable'), mapping = aes(ymin = final_low, ymax = final_upp, fill = 'Probable'), lty = 0)+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Non-migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Non-migrants'))+
  #geom_line(data = filter(outcome_plot_nb, migrant_status == 'Migrant'), aes(y = final_pred_noLdn, color = 'Counterfactual scenario - Migrants')) +
  labs(y = 'Consultation rate (per person-year)', x = 'Year', title = 'Phone consultations: 65 years and over', color = 'Legend') + 
  scale_fill_manual(name = 'Legend', values = c('Non-migrant' = '#F8766D','Definite' = '#00BFC4', 'Probable' = '#7CAE00'), breaks = c('Non-migrant', 'Definite', 'Probable'))+
  geom_vline(data = lockdowns, mapping = aes(xintercept = date), color = 'black') + 
  geom_text(data =lockdowns, mapping = aes(x=date, y = 0, label = event), size = 4, angle = 90, vjust = -0.4, hjust =0)+
  scale_x_date(name = 'Date', date_breaks = '3 month', limits = as.Date(c('01/01/2015', '01/01/2021'), format = '%d/%m/%Y'), expand = c(0,0), minor_breaks = NULL, date_labels =  '%b-%Y')+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
phone_cons_nb_migcertainty_65plusyears
ggsave('filepath/phone_cons_nb_migcertainty_65plusyears.png')

# Get rate ratios 

# Get RRs for effect of lockdown
parameter_estimates <- as.data.frame(ci.exp(final_model_nb_all_lr))
effect_of_lockdown_cons_phone_migcertainty <- parameter_estimates %>%
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
phonecons_RRs_migcertainty_65plusyears <- bind_rows(effect_of_lockdown_cons_phone_migcertainty,
                                                    post_lockdown_migcertaintyDefinite, post_lockdown_migcertaintyProbable)
row.names(phonecons_RRs_migcertainty_65plusyears) = NULL
write.csv(phonecons_RRs_migcertainty_65plusyears, "filepath/phonecons_RRs_migcertainty_65plusyears.csv")

# 6_Model choice -----

# Poisson vs negative binomial --

## Remove studyweek1 (as it's not a full week)
covid_weekly_conscounts_aggregate_England <- England_2015_2020_aggregate_weekly_conscounts %>%
  filter(studyweek >1) %>%
  mutate(studyweek = studyweek -1)

## Create dataset to account for the lockdown1 adjustment period (i.e. change lockdown1 values to NA from 8th to 28th March 2020 [Mansfield et al. 2021])
adjustment_period <- c(as_date('2020-03-08'), as_date('2020-03-15'), as_date('2020-03-22'))
covid_weekly_conscounts_aggregate_England_ap <- covid_weekly_conscounts_aggregate_England
covid_weekly_conscounts_aggregate_England_ap$lockdown1[covid_weekly_conscounts_aggregate_England_ap$date %in% adjustment_period] <- NA

# 4b_Calculate IRs

## Without accounting for an adjustment period 
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England, 'covid_weekly_conscounts_phone')
# Accounting for an adjustment period (8th-28th March)
function_to_calculate_IR(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_all')
function_to_calculate_IR_F2F(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_F2F')
function_to_calculate_IR_phone(covid_weekly_conscounts_aggregate_England_ap, 'covid_weekly_conscounts_aggregate_England_ap_phone')

# Poisson
poisson_model <- glm(conscount ~ offset(log(person_years)) + as.factor(studymonth) + migrant_status*lockdown1, family = 'poisson',
                     data = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))

# Negative binomial 

nb_model <- glm.nb(conscount ~ offset(log(person_years)) + as.factor(studymonth) + migrant_status*lockdown1,
                   data = filter(covid_weekly_conscounts_aggregate_England_ap_all, !is.na(lockdown1)))

lrtest(poisson_model, nb_model) #p<0.001 NB better

