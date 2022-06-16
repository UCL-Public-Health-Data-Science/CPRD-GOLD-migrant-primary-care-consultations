
####---- Description -------------------------------------------------------------------------

## Forest plots for consultations annual analyses
## Date started: 02/07/2021
## Authors: Yamina Boukari / Claire Zhang
## QC (date): Claire Zhang (05/03/2022)

# 01_Load packages -------------------------------------------------------------------------

if (!require('pacman')) install.packages('pacman')
pacman::p_load(tidyverse,lubridate, Hmisc, epitools, data.table, epiR, MASS, psych, forestplot, gtable, gridExtra)

# 02_Set working directory ------------------------------------------------------------------

setwd("filepath")


# 03_FOREST PLOT 1 - MAIN ANALYSIS  -----------------

# England 

load(file = "filepath") 
full_cons_glm_unadj <- glm_mig %>% filter(names=="Migrant")
full_cons_glm_unadj$names <- recode(full_cons_glm_unadj$names,
                                    "Migrant" = "Unadjusted")

load(file = "filepath") 
full_cons_glm_adj_noIMD <- multivariable_allages %>% filter(names=="Migrant") 
full_cons_glm_adj_noIMD$names <- recode(full_cons_glm_adj_noIMD$names,
                                        "Migrant" = "Adjusted, without IMD")

load(file = "filepath")
full_cons_glm_adj_withIMD <- multivariable_allages %>% filter(names=="Migrant") 
full_cons_glm_adj_withIMD$names <- recode(full_cons_glm_adj_withIMD$names,
                                          "Migrant" = "Adjusted, with IMD")

load(file = "filepath")
full_cons_glm_adj_withIMD_0to15 <- glm_multivariable_0to15 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_0to15$names <- recode(full_cons_glm_adj_withIMD_0to15$names,
                                                "Migrant" = "    0-15 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_16to24 <- glm_multivariable_16to24 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_16to24$names <- recode(full_cons_glm_adj_withIMD_16to24$names,
                                                 "Migrant" = "    16-24 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_25to34 <- glm_multivariable_25to34 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_25to34$names <- recode(full_cons_glm_adj_withIMD_25to34$names,
                                                 "Migrant" = "    25-34 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_35to49 <- glm_multivariable_35to49 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_35to49$names <- recode(full_cons_glm_adj_withIMD_35to49$names,
                                                 "Migrant" = "    35-49 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_50to64 <- glm_multivariable_50to64 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_50to64$names <- recode(full_cons_glm_adj_withIMD_50to64$names,
                                                 "Migrant" = "    50-64 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_65plus <- glm_multivariable_65plus %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_65plus$names <- recode(full_cons_glm_adj_withIMD_65plus$names, 
                                                 "Migrant" = "    65 years and older")

## Join all subgroups 

England_cohort_plot_data <- bind_rows(full_cons_glm_unadj, full_cons_glm_adj_noIMD, full_cons_glm_adj_withIMD,
                                   full_cons_glm_adj_withIMD_0to15, full_cons_glm_adj_withIMD_16to24,
                                   full_cons_glm_adj_withIMD_25to34, full_cons_glm_adj_withIMD_35to49,
                                   full_cons_glm_adj_withIMD_50to64,
                                   full_cons_glm_adj_withIMD_65plus) %>%
  add_row(.before = 4, names = 'Adjusted, with IMD, by age group')

## London only  -------

load(file = "filepath")
london_cons_glm_unadj <- London_glm_mig %>% filter(names=="Migrant") 
london_cons_glm_unadj$names <- recode(london_cons_glm_unadj$names,
                                      "Migrant" = "Unadjusted")

load(file = "filepath")
london_cons_glm_adj_withoutIMD <- London_multivariable_allages_noIMD %>% filter(names=="Migrant") 
london_cons_glm_adj_withoutIMD$names <- recode(london_cons_glm_adj_withoutIMD$names,
                                               "Migrant" = "Adjusted, without IMD")

load(file = "filepath")
london_cons_glm_adj_withIMD <- London_multivariable_allages %>% filter(names=="Migrant") 
london_cons_glm_adj_withIMD$names <- recode(london_cons_glm_adj_withIMD$names,
                                            "Migrant" = "Adjusted, with IMD")

load(file = "filepath")
london_cons_glm_adj_withIMD_0to15 <- London_glm_multivariable_0to15 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
london_cons_glm_adj_withIMD_0to15$names <- recode(london_cons_glm_adj_withIMD_0to15$names,
                                                  "Migrant" = "    0-15 years")

load(file = "filepath")
london_cons_glm_adj_withIMD_16to24 <- London_glm_multivariable_16to24 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
london_cons_glm_adj_withIMD_16to24$names <- recode(london_cons_glm_adj_withIMD_16to24$names,
                                                   "Migrant" = "    16-24 years")

load(file = "filepath")
london_cons_glm_adj_withIMD_25to34 <- London_glm_multivariable_25to34 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
london_cons_glm_adj_withIMD_25to34$names <- recode(london_cons_glm_adj_withIMD_25to34$names,
                                                   "Migrant" = "    25-34 years")

load(file = "filepath")
london_cons_glm_adj_withIMD_35to49 <- London_glm_multivariable_35to49 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
london_cons_glm_adj_withIMD_35to49$names <- recode(london_cons_glm_adj_withIMD_35to49$names,
                                                   "Migrant" = "    35-49 years")

load(file = "filepath")
london_cons_glm_adj_withIMD_50to64 <- London_glm_multivariable_50to64 %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
london_cons_glm_adj_withIMD_50to64$names <- recode(london_cons_glm_adj_withIMD_50to64$names,
                                                   "Migrant" = "    50-64 years")

load(file = "filepath")
london_cons_glm_adj_withIMD_65plus <- London_glm_multivariable_65plus %>% filter(names=="Migrant") %>% dplyr::select(-c(age_subcohort))
london_cons_glm_adj_withIMD_65plus$names <- recode(london_cons_glm_adj_withIMD_65plus$names, 
                                                   "Migrant" = "    65 years and older")

london_cohort_plot_data <- bind_rows(london_cons_glm_unadj, london_cons_glm_adj_withoutIMD, london_cons_glm_adj_withIMD,
                                     london_cons_glm_adj_withIMD_0to15, london_cons_glm_adj_withIMD_16to24,
                                     london_cons_glm_adj_withIMD_25to34, london_cons_glm_adj_withIMD_35to49,
                                     london_cons_glm_adj_withIMD_50to64,london_cons_glm_adj_withIMD_65plus) %>%
  add_row(.before = 4, names = 'Adjusted, with IMD by age group') 


## Combine England and London into wide format for forest plot -----

forest_plot_main_all <- England_cohort_plot_data %>%
  dplyr::select(-c(p, irr_ci))
forest_plot_main_all <- forest_plot_main_all %>%
  rename(estimateEngland = estimate,
         lowerEngland = lower,
         upperEngland = upper,
         ciEngland = ci)

forest_plot_main_all$estimateLondon <- london_cohort_plot_data$estimate
forest_plot_main_all$lowerLondon <- london_cohort_plot_data$lower
forest_plot_main_all$upperLondon <- london_cohort_plot_data$upper
forest_plot_main_all$ciLondon <- london_cohort_plot_data$ci

## Create forest plot -----
# Function to make England labels above London labels 
sfrac <- function(top, bottom, data=NULL){
  with(data,lapply(paste0('atop(', top, ',', bottom, ')'), str2expression))}

tabletext <- list(
  c("Migrants vs. non-migrants", forest_plot_main_all$names),
  c("RR", sfrac(forest_plot_main_all$estimateEngland, forest_plot_main_all$estimateLondon, data = forest_plot_main_all)),
  c("95% CI", sfrac(forest_plot_main_all$ciEngland, forest_plot_main_all$ciLondon, data = forest_plot_main_all)))

dev.new()
png(file = "filepath", width = 2500, height = 2000)
forest_plot_main <- forestplot(tabletext, 
                             mean = cbind(c(NA, forest_plot_main_all$estimateEngland), c(NA, forest_plot_main_all$estimateLondon)),
                             lower = cbind(c(NA, forest_plot_main_all$lowerEngland), c(NA, forest_plot_main_all$lowerLondon)),
                             upper = cbind(c(NA, forest_plot_main_all$upperEngland), c(NA, forest_plot_main_all$upperLondon)),
                             new_page = TRUE,
                             zero = 1,
                             lineheight = unit(50, 'mm'),
                             line.margin = .1,
                             xlog = TRUE, 
                             xlab = 'Rate ratio with 95% CI',
                             col = fpColors(box = c('red4', 'skyblue3'),
                                            lines = c('red3', 'skyblue2')),
                             fn.ci.norm = c(fpDrawNormalCI, fpDrawDiamondCI),
                             is.summary = c(TRUE, rep(FALSE, 10)),
                             boxsize = 0.1,
                             xticks = c(0.5, 1.0, 1.5),
                             legend = c('England', 'London'),
                             ci.vertices = TRUE,
                             ci.vertices.height = 0.1, 
                             colgap = unit(10,'mm'),
                             graphwidth = unit(300, 'mm'),
                             graph.pos =2,
                             alim = c(0.5, 1.5),
                             align = c('l', 'c', 'c', 'c'),
                             title = ' ',
                             txt_gp = fpTxtGp(ticks=gpar(cex=4), xlab=gpar(cex=4), cex = 4, summary = gpar(fontface = 'bold'),
                                              label = list(gpar(fontface = 'plain'),
                                                           gpar(fontface = 'bold', cex = 4),
                                                           gpar(fontface = 'bold', cex =4))))

dev.off()


# 04_FOREST PLOT 2 - Ethnicity analysis ----
## Migrants vs White British non-migrants (adjusted, with IMD)------
load(file = 'filepath')
load(file = "filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")

ethnicity_vsWBNM <- bind_rows(white_migrant_vsWBNM, white_nonbritish_migrant_vsWBNM, mixed_migrant_vsWBNM, 
                              asian_migrant_vsWBNM, black_migrant_vsWBNM, other_migrant_vsWBNM,unknown_migrant_vsWBNM, all_migrant_vsWBNM)
names <- c('White British', 'White non-British', 'Mixed/Multiple ethnic groups', 'Asian/Asian British', 'Black/African/Caribbean/Black British',
           'Other ethnic group', 'Unknown', 'All ethnic groups*')
ethnicity_vsWBNM$names <- names
ethnicity_vsWBNM <- ethnicity_vsWBNM %>%
  mutate(p = ' ') %>%
  dplyr::select(names, Estimate, Lower.CI, Upper.CI, p, ci, RR) %>%
  rename(estimate = Estimate,
         lower = Lower.CI,
         upper = Upper.CI,
         irr_ci = RR)

row.names(ethnicity_vsWBNM) <- NULL

## Migrants vs. non-migrants of the same ethnicity (adjusted, with IMD) ----

load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")
load(file="filepath")

ethnicity_vsNM <- bind_rows(white_strata_vsNM, white_nonbritish_strata_vsNM, mixed_strata_vsNM, 
                            asian_strata_vsNM, black_strata_vsNM, other_strata_vsNM,unknown_strata_vsNM)
names <- c('White British', 'White non-British', 'Mixed/Multiple ethnic groups', 'Asian/Asian British', 'Black/African/Caribbean/Black British',
           'Other ethnic group', 'Unknown')
ethnicity_vsNM$names <- names
ethnicity_vsNM <- ethnicity_vsNM %>%
  mutate(p = ' ') %>%
  dplyr::select(names, Estimate, Lower.CI, Upper.CI, p, ci, RR) %>%
  rename(estimate = Estimate,
         lower = Lower.CI,
         upper = Upper.CI,
         irr_ci = RR) %>%
  add_row(.before = 1, names = 'Migrants vs non-migrants of the same ethnicity')

row.names(ethnicity_vsNM) <- NULL

ethnicity_plot_data <- bind_rows(ethnicity_vsWBNM, ethnicity_vsNM)
ethnicity_plot_data$p <- as.numeric(ethnicity_plot_data$p)

## Make forest plot using forest plot function -----
tabletext <- cbind(c("Migrants vs. White British non-migrants", ethnicity_plot_data$names),
                   c("RR", ethnicity_plot_data$estimate),
                   c("95% CI", ethnicity_plot_data$ci))

dev.new()
png(file = "filepath", width = 3000, height = 2600)
forest_plots_main_analysis_ethnicity <- forestplot(tabletext, 
                             mean = cbind(c(NA, ethnicity_plot_data$estimate)),
                             lower = cbind(c(NA, ethnicity_plot_data$lower)),
                             upper = cbind(c(NA, ethnicity_plot_data$upper)),
                             new_page = TRUE,
                             is.summary = c(TRUE, rep(FALSE, 8), TRUE, rep(FALSE, 7)),
                             zero = 1,
                             lineheight = unit(50, 'mm'),
                             line.margin = .1,
                             xlog = TRUE, 
                             xlab = 'Rate ratio with 95% CI',
                             col = fpColors(box = 'red4',
                                            lines = 'red3'),
                             fn.ci.norm = c(fpDrawNormalCI, fpDrawDiamondCI),
                             boxsize = 0.1,
                             xticks = c(0.5, 1.0, 1.5),
                             ci.vertices = TRUE,
                             ci.vertices.height = 0.1, 
                             colgap = unit(10,'mm'),
                             graphwidth = unit(300, 'mm'),
                             graph.pos =2,
                             alim = c(0.5, 1.5),
                             align = c('l', 'c', 'c', 'c'),
                             title = ' ',
                             txt_gp = fpTxtGp(ticks=gpar(cex=4), xlab=gpar(cex=4), cex = 4, summary = gpar(fontface = 'bold'),
                                              label = list(gpar(fontface = 'plain'),
                                                           gpar(fontface = 'plain', cex = 4),
                                                           gpar(fontface = 'plain', cex =4))))

dev.off()

# 05_Sensitivity analyses 1 and 2 ------

## SA 1: Matched on age_data_start, year_data_start and prac_region -------

load(file = "filepath") 
full_cons_glm_unadj <- glm_mig %>% filter(names=="Migrant") %>% mutate(outcome = 'Migrants vs. non-migrants')
full_cons_glm_unadj$names <- recode(full_cons_glm_unadj$names,
                                    "Migrant" = "Unadjusted")

load(file = "filepath") 
full_cons_glm_adj_noIMD <- multivariable_matched_allages %>% filter(names=="Migrant") %>% mutate(outcome = ' ')
full_cons_glm_adj_noIMD$names <- recode(full_cons_glm_adj_noIMD$names,
                                        "Migrant" = "Adjusted, without IMD")

load(file = "filepath")
full_cons_glm_adj_withIMD <- multivariable_matched_allages %>% filter(names=="Migrant") %>% mutate(outcome = '')
full_cons_glm_adj_withIMD$names <- recode(full_cons_glm_adj_withIMD$names,
                                          "Migrant" = "Adjusted, with IMD")

load(file = "filepath")
full_cons_glm_adj_withIMD_0to15 <- glm_multivariable_matched_0to15 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_0to15$names <- recode(full_cons_glm_adj_withIMD_0to15$names,
                                                "Migrant" = "    0-15 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_16to24 <- glm_multivariable_matched_16to24 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_16to24$names <- recode(full_cons_glm_adj_withIMD_16to24$names,
                                                 "Migrant" = "    16-24 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_25to34 <- glm_multivariable_matched_25to34 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_25to34$names <- recode(full_cons_glm_adj_withIMD_25to34$names,
                                                 "Migrant" = "    25-34 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_35to49 <- glm_multivariable_matched_35to49 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_35to49$names <- recode(full_cons_glm_adj_withIMD_35to49$names,
                                                 "Migrant" = "    35-49 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_50to64 <- glm_multivariable_matched_50to64 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_50to64$names <- recode(full_cons_glm_adj_withIMD_50to64$names,
                                                 "Migrant" = "    50-64 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_65plus <- glm_multivariable_matched_65plus %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_65plus$names <- recode(full_cons_glm_adj_withIMD_65plus$names, 
                                                 "Migrant" = "    65 years and older")

## Join all subgroups 

exact_matched_plot_data_SA1 <- bind_rows(full_cons_glm_unadj, full_cons_glm_adj_noIMD, full_cons_glm_adj_withIMD,
                                     full_cons_glm_adj_withIMD_0to15, full_cons_glm_adj_withIMD_16to24,
                                     full_cons_glm_adj_withIMD_25to34, full_cons_glm_adj_withIMD_35to49,
                                     full_cons_glm_adj_withIMD_50to64, full_cons_glm_adj_withIMD_65plus) %>%
  add_row(.before = 4, names = 'Adjusted, with IMD, by age group') 


## SA 2: Matched on pyears and prac_region ----------------

load(file = "filepath") 
full_cons_glm_unadj <- glm_mig %>% filter(names=="Migrant") %>% mutate(outcome = 'Migrants vs. non-migrants')
full_cons_glm_unadj$names <- recode(full_cons_glm_unadj$names,
                                    "Migrant" = "Unadjusted")

load(file = "filepath") 
full_cons_glm_adj_noIMD <- multivariable_matched_allages %>% filter(names=="Migrant") %>% mutate(outcome = ' ')
full_cons_glm_adj_noIMD$names <- recode(full_cons_glm_adj_noIMD$names,
                                        "Migrant" = "Adjusted, without IMD")

load(file = "filepath")
full_cons_glm_adj_withIMD <- multivariable_matched_allages %>% filter(names=="Migrant") %>% mutate(outcome = '')
full_cons_glm_adj_withIMD$names <- recode(full_cons_glm_adj_withIMD$names,
                                          "Migrant" = "Adjusted, with IMD")

load(file = "filepath")
full_cons_glm_adj_withIMD_0to15 <- glm_multivariable_matched_0to15 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_0to15$names <- recode(full_cons_glm_adj_withIMD_0to15$names,
                                                "Migrant" = "    0-15 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_16to24 <- glm_multivariable_matched_16to24 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_16to24$names <- recode(full_cons_glm_adj_withIMD_16to24$names,
                                                 "Migrant" = "    16-24 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_25to34 <- glm_multivariable_matched_25to34 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_25to34$names <- recode(full_cons_glm_adj_withIMD_25to34$names,
                                                 "Migrant" = "    25-34 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_35to49 <- glm_multivariable_matched_35to49 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_35to49$names <- recode(full_cons_glm_adj_withIMD_35to49$names,
                                                 "Migrant" = "    35-49 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_50to64 <- glm_multivariable_matched_50to64 %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_50to64$names <- recode(full_cons_glm_adj_withIMD_50to64$names,
                                                 "Migrant" = "    50-64 years")

load(file = "filepath")
full_cons_glm_adj_withIMD_65plus <- glm_multivariable_matched_65plus %>% filter(names=="Migrant") %>% mutate(outcome = '') %>% dplyr::select(-c(age_subcohort))
full_cons_glm_adj_withIMD_65plus$names <- recode(full_cons_glm_adj_withIMD_65plus$names, 
                                                 "Migrant" = "    65 years and older")

## Join all subgroups 

exact_matched_plot_data_SA2 <- bind_rows(full_cons_glm_unadj, full_cons_glm_adj_noIMD, full_cons_glm_adj_withIMD,
                                     full_cons_glm_adj_withIMD_0to15, full_cons_glm_adj_withIMD_16to24,
                                     full_cons_glm_adj_withIMD_25to34, full_cons_glm_adj_withIMD_35to49,
                                     full_cons_glm_adj_withIMD_50to64, full_cons_glm_adj_withIMD_65plus) %>%
  add_row(.before = 4, names = 'Adjusted, with IMD, by age group') 

## Combine SA 1 and 2 into wide format for forest plot -----

forest_plot_sa_all <- exact_matched_plot_data_SA1 %>%
  dplyr::select(-c(p, irr_ci))
forest_plot_sa_all <- forest_plot_sa_all %>%
  rename(estimateSA1 = estimate,
         lowerSA1 = lower,
         upperSA1 = upper,
         ciSA1 = ci)

forest_plot_sa_all$estimateSA2 <- exact_matched_plot_data_SA2$estimate
forest_plot_sa_all$lowerSA2 <- exact_matched_plot_data_SA2$lower
forest_plot_sa_all$upperSA2 <- exact_matched_plot_data_SA2$upper
forest_plot_sa_all$ciSA2 <- exact_matched_plot_data_SA2$ci

## SA1 and SA2 forest plot -----

sfrac <- function(top, bottom, data=NULL){
  with(data,lapply(paste0('atop(', top, ',', bottom, ')'), str2expression))}

tabletext <- list(
  c("Migrants vs. non-migrants", forest_plot_sa_all$names),
  c("RR", sfrac(forest_plot_sa_all$estimateSA1, forest_plot_sa_all$estimateSA2, data = forest_plot_sa_all)),
  c("95% CI", sfrac(forest_plot_sa_all$ciSA1, forest_plot_sa_all$ciSA2, data = forest_plot_sa_all)))

dev.new()
png(file = "filepath", width = 2500, height = 2000)
forest_plot_SA <- forestplot(tabletext, 
                             mean = cbind(c(NA, forest_plot_sa_all$estimateSA1), c(NA, forest_plot_sa_all$estimateSA2)),
                             lower = cbind(c(NA, forest_plot_sa_all$lowerSA1), c(NA, forest_plot_sa_all$lowerSA2)),
                             upper = cbind(c(NA, forest_plot_sa_all$upperSA1), c(NA, forest_plot_sa_all$upperSA2)),
                             new_page = TRUE,
                             zero = 1,
                             lineheight = unit(50, 'mm'),
                             line.margin = .1,
                             xlog = TRUE, 
                             xlab = 'Rate ratio with 95% CI',
                             col = fpColors(box = c('red4', 'skyblue3'),
                                            lines = c('red3', 'skyblue2')),
                             fn.ci.norm = c(fpDrawNormalCI, fpDrawDiamondCI),
                             boxsize = 0.1,
                             is.summary = c(TRUE, rep(FALSE, 10)),
                             xticks = c(0.5, 1.0, 1.5),
                             legend = c('Matched on age and year of joining CPRD Gold and practice region', 'Matched on follow-up time and practice region'),
                             ci.vertices = TRUE,
                             ci.vertices.height = 0.1, 
                             colgap = unit(10,'mm'),
                             graphwidth = unit(300, 'mm'),
                             graph.pos =2,
                             alim = c(0.5, 1.5),
                             align = c('l', 'c', 'c', 'c'),
                             title = ' ',
                             txt_gp = fpTxtGp(ticks=gpar(cex=4), xlab=gpar(cex=4), cex = 4, summary = gpar(fontface = 'bold'),
                                              label = list(gpar(fontface = 'plain'),
                                                           gpar(fontface = 'bold', cex = 4),
                                                           gpar(fontface = 'bold', cex =4))))
        
dev.off()

# 06_Sensitivity analysis: Migration certainty (full cohort) forest plot -----

load(file="filepath")
load(file="filepath")

# Definite migrants vs non-migrants
# Unadjusted
full_cons_migcertainty_unadj_definite <- glm_migcertainty %>% filter(names=="Definite") 
full_cons_migcertainty_unadj_definite$names <- recode(full_cons_migcertainty_unadj_definite$names,
                                                      "Definite" = "Unadjusted")
# Adjusted 
full_cons_migcertainty_adj_definite <- multivariable_migcertainty_allages %>% filter(names=="Definite") %>% mutate(names = 'Adjusted, with IMD')

# Probable migrants vs. non-migrants
# Unadjusted
full_cons_migcertainty_unadj_probable <- glm_migcertainty %>% filter(names=="Probable") 
full_cons_migcertainty_unadj_probable$names <- recode(full_cons_migcertainty_unadj_probable$names,
                                                      "Probable" = "Unadjusted")
# Adjusted 
full_cons_migcertainty_adj_probable <- multivariable_migcertainty_allages %>% filter(names=="Probable") %>% mutate(names = 'Adjusted, with IMD') 

migcertainty_plot_data_def <- bind_rows(full_cons_migcertainty_unadj_definite, full_cons_migcertainty_adj_definite)
                                    
migcertainty_plot_data_prob <-bind_rows(full_cons_migcertainty_unadj_probable, full_cons_migcertainty_adj_probable) 

# Combine definite and probable into wide format for forest plot

forest_plot_sa_all_migcert <- migcertainty_plot_data_def%>%
  dplyr::select(-c(p, irr_ci))
forest_plot_sa_all_migcert <- forest_plot_sa_all_migcert %>%
  rename(estimateDef = estimate,
         lowerDef = lower,
         upperDef = upper,
         ciDef = ci)

forest_plot_sa_all_migcert$estimateProb <- migcertainty_plot_data_prob$estimate
forest_plot_sa_all_migcert$lowerProb <- migcertainty_plot_data_prob$lower
forest_plot_sa_all_migcert$upperProb <- migcertainty_plot_data_prob$upper
forest_plot_sa_all_migcert$ciProb <- migcertainty_plot_data_prob$ci


# Make forest plot

tabletext <- list(
  c("Migrants vs. non-migrants", forest_plot_sa_all_migcert$names),
  c("RR", sfrac(forest_plot_sa_all_migcert$estimateDef, forest_plot_sa_all_migcert$estimateProb, data = forest_plot_sa_all_migcert)),
  c("95% CI", sfrac(forest_plot_sa_all_migcert$ciDef, forest_plot_sa_all_migcert$ciProb, data = forest_plot_sa_all_migcert)))

dev.new()
png(file = "filepath", width = 2000, height = 800)
forest_plot_SA <- forestplot(tabletext, 
                             mean = cbind(c(NA, forest_plot_sa_all_migcert$estimateDef), c(NA, forest_plot_sa_all_migcert$estimateProb)),
                             lower = cbind(c(NA, forest_plot_sa_all_migcert$lowerDef), c(NA, forest_plot_sa_all_migcert$lowerProb)),
                             upper = cbind(c(NA, forest_plot_sa_all_migcert$upperDef), c(NA, forest_plot_sa_all_migcert$upperProb)),
                             new_page = TRUE,
                             zero = 1,
                             lineheight = unit(60, 'mm'),
                             line.margin = .1,
                             xlog = TRUE, 
                             xlab = 'Rate ratio with 95% CI',
                             col = fpColors(box = c('red4', 'skyblue3'),
                                            lines = c('red3', 'skyblue2')),
                             fn.ci.norm = c(fpDrawNormalCI, fpDrawDiamondCI),
                             boxsize = 0.1,
                             is.summary = c(TRUE, rep(FALSE, )),
                             xticks = c(0.5, 1.0, 1.5),
                             legend = c('Definite migrants', 'Probable migrants'),
                             ci.vertices = TRUE,
                             ci.vertices.height = 0.1, 
                             colgap = unit(10,'mm'),
                             graphwidth = unit(300, 'mm'),
                             graph.pos =2,
                             alim = c(0.5, 1.5),
                             align = c('l', 'c', 'c', 'c'),
                             title = ' ',
                             txt_gp = fpTxtGp(ticks=gpar(cex=4), xlab=gpar(cex=4), cex = 4, summary = gpar(fontface = 'bold'),
                                              label = list(gpar(fontface = 'plain'),
                                                           gpar(fontface = 'bold', cex = 4),
                                                           gpar(fontface = 'bold', cex =4))))

dev.off()
