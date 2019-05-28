###########################################################
### Imperial HBV model:                                 ###
### Clean input natural history data for The Gambia     ###
### Source: mapping review                              ###
###########################################################
## Load packages and set directories ----
require(tidyr)  # for data processing
require(dplyr)  # for data processing
require(here)  # for setting working directory
require(ggplot2)

inpath_hbvdata <- "data_raw"
outpath_hbvdata <- "data_clean"

## Load datasets ----

# Gambia age- and sex-specific HBsAg and anti-HBc prevalence in the general population
input_hbsag_antihbc_prev <- read.csv(here(inpath_hbvdata,
                                          "hbsag_prevalence.csv"),
                                     header = TRUE, check.names = FALSE,
                                     stringsAsFactors = FALSE)

# Gambia age- and sex-specific HBeAg prevalence in HBV carriers
input_hbeag_prev <- read.csv(here(inpath_hbvdata,
                                  "hbeag_prevalence.csv"),
                             header = TRUE, check.names = FALSE,
                             stringsAsFactors = FALSE)

# Gambia other (natural history) prevalence data in HBV carriers
input_natural_history_prev <- read.csv(here(inpath_hbvdata,
                                            "natural_history_prevalence.csv"),
                                       header = TRUE, check.names = FALSE,
                                       stringsAsFactors = FALSE)

# Natural history progression rates in West Africa
input_progression_rates <- read.csv(here(inpath_hbvdata,
                                         "natural_history_progression_rates.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

# Mortality curves in West Africa (time-to-event analysis)
input_mortality_curves <- read.csv(here(inpath_hbvdata,
                                         "natural_history_survival_curves.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)


# Mother-to-child transmission risk in West Africa
input_mtct_risk <- read.csv(here(inpath_hbvdata,
                                  "mtct_risk.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

# Risk of chronic carriage in West Africa
input_chronic_carriage_risk <- read.csv(here(inpath_hbvdata,
                                 "chronic_carriage_risk.csv"),
                            header = TRUE, check.names = FALSE,
                            stringsAsFactors = FALSE)


# Gambia age- and sex-specific HBeAg prevalence in cirrhosis and HCC patients
input_hbeag_prev_ld_patients <- read.csv(here(inpath_hbvdata,
                                  "hbeag_prevalence_ld_patients.csv"),
                             header = TRUE, check.names = FALSE,
                             stringsAsFactors = FALSE)


# Odds ratios
input_odds_ratios <- read.csv(here(inpath_hbvdata,
                                   "odds_ratios.csv"),
                                         header = TRUE, check.names = FALSE,
                                         stringsAsFactors = FALSE)

# Infant vaccine efficacy against chronic carriage
input_infant_vaccine_efficacy <- read.csv(here(inpath_hbvdata,
                                          "infant_vaccine_efficacy.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE)

# PAF for association of HBsAg and HCC/cirrhosis
input_paf_liver_disease <- read.csv(here(inpath_hbvdata,
                                         "paf_liver_disease.csv"),
                                          header = TRUE, check.names = FALSE,
                                          stringsAsFactors = FALSE)

# Mean age and proportion male in HBV-related liver disease patient samples
input_liver_disease_demography <- read.csv(here(inpath_hbvdata,
                                         "liver_disease_demography.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)


#### CLEAN CALIBRATION DATASETS

## HBsAg and anti-HBc prevalence in The Gambia dataset ----
subset_hbsag_prev <- select(input_hbsag_antihbc_prev,
                            id_paper,
                            id_group,
                            id_proc,
                            pop_group_clinical,
                            pop_group_demographic,
                            geographic_scope,
                            location,
                            recruitment_setting,
                            study_link,
                            dp_period,
                            starts_with("age"),
                            sex,
                            proportion_male,
                            vaccinated,
                            hbsag_positive_prop,
                            hbsag_positive_prop_ci_lower,
                            hbsag_positive_prop_ci_upper,
                            antihbc_positive_prop,
                            antihbc_positive_prop_ci_lower,
                            antihbc_positive_prop_ci_upper,
                            sample_size) 

## Processing

# 1) Assign a specific age to each data point
# Use mean age if available
subset_hbsag_prev$age_assign_years <- subset_hbsag_prev$age_mean_years
# Else use median age (may be estimated from frequency distribution)
subset_hbsag_prev$age_assign_years[subset_hbsag_prev$age_assign_years == "NR"] <-
  subset_hbsag_prev$age_median_years[subset_hbsag_prev$age_assign_years == "NR"]
# Else use mid-point of age range
subset_hbsag_prev$age_assign_years[subset_hbsag_prev$age_assign_years == "NR"] <-
  (as.numeric(subset_hbsag_prev$age_min_years[subset_hbsag_prev$age_assign_years == "NR"]) +
  as.numeric(subset_hbsag_prev$age_max_years[subset_hbsag_prev$age_assign_years == "NR"]) + 1)/2

# 2) Assign a specific year to each data point
# Split datapoint collection period column into minimum and maximum year
subset_hbsag_prev <- subset_hbsag_prev %>%
  separate(col = dp_period, into = c("dp_period_min", "dp_period_max"),
           sep = "-", remove = FALSE) %>%
  mutate(dp_period_max = coalesce(dp_period_max, dp_period_min)) 
# Take mean of minimum and maximum time of data collection
subset_hbsag_prev$dp_period_assign_years <- (as.numeric(subset_hbsag_prev$dp_period_min) +
                                             as.numeric(subset_hbsag_prev$dp_period_max))/2


# In Keneba and Manduar, vaccination was introduced in 1985 rather than in 1991 like
# the rest of The Gambia. For use in model, therefore need to assign an appropriate 
# year of data collection with regards to vaccination introduction
subset_hbsag_prev[(subset_hbsag_prev$location == "Keneba" | 
                     subset_hbsag_prev$location == "Manduar") &
                    subset_hbsag_prev$dp_period_assign_years >= 1985,] <- filter(subset_hbsag_prev, 
                             (location == "Keneba" | location == "Manduar") &
                               dp_period_assign_years >= 1985) %>%
  mutate(dp_period_assign_years = replace(dp_period_assign_years, values = (1991+(dp_period_assign_years-1985))))

# 3) Assign to the pre- or post-vaccination period (1991)
# Assign post-vaccination status to data points collected in 1991 or after
# and pre-vaccination status to data points collected before 1991
subset_hbsag_prev$dp_period_vacc <- "post-vacc"
subset_hbsag_prev$dp_period_vacc[subset_hbsag_prev$dp_period_assign_years < 1991] <- "pre-vacc"

# 4) Assign data points to a series based on the study link (if available) or the IDs
subset_hbsag_prev$series <- subset_hbsag_prev$study_link
subset_hbsag_prev$series[is.na(subset_hbsag_prev$series) == TRUE] <-
  paste0(subset_hbsag_prev$id_paper[is.na(subset_hbsag_prev$series) == TRUE], "-",
        subset_hbsag_prev$id_group[is.na(subset_hbsag_prev$series) == TRUE], "-",
        subset_hbsag_prev$id_proc[is.na(subset_hbsag_prev$series) == TRUE])
# For Keneba and Manduar, separate by village
subset_hbsag_prev$series[subset_hbsag_prev$location == "Keneba"] <- "Keneba"
subset_hbsag_prev$series[subset_hbsag_prev$location == "Manduar"] <- "Manduar"
# For GMB12 (Ryder) combine group IDs 3 and 4 (these are all family contacts with different age/sex)
subset_hbsag_prev$series[subset_hbsag_prev$series == "GMB12-3-x" |
                         subset_hbsag_prev$series == "GMB12-4-a" |
                         subset_hbsag_prev$series == "GMB12-4-b"] <- "GMB12-3+4"

# Save a subset of HBsAg prevalence dataset for use in model
hbsag_dataset_for_fitting <- select(subset_hbsag_prev,
                                    id_paper,
                                    id_group,
                                    id_proc,
                                    pop_group_clinical,
                                    dp_period_assign_years,
                                    sex,
                                    age_assign_years,
                                    hbsag_positive_prop, 
                                    hbsag_positive_prop_ci_lower,
                                    hbsag_positive_prop_ci_upper,
                                    sample_size) %>%
  filter(hbsag_positive_prop != "DUPLICATE" & hbsag_positive_prop != "NR")

# Turn number columns into numeric format
hbsag_dataset_for_fitting[,c(5,7:11)] <- apply(hbsag_dataset_for_fitting[,c(5,7:11)], 2, 
                                               function(x) as.numeric(x))

hbsag_dataset_for_fitting <- cbind(outcome = "HBsAg_prevalence",
                                   hbsag_dataset_for_fitting) %>%
  arrange(sex, dp_period_assign_years, id_paper, id_group, age_assign_years)
# Replace spaces and dashes with underscores
hbsag_dataset_for_fitting$pop_group_clinical <- gsub(" - |-|[[:space:]]", "_", 
                                                     hbsag_dataset_for_fitting$pop_group_clinical)

# Round ages down to the nearest 0.5 to match model output
hbsag_dataset_for_fitting$age_assign_years <- floor(hbsag_dataset_for_fitting$age_assign_years/0.5)*0.5
# Rename columns to match model output
names(hbsag_dataset_for_fitting)[names(hbsag_dataset_for_fitting)=="dp_period_assign_years"] <- "time"
names(hbsag_dataset_for_fitting)[names(hbsag_dataset_for_fitting)=="age_assign_years"] <- "age"
names(hbsag_dataset_for_fitting)[names(hbsag_dataset_for_fitting)=="hbsag_positive_prop"] <- "data_value"
names(hbsag_dataset_for_fitting)[names(hbsag_dataset_for_fitting)=="hbsag_positive_prop_ci_lower"] <- "ci_lower"
names(hbsag_dataset_for_fitting)[names(hbsag_dataset_for_fitting)=="hbsag_positive_prop_ci_upper"] <- "ci_upper"

#write.csv(hbsag_dataset_for_fitting, file = here(outpath_hbvdata, "hbsag_prevalence.csv"), row.names = FALSE)

# Save a subset of anti-HBc prevalence dataset for use in model
# We can only fit anti-HBc data as carriers+immunes until vaccination is introduced
# Therefore need to remove datapoints after 1990 and after 1984 in Keneba or Manduar
antihbc_dataset_for_fitting <- select(subset_hbsag_prev,
                                    id_paper,
                                    id_group,
                                    id_proc,
                                    pop_group_clinical,
                                    location,
                                    dp_period_assign_years,
                                    sex,
                                    age_assign_years,
                                    antihbc_positive_prop, 
                                    antihbc_positive_prop_ci_lower,
                                    antihbc_positive_prop_ci_upper,
                                    sample_size) %>%
  filter(antihbc_positive_prop != "DUPLICATE" & antihbc_positive_prop != "NR") %>%
  filter(dp_period_assign_years < 1991 | 
         (location == "Keneba" & dp_period_assign_years < 1985) |
         (location == "Manduar" & dp_period_assign_years < 1985)) %>%
  select(-location)
           

# Turn number columns into numeric format
antihbc_dataset_for_fitting[,c(5,7:11)] <- apply(antihbc_dataset_for_fitting[,c(5,7:11)], 2, 
                                               function(x) as.numeric(x))

antihbc_dataset_for_fitting <- cbind(outcome = "Anti_HBc_prevalence",
                                   antihbc_dataset_for_fitting) %>%
  arrange(sex, dp_period_assign_years, id_paper, id_group, age_assign_years)
# Replace spaces and dashes with underscores
antihbc_dataset_for_fitting$pop_group_clinical <- gsub(" - |-|[[:space:]]", "_", 
                                                     antihbc_dataset_for_fitting$pop_group_clinical)

# Round ages down to the nearest 0.5 to match model output
antihbc_dataset_for_fitting$age_assign_years <- floor(antihbc_dataset_for_fitting$age_assign_years/0.5)*0.5
# Rename columns to match model output
names(antihbc_dataset_for_fitting)[names(antihbc_dataset_for_fitting)=="dp_period_assign_years"] <- "time"
names(antihbc_dataset_for_fitting)[names(antihbc_dataset_for_fitting)=="age_assign_years"] <- "age"
names(antihbc_dataset_for_fitting)[names(antihbc_dataset_for_fitting)=="antihbc_positive_prop"] <- "data_value"
names(antihbc_dataset_for_fitting)[names(antihbc_dataset_for_fitting)=="antihbc_positive_prop_ci_lower"] <- "ci_lower"
names(antihbc_dataset_for_fitting)[names(antihbc_dataset_for_fitting)=="antihbc_positive_prop_ci_upper"] <- "ci_upper"

#write.csv(antihbc_dataset_for_fitting, file = here(outpath_hbvdata, "antihbc_prevalence.csv"), row.names = FALSE)


## HBsAg prevalence plots ----
## PRE-VACCINATION PLOTS
# Plot all data points
ggplot() +
  geom_point(data = filter(subset_hbsag_prev, dp_period_vacc == "pre-vacc"),
             aes(x = as.numeric(age_assign_years), y = as.numeric(hbsag_positive_prop)), size = 3) +
  theme_bw() + ylim(0, 0.4) + xlim(0,80) + ggtitle("Pre-1991 data points") +
  ylab("HBsAg prevalence (proportion)") + xlab("Age (years)")

# Plot data points in Keneba and Manduar
ggplot(data = filter(subset_hbsag_prev, dp_period_vacc == "pre-vacc",
                     study_link == "KM vaccine cohort"),
       aes(x = as.numeric(age_assign_years), y = as.numeric(hbsag_positive_prop),
           color = series)) +
  geom_text(aes(label = dp_period), vjust=1.5) +
  geom_point(size = 3) +
  theme_bw() + ylim(0, 0.4) + xlim(0,80) + ggtitle("Pre-1991 data points") +
  ylab("HBsAg prevalence (proportion)") + xlab("Age (years)")
# Note at data points from 1989, vaccination had already been introduced but not in the
# included age groups here (not in those >9 years old)

# Plot points for all other studies
ggplot() +
  geom_point(data = filter(subset_hbsag_prev, dp_period_vacc == "pre-vacc",
                           series != "Keneba", series != "Manduar"),
             aes(x = as.numeric(age_assign_years), y = as.numeric(hbsag_positive_prop),
                 color = series), size = 3) +
  theme_bw() + ylim(0, 0.4) + xlim(0,80) + ggtitle("Pre-1991 data points") +
  ylab("HBsAg prevalence (proportion)") + xlab("Age (years)") +
  scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA",
                              "#3F4921", "#D1A33D", "#8569D5", "#5F7FC7",
                              "#673770", "#38333E"))


## POST-VACCINATION PLOTS
# Plot all data points
ggplot(data = filter(subset_hbsag_prev, dp_period_vacc == "post-vacc"),
       aes(x = as.numeric(age_assign_years), y = as.numeric(hbsag_positive_prop),
           color = series)) +
  geom_point(size = 3) +
  geom_text(aes(label = dp_period_assign_years), vjust=1.5) +
  theme_bw() + ylim(0, 0.4) + xlim(0,80) + ggtitle("Post-1991 data points") +
  ylab("HBsAg prevalence (proportion)") + xlab("Age (years)")  +
  scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA",
                              "#3F4921", "#D1A33D", "#8569D5", "#5F7FC7",
                              "#673770", "#38333E"))


## HBeAg prevalence in The Gambia dataset ----
# This is only for HBeAg prevalence in carriers overall
subset_hbeag_prev <- select(input_hbeag_prev,
                            id_paper,
                            id_group,
                            id_proc,
                            pop_group_clinical,
                            pop_group_demographic,
                            geographic_scope,
                            location,
                            recruitment_setting,
                            study_link,
                            dp_period,
                            starts_with("age"),
                            sex,
                            proportion_male,
                            vaccinated,
                            hbeag_positive_prop,
                            hbeag_positive_prop_ci_lower,
                            hbeag_positive_prop_ci_upper,
                            sample_size)

## Processing

# 1) Assign a specific age to each data point
# Use mean age if available
subset_hbeag_prev$age_assign_years <- subset_hbeag_prev$age_mean_years
# Else use median age (may be estimated from frequency distribution)
subset_hbeag_prev$age_assign_years[subset_hbeag_prev$age_assign_years == "NR"] <-
  subset_hbeag_prev$age_median_years[subset_hbeag_prev$age_assign_years == "NR"]
# Else use mid-point of age range
subset_hbeag_prev$age_assign_years[subset_hbeag_prev$age_assign_years == "NR"] <-
  (as.numeric(subset_hbeag_prev$age_min_years[subset_hbeag_prev$age_assign_years == "NR"]) +
     as.numeric(subset_hbeag_prev$age_max_years[subset_hbeag_prev$age_assign_years == "NR"]) + 1)/2

# 2) Assign a specific year to each data point
# Split datapoint collection period column into minimum and maximum year
subset_hbeag_prev <- subset_hbeag_prev %>%
  separate(col = dp_period, into = c("dp_period_min", "dp_period_max"),
           sep = "-", remove = FALSE) %>%
  mutate(dp_period_max = coalesce(dp_period_max, dp_period_min)) 
# Take mean of minimum and maximum time of data collection
subset_hbeag_prev$dp_period_assign_years <- (as.numeric(subset_hbeag_prev$dp_period_min) +
                                               as.numeric(subset_hbeag_prev$dp_period_max))/2

# Save a subset of dataset for use in model
hbeag_dataset_for_fitting <- select(subset_hbeag_prev,
                                    id_paper,
                                    id_group,
                                    id_proc,
                                    pop_group_clinical,
                                    dp_period_assign_years,
                                    sex,
                                    age_assign_years,
                                    hbeag_positive_prop, 
                                    hbeag_positive_prop_ci_lower,
                                    hbeag_positive_prop_ci_upper,
                                    sample_size)

# Turn number columns into numeric format
hbeag_dataset_for_fitting[,c(5,7:11)] <- apply(hbeag_dataset_for_fitting[,c(5,7:11)], 2, 
                                               function(x) as.numeric(x))

hbeag_dataset_for_fitting <- cbind(outcome = "HBeAg_prevalence",
                                   hbeag_dataset_for_fitting) %>%
  arrange(sex, dp_period_assign_years, id_paper, id_group, age_assign_years)
# Replace spaces and dashes with underscores
hbeag_dataset_for_fitting$pop_group_clinical <- gsub(" - |-|[[:space:]]", "_", 
                                                     hbeag_dataset_for_fitting$pop_group_clinical)

# Round ages down to the nearest 0.5 to match model output
hbeag_dataset_for_fitting$age_assign_years <- floor(hbeag_dataset_for_fitting$age_assign_years/0.5)*0.5
# Rename columns to match model output
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="dp_period_assign_years"] <- "time"
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="age_assign_years"] <- "age"
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="hbeag_positive_prop"] <- "data_value"
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="hbeag_positive_prop_ci_lower"] <- "ci_lower"
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="hbeag_positive_prop_ci_upper"] <- "ci_upper"

#write.csv(hbeag_dataset_for_fitting, file = here(outpath_hbvdata, "hbeag_prevalence.csv"), row.names = FALSE)



# Natural history prevalence dataset ----
# Here the age and timepoint assignment is different from sAg/anti-HBc/eAg datasets!
subset_natural_history_prev <- select(input_natural_history_prev,
                            id_paper,
                            id_group,
                            id_proc,
                            pop_group_clinical,
                            pop_group_demographic,
                            geographic_scope,
                            location,
                            recruitment_setting,
                            study_link,
                            recruitment_period,
                            dp_period,
                            starts_with("age"),
                            sex,
                            proportion_male,
                            vaccinated,
                            dp_description,
                            model_numerator,
                            model_denominator,
                            prevalence,
                            prevalence_ci_lower,
                            prevalence_ci_upper,
                            sample_size) 

## Processing

# 1) Assign an age range to each datapoint
# Use given age range if available
subset_natural_history_prev$age_assign_min_years <- subset_natural_history_prev$age_min_years
subset_natural_history_prev$age_assign_max_years <- subset_natural_history_prev$age_max_years
# Else use interquartile range
# Splitting IQR column into the minimum and maximum value
subset_natural_history_prev$age_assign_min_years[subset_natural_history_prev$age_assign_min_years == "NR"] <-
  unlist(lapply(strsplit(subset_natural_history_prev$age_iqr_years[subset_natural_history_prev$age_assign_min_years == "NR"], "-"), "[[",1))
subset_natural_history_prev$age_assign_max_years[subset_natural_history_prev$age_assign_max_years == "NR"] <-
sapply(strsplit(subset_natural_history_prev$age_iqr_years[subset_natural_history_prev$age_assign_max_years == "NR"], "-"), "[",2)
# Else use mean +/- 2*SD (in a normal distribution, 95% of data falls within this range)
subset_natural_history_prev$age_assign_min_years[is.na(subset_natural_history_prev$age_assign_min_years) == TRUE] <-
  (as.numeric(subset_natural_history_prev$age_mean_years[is.na(subset_natural_history_prev$age_assign_min_years) == TRUE])-2*
  subset_natural_history_prev$age_sd_years[is.na(subset_natural_history_prev$age_assign_min_years) == TRUE])
subset_natural_history_prev$age_assign_max_years[is.na(subset_natural_history_prev$age_assign_max_years) == TRUE] <-
  (as.numeric(subset_natural_history_prev$age_mean_years[is.na(subset_natural_history_prev$age_assign_max_years) == TRUE])+2*
     subset_natural_history_prev$age_sd_years[is.na(subset_natural_history_prev$age_assign_max_years) == TRUE])


# 2) Assign a specific year to each data point WITH EXCEPTION
# Split datapoint collection period column into minimum and maximum year
subset_natural_history_prev <- subset_natural_history_prev %>%
  separate(col = dp_period, into = c("dp_period_min", "dp_period_max"),
           sep = "-", remove = FALSE) %>%
  mutate(dp_period_max = coalesce(dp_period_max, dp_period_min)) 
# Take mean of minimum and maximum time of data collection
# Round to nearest integer for simplicity
subset_natural_history_prev$dp_period_assign_years <- round((as.numeric(subset_natural_history_prev$dp_period_min) +
                                               as.numeric(subset_natural_history_prev$dp_period_max))/2,0)
# EXCEPTION: Shimakawa baseline data depends on the recruitment period (1974-2008)
# We know the follow-up lasted a median of 27 years and the later timepoint was 2013
# Therefore, overwrite with assumed baseline dp_period of 1986 (2013-27)
subset_natural_history_prev$dp_period_assign_years[subset_natural_history_prev$id_paper == "1" &
                                                     subset_natural_history_prev$dp_period_min == "1974"] <- 1986


# Save a subset of natural history prevalence dataset for use in model
# NOTE: filtering out datapoints that are not directly applicable prevalence measurements
# Need to come back to those later
natural_history_prev_for_fitting <- select(subset_natural_history_prev,
                                    id_paper,
                                    id_group,
                                    id_proc,
                                    model_numerator,
                                    model_denominator,
                                    dp_period_assign_years,
                                    sex,
                                    age_assign_min_years,
                                    age_assign_max_years,
                                    prevalence,
                                    prevalence_ci_lower,
                                    prevalence_ci_upper,
                                    sample_size) 

# Turn number columns into numeric format
natural_history_prev_for_fitting[,8:13] <- apply(natural_history_prev_for_fitting[,8:13], 2, 
                                                 function(x) as.numeric(x))

#natural_history_prev_for_fitting$model_denominator[
#  natural_history_prev_for_fitting$id_paper == "GMB2"] <- "HCC"
#natural_history_prev_for_fitting$model_numerator[
#  natural_history_prev_for_fitting$id_paper == "GMB2"][1] <- "CC"
#natural_history_prev_for_fitting$model_numerator[
#  natural_history_prev_for_fitting$id_paper == "GMB2"][2] <- "DCC"

# Define prevalence outcome for each row
prevalence_outcome <- paste0(natural_history_prev_for_fitting$model_numerator, "_prevalence_in_", natural_history_prev_for_fitting$model_denominator)


natural_history_prev_for_fitting <- cbind(prevalence_outcome = prevalence_outcome,
                                          natural_history_prev_for_fitting) %>%
  arrange(sex, dp_period_assign_years, id_paper, id_group, age_assign_min_years) %>%
  select(-model_numerator, -model_denominator)

# Abbreviate model numerator and denominators
natural_history_prev_for_fitting$prevalence_outcome <- gsub("Immune tolerant", "IT", natural_history_prev_for_fitting$prevalence_outcome, ignore.case = TRUE)
natural_history_prev_for_fitting$prevalence_outcome <- gsub("Immune reactive", "IR", natural_history_prev_for_fitting$prevalence_outcome, ignore.case = TRUE)
natural_history_prev_for_fitting$prevalence_outcome <- gsub("Inactive carrier", "IC", natural_history_prev_for_fitting$prevalence_outcome, ignore.case = TRUE)
natural_history_prev_for_fitting$prevalence_outcome <- gsub("Decompensated cirrhosis", "DCC", natural_history_prev_for_fitting$prevalence_outcome, ignore.case = TRUE)
natural_history_prev_for_fitting$prevalence_outcome <- gsub("Compensated cirrhosis", "CC", natural_history_prev_for_fitting$prevalence_outcome, ignore.case = TRUE)

natural_history_prev_for_fitting$prevalence_outcome <- tolower(natural_history_prev_for_fitting$prevalence_outcome)

# Replace spaces with underscores
natural_history_prev_for_fitting$prevalence_outcome <- gsub("[[:space:]]", "_", 
                                                         natural_history_prev_for_fitting$prevalence_outcome)

# Add output IDs based on numerator to ensure correct mapping:
# Only separate compartment names by underscores 
natural_history_prev_for_fitting$id_output <- gsub("_prevalence.*", "", natural_history_prev_for_fitting$prevalence_outcome)
natural_history_prev_for_fitting$id_output <- gsub("and_", "", natural_history_prev_for_fitting$id_output )
natural_history_prev_for_fitting$id_output <- gsub(",", "", natural_history_prev_for_fitting$id_output)
natural_history_prev_for_fitting$id_output <- gsub("originating_", "", natural_history_prev_for_fitting$id_output)
natural_history_prev_for_fitting$id_output[natural_history_prev_for_fitting$id_paper == "A4"] <- "shadow_incident_deaths"

# Create a unique identifier for each row from paper ID, group ID, timepoint and output ID
# Convert to lowercase
natural_history_prev_for_fitting$id_unique <- paste0("id_",
                                                     natural_history_prev_for_fitting$id_paper, "_",
                                                     natural_history_prev_for_fitting$id_group, "_",
                                                     natural_history_prev_for_fitting$dp_period_assign_years, "_",
                                                     natural_history_prev_for_fitting$id_output)
natural_history_prev_for_fitting$id_unique <- tolower(natural_history_prev_for_fitting$id_unique)
natural_history_prev_for_fitting$id_output <- NULL

# Round ages down to the nearest 0.5 to match model output
natural_history_prev_for_fitting$age_assign_min_years <- floor(natural_history_prev_for_fitting$age_assign_min_years/0.5)*0.5
natural_history_prev_for_fitting$age_assign_max_years <- floor(natural_history_prev_for_fitting$age_assign_max_years/0.5)*0.5

# Rename columns to match model output
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="prevalence_outcome"] <- "outcome"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="dp_period_assign_years"] <- "time"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="age_assign_min_years"] <- "age_min"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="age_assign_max_years"] <- "age_max"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="prevalence"] <- "data_value"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="prevalence_ci_lower"] <- "ci_lower"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="prevalence_ci_upper"] <- "ci_upper"

# HBeAg prevalence in Gambian liver disease patients ----
subset_hbeag_prev_ld_patients <- select(input_hbeag_prev_ld_patients,
                                        id_paper,
                                        id_group,
                                        id_proc,
                                        pop_group_clinical,
                                        pop_group_demographic,
                                        geographic_scope,
                                        location,
                                        recruitment_setting,
                                        study_link,
                                        dp_period,
                                        starts_with("age"),
                                        sex,
                                        proportion_male,
                                        vaccinated,
                                        hbeag_positive_prop,
                                        hbeag_positive_prop_ci_lower,
                                        hbeag_positive_prop_ci_upper,
                                        sample_size)


# 2) Assign a specific year to each data point
# Split datapoint collection period column into minimum and maximum year
subset_hbeag_prev_ld_patients <- subset_hbeag_prev_ld_patients%>%
  separate(col = dp_period, into = c("dp_period_min", "dp_period_max"),
           sep = "-", remove = FALSE) %>%
  mutate(dp_period_max = coalesce(dp_period_max, dp_period_min)) 
# Take mean of minimum and maximum time of data collection
subset_hbeag_prev_ld_patients$dp_period_assign_years <- (as.numeric(subset_hbeag_prev_ld_patients$dp_period_min) +
                                                           as.numeric(subset_hbeag_prev_ld_patients$dp_period_max))/2



# Save a subset of dataset for use in model
hbeag_ld_for_fitting <- select(subset_hbeag_prev_ld_patients,
                               id_paper,
                               id_group,
                               id_proc,
                               pop_group_clinical,
                               dp_period_assign_years,
                               sex,
                               age_min_years,
                               age_max_years,
                               hbeag_positive_prop, 
                               hbeag_positive_prop_ci_lower,
                               hbeag_positive_prop_ci_upper,
                               sample_size)



# Turn number columns into numeric format
hbeag_ld_for_fitting[,c(5,7:12)] <- apply(hbeag_ld_for_fitting[,c(5,7:12)], 2, 
                                          function(x) as.numeric(x))

# Assign outcome and remove pop group
hbeag_ld_outcome <- c(rep("hbeag_prevalence_in_hcc",7),
                      rep("hbeag_prevalence_in_cirrhosis",4))

hbeag_ld_for_fitting <- cbind(outcome = hbeag_ld_outcome,
                              hbeag_ld_for_fitting) %>%
  arrange(sex, dp_period_assign_years, id_paper, id_group, age_min_years) %>%
  select(-pop_group_clinical)


# Add output IDs based on outcome to ensure correct mapping:
# Only separate compartment names by underscores 
hbeag_ld_for_fitting$id_output <- gsub("_prevalence_in", "", hbeag_ld_for_fitting$outcome)

# Create a unique identifier for each row from paper ID, group ID, timepoint and output ID
# Convert to lowercase
hbeag_ld_for_fitting$id_unique <- paste0("id_",
                                         hbeag_ld_for_fitting$id_paper, "_",
                                         hbeag_ld_for_fitting$id_group, "_",
                                         hbeag_ld_for_fitting$dp_period_assign_years, "_",
                                         hbeag_ld_for_fitting$id_output)
hbeag_ld_for_fitting$id_unique <- tolower(hbeag_ld_for_fitting$id_unique)
hbeag_ld_for_fitting$id_output <- NULL

# Rename columns to match model output
names(hbeag_ld_for_fitting)[names(hbeag_ld_for_fitting)=="dp_period_assign_years"] <- "time"
names(hbeag_ld_for_fitting)[names(hbeag_ld_for_fitting)=="age_min_years"] <- "age_min"
names(hbeag_ld_for_fitting)[names(hbeag_ld_for_fitting)=="age_max_years"] <- "age_max"
names(hbeag_ld_for_fitting)[names(hbeag_ld_for_fitting)=="hbeag_positive_prop"] <- "data_value"
names(hbeag_ld_for_fitting)[names(hbeag_ld_for_fitting)=="hbeag_positive_prop_ci_lower"] <- "ci_lower"
names(hbeag_ld_for_fitting)[names(hbeag_ld_for_fitting)=="hbeag_positive_prop_ci_upper"] <- "ci_upper"



## Append HBeAg prevalence in LD patients to natural history prevalence dataset ----
natural_history_prev_for_fitting <- rbind(natural_history_prev_for_fitting,
                                          hbeag_ld_for_fitting)

#write.csv(natural_history_prev_for_fitting, file = here(outpath_hbvdata, "natural_history_prevalence.csv"), row.names = FALSE)

## Natural history progression rates in West Africa ----
subset_progression_rates <- input_progression_rates %>%
                            select(id_paper,
                            id_group,
                            id_proc,
                            pop_group_clinical,
                            pop_group_demographic,
                            recruitment_period,
                            dp_period,
                            country,
                            study_link,
                            study_details,
                            starts_with("bl_age"),
                            starts_with("current_age"),
                            sex,
                            proportion_male,
                            vaccinated,
                            model_prog_from,
                            model_prog_to,
                            rate_100py,
                            rate_100py_ci_lower,
                            rate_100py_ci_upper,
                            py_at_risk,
                            sample_size,
                            starts_with("follow_up"),
                            dp_details,
                            modelling_notes) %>%
  mutate(rate_py = rate_100py/100,  # Convert rate per 100 person-years to per person-year
         rate_py_ci_lower = as.numeric(rate_100py_ci_lower)/100,
         rate_py_ci_upper = as.numeric(rate_100py_ci_upper)/100)

# Split into 2 datasets based on use as input or output within model
prog_rates_for_input <- subset(subset_progression_rates, grepl("Input.*", subset_progression_rates$modelling_notes))
prog_rates_for_output <- subset(subset_progression_rates, !grepl("Input.*", subset_progression_rates$modelling_notes))

# Output dataset (to fit to)
# 1) Assign a specific age to each data point
# Use mean age if available
prog_rates_for_output$age_assign_years <- prog_rates_for_output$bl_age_mean_years
# Else use median age (may be estimated from frequency distribution)
prog_rates_for_output$age_assign_years[prog_rates_for_output$age_assign_years == "NR"] <-
  prog_rates_for_output$bl_age_median_years[prog_rates_for_output$age_assign_years == "NR"]
# Else use mid-point of age range
prog_rates_for_output$age_assign_years[prog_rates_for_output$age_assign_years == "NR"] <-
  (as.numeric(prog_rates_for_output$bl_age_min_years[prog_rates_for_output$age_assign_years == "NR"]) +
     as.numeric(prog_rates_for_output$bl_age_max_years[prog_rates_for_output$age_assign_years == "NR"]) + 1)/2
# Round ages
prog_rates_for_output$age_assign_years <- round(as.numeric(prog_rates_for_output$age_assign_years))


# 2) Assign a specific follow-up time to each data point
# Use mean age if available
prog_rates_for_output$fu_assign_years <- prog_rates_for_output$follow_up_mean_years
# Else use median age (may be estimated from frequency distribution)
prog_rates_for_output$fu_assign_years[prog_rates_for_output$fu_assign_years == "NR"] <-
  prog_rates_for_output$follow_up_median_years[prog_rates_for_output$fu_assign_years == "NR"]
# Else use mid-point of age range
prog_rates_for_output$fu_assign_years[prog_rates_for_output$fu_assign_years == "NR"] <-
  (as.numeric(prog_rates_for_output$follow_up_min_years[prog_rates_for_output$fu_assign_years == "NR"]) +
     as.numeric(prog_rates_for_output$follow_up_max_years[prog_rates_for_output$fu_assign_years == "NR"]) + 1)/2
# One study only has maximum follow-up time, use this instead
prog_rates_for_output$fu_assign_years[is.na(prog_rates_for_output$fu_assign_years == TRUE)] <-
  as.numeric(prog_rates_for_output$follow_up_max_years[is.na(prog_rates_for_output$fu_assign_years == TRUE)])
# Round follow-up times
prog_rates_for_output$fu_assign_years <- round(as.numeric(prog_rates_for_output$fu_assign_years))

# 3) Assign a specific start time for the cohort (first year of recruitment)
prog_rates_for_output$start_period_assign_years <- substr(prog_rates_for_output$recruitment_period,1,4)
# EXCEPTION: for Shimakawa cohort, assign 1985 as start period (2013-28)
prog_rates_for_output$start_period_assign_years[prog_rates_for_output$id_paper == "1"] <- 1985

# Subset for use in model
prog_rates_for_fitting <- prog_rates_for_output %>%
  select(id_paper,
         id_group,
         id_proc,
         pop_group_clinical,
         start_period_assign_years,
         age_assign_years,
         fu_assign_years,
         starts_with("bl_age"),
         sex,
         rate_py,
         rate_py_ci_lower,
         rate_py_ci_upper,
         py_at_risk,
         sample_size)

# Assign an outcome description which also serves as a unique ID
prog_rates_for_fitting_outcome <- c("shadow1a_eag_loss_m",
                                    "shadow1b_eag_loss_m",
                                    "shadow1a_eag_loss_f",
                                    "shadow1b_eag_loss_f",
                                    "shadow1a_hcc_incidence_m",
                                    "shadow1b_hcc_incidence_m",
                                    "shadow1a_hcc_incidence_f",
                                    "shadow1b_hcc_incidence_f",
                                    "shadow1_dcc_incidence",
                                    "shadow1_mortality_m",
                                    "shadow1_mortality_f",
                                    "shadow3_mortality",
                                    "shadow2_sag_loss",
                                    "gmb6_1_a_foi",
                                    "gmb6_1_b_foi",
                                    "gmb7_1_chronic_infection_incidence")

length(prog_rates_for_fitting_outcome) == nrow(prog_rates_for_fitting)

prog_rates_for_fitting <- cbind(outcome = prog_rates_for_fitting_outcome, prog_rates_for_fitting)

# Turn number columns into numeric format
prog_rates_for_fitting[,c(6:15,17:21)] <- apply(prog_rates_for_fitting[,c(6:15,17:21)], 2, 
                                               function(x) as.numeric(x))

# Change column names
names(prog_rates_for_fitting)[names(prog_rates_for_fitting)=="rate_py"] <- "data_value"
names(prog_rates_for_fitting)[names(prog_rates_for_fitting)=="rate_py_ci_lower"] <- "ci_lower"
names(prog_rates_for_fitting)[names(prog_rates_for_fitting)=="rate_py_ci_upper"] <- "ci_upper"

#write.csv(prog_rates_for_fitting, file = here(outpath_hbvdata, "progression_rates.csv"), row.names = FALSE)

## Mortality curves ----
subset_mortality_curves <- input_mortality_curves %>%
                            select(id_paper,
                            id_group,
                            id_proc,
                            pop_group_clinical,
                            recruitment_period,
                            dp_period,
                            country,
                            study_link,
                            study_details,
                            starts_with("bl_age"),
                            bl_hbsag_positive_prop,
                            bl_hcc_proportion,
                            bl_cirrhosis_proportion,
                            sex,
                            proportion_male,
                            dp_description,
                            model_prog_from,
                            model_prog_to,
                            time_interval_years,
                            cum_risk,
                            cum_risk_ci_lower,
                            cum_risk_ci_upper,
                            number_at_risk,
                            sample_size,
                            dp_details) 

# Assign a specific start time for the cohort (first year of recruitment)
# For A4 and A3, there is only 1 start year anyway
# For A1, assign the midpoint of the recruitment period (2012)
subset_mortality_curves$start_period_assign_years <- substr(subset_mortality_curves$recruitment_period,1,4)
subset_mortality_curves$start_period_assign_years[subset_mortality_curves$id_paper == "A1"] <- 2012

# No need to assign age because I'm assuming these are representative samples of the patient groups

# For fitting, need to filter out time intervals that aren't multiples of 0.5
# Subset for use in model
mortality_curves_for_fitting <- subset_mortality_curves %>%
  select(id_paper,
         id_group,
         id_proc,
         pop_group_clinical,
         start_period_assign_years,
         starts_with("bl_age"),
         sex,
         dp_description,
         time_interval_years,
         cum_risk,
         cum_risk_ci_lower,
         cum_risk_ci_upper,
         number_at_risk,
         sample_size) %>%
  filter(time_interval_years%%0.5 == 0)

# Assign an outcome description which also serves as a unique ID
mortality_curves_for_fitting_outcome <- c("shadow4_cum_mortality",
                                          "shadow5_cum_mortality",
                                          "shadow5_cum_mortality",
                                          "shadow5_cum_mortality",
                                          "shadow6_cum_mortality",
                                          "shadow6_cum_mortality",
                                          "shadow6_cum_hcc_incidence",
                                          "shadow6_cum_hcc_incidence")

length(mortality_curves_for_fitting_outcome) == nrow(mortality_curves_for_fitting)

mortality_curves_for_fitting <- cbind(outcome = mortality_curves_for_fitting_outcome,
                                      mortality_curves_for_fitting) %>%
  select(-dp_description, - bl_age_distribution_years)

# Turn number columns into numeric format
mortality_curves_for_fitting[,c(6:12,14:19)] <- apply(mortality_curves_for_fitting[,c(6:12,14:19)], 2, 
                                                function(x) as.numeric(x))

names(mortality_curves_for_fitting)[names(mortality_curves_for_fitting)=="cum_risk"] <- "data_value"
names(mortality_curves_for_fitting)[names(mortality_curves_for_fitting)=="cum_risk_ci_lower"] <- "ci_lower"
names(mortality_curves_for_fitting)[names(mortality_curves_for_fitting)=="cum_risk_ci_upper"] <- "ci_upper"

#write.csv(mortality_curves_for_fitting, file = here(outpath_hbvdata, "mortality_curves.csv"), row.names = FALSE)

## MTCT risk in West Africa ----
# Split up into input and output datasets
input_mtct_risk_for_input <- filter(input_mtct_risk, !grepl("output", modelling_notes))
input_mtct_risk_for_output <- filter(input_mtct_risk, grepl("output", modelling_notes))

mtct_risk_for_fitting <- select(input_mtct_risk_for_output,
                                id_paper,
                                id_group,
                                id_proc,
                                dp_period,
                                vaccinated,
                                mtct_risk_prop,
                                mtct_risk_prop_ci_lower,
                                mtct_risk_prop_ci_upper,
                                sample_size)


# Turn number columns into numeric format
mtct_risk_for_fitting[,c("dp_period", "mtct_risk_prop", "mtct_risk_prop_ci_lower",
                         "mtct_risk_prop_ci_upper", "sample_size")] <-
  apply(mtct_risk_for_fitting[,c("dp_period", "mtct_risk_prop", "mtct_risk_prop_ci_lower",
                                 "mtct_risk_prop_ci_upper", "sample_size")], 2, 
        function(x) as.numeric(x))

# Add outcome column
mtct_risk_for_fitting <- cbind(outcome = "mtct_risk",
                               mtct_risk_for_fitting) %>%
  arrange(dp_period)

# Rename columns to match model output
names(mtct_risk_for_fitting)[names(mtct_risk_for_fitting)=="dp_period"] <- "time"
names(mtct_risk_for_fitting)[names(mtct_risk_for_fitting)=="mtct_risk_prop"] <- "data_value"
names(mtct_risk_for_fitting)[names(mtct_risk_for_fitting)=="mtct_risk_prop_ci_lower"] <- "ci_lower"
names(mtct_risk_for_fitting)[names(mtct_risk_for_fitting)=="mtct_risk_prop_ci_upper"] <- "ci_upper"

#write.csv(mtct_risk_for_fitting, file = here(outpath_hbvdata, "mtct_risk.csv"), row.names = FALSE)

                              
## Risk of chronic carriage in West Africa ----
prop_chronic_for_fitting <- select(input_chronic_carriage_risk,
                                id_paper,
                                id_group,
                                id_proc,
                                age_at_infection,
                                p_chronic,
                                p_chronic_ci_lower,
                                p_chronic_ci_upper,
                                sample_size)

# Turn number columns into numeric format
prop_chronic_for_fitting[,c(4:8)] <- apply(prop_chronic_for_fitting[,c(4:8)], 2, 
                                           function(x) as.numeric(x))

# Add outcome column
prop_chronic_for_fitting <- cbind(outcome = "p_chronic",
                                  prop_chronic_for_fitting)

# Rename columns to match model output
names(prop_chronic_for_fitting)[names(prop_chronic_for_fitting)=="age_at_infection"] <- "age"
names(prop_chronic_for_fitting)[names(prop_chronic_for_fitting)=="p_chronic"] <- "data_value"
names(prop_chronic_for_fitting)[names(prop_chronic_for_fitting)=="p_chronic_ci_lower"] <- "ci_lower"
names(prop_chronic_for_fitting)[names(prop_chronic_for_fitting)=="p_chronic_ci_upper"] <- "ci_upper"

# Round ages to the nearest 0.5
prop_chronic_for_fitting$age <- floor(prop_chronic_for_fitting$age/0.5)*0.5

#write.csv(prop_chronic_for_fitting, file = here(outpath_hbvdata, "risk_of_chronic_carriage.csv"), row.names = FALSE)




## Odds ratios ----
odds_ratios <- input_odds_ratios

# Assign year as mid-point for GLCS
odds_ratios$dp_period_assign_years <- c(1999,1999,2013)

# Assign minimum and maximum years
# For GLCS this was calculated as mean +- 2*SD using controls and cases ages
odds_ratios$age_assign_min_years <- c(15,15,8)
odds_ratios$age_assign_max_years <- c(83.5,83.5,95.5)

# Save a dataset for use in model
odds_ratios_for_fitting <- select(odds_ratios,
                                id_paper,
                                id_group,
                                id_proc,
                                dp_period_assign_years,
                                age_assign_min_years,
                                age_assign_max_years,
                                proportion_male,
                                odds_ratio,
                                odds_ratio_ci_lower,
                                odds_ratio_ci_upper,
                                sample_size)

# Assign outcome description 
odds_ratios_outcome <- tolower(paste0("odds_ratio_", odds_ratios$exposure, 
                                      "_and_", odds_ratios$outcome))
odds_ratios_outcome <- gsub(" ", "_", odds_ratios_outcome)

odds_ratios_for_fitting <- cbind(outcome = odds_ratios_outcome,
                                 odds_ratios_for_fitting)

# Rename columns to match model output
names(odds_ratios_for_fitting)[names(odds_ratios_for_fitting)=="dp_period_assign_years"] <- "time"
names(odds_ratios_for_fitting)[names(odds_ratios_for_fitting)=="age_assign_min_years"] <- "age_min"
names(odds_ratios_for_fitting)[names(odds_ratios_for_fitting)=="age_assign_max_years"] <- "age_max"
names(odds_ratios_for_fitting)[names(odds_ratios_for_fitting)=="odds_ratio"] <- "data_value"
names(odds_ratios_for_fitting)[names(odds_ratios_for_fitting)=="odds_ratio_ci_lower"] <- "ci_lower"
names(odds_ratios_for_fitting)[names(odds_ratios_for_fitting)=="odds_ratio_ci_upper"] <- "ci_upper"

#write.csv(odds_ratios_for_fitting, file = here(outpath_hbvdata, "odds_ratios.csv"), row.names = FALSE)



#### INPUT DATASETS
input_prog_rates_for_input <- filter(input_progression_rates, modelling_use == "input")
input_mtct_risk_for_input
input_infant_vaccine_efficacy
input_paf_liver_disease

#write.csv(input_prog_rates_for_input, file = here(outpath_hbvdata, "input_progression_rates.csv"), row.names = FALSE)
#write.csv(input_mtct_risk_for_input, file = here(outpath_hbvdata, "input_mtct_risk.csv"), row.names = FALSE)



## Mean age and proportion male in HBV-related liver disease ----
# Assign midpoint of study period: 1999 (GLCS)
input_liver_disease_demography$dp_period_assign_years <- as.numeric(1999)

# Split proportion male and mean age into different rows
liver_disease_demography_age <- select(input_liver_disease_demography,
                                       -proportion_male)
liver_disease_demography_sex <- select(input_liver_disease_demography,
                                       -age_mean_years, 
                                       -age_mean_years_ci_lower,
                                       -age_mean_years_ci_upper)
# Add outcome for identification
liver_disease_demography_age$outcome <- c("hcc_mean_age", "cirrhosis_mean_age")
liver_disease_demography_sex$outcome <- c("hcc_prop_male", "cirrhosis_prop_male")

# Rename columns
names(liver_disease_demography_age)[
  names(liver_disease_demography_age)=="age_mean_years"] <- "data_value"
names(liver_disease_demography_age)[
  names(liver_disease_demography_age)=="age_mean_years_ci_lower"] <- "ci_lower"
names(liver_disease_demography_age)[
  names(liver_disease_demography_age)=="age_mean_years_ci_upper"] <- "ci_upper"
names(liver_disease_demography_sex)[
  names(liver_disease_demography_sex)=="proportion_male"] <- "data_value"

# Add empty confidence interval columns for proportion male and reorder into right position
liver_disease_demography_sex$ci_lower <- NA
liver_disease_demography_sex$ci_upper <- NA
liver_disease_demography_sex <- liver_disease_demography_sex[,c(1:23,29,30,24:28)]

liver_disease_demography <- rbind(liver_disease_demography_sex, liver_disease_demography_age)

# Save smaller dataset for fitting
liver_disease_demography_for_fitting <- select(liver_disease_demography,
                                               outcome,
                                               id_paper,
                                               id_group,
                                               id_proc,
                                               pop_group_clinical,
                                               dp_period_assign_years,
                                               data_value,
                                               ci_lower,
                                               ci_upper,
                                               sample_size)


# Rename columns to match model output
names(liver_disease_demography_for_fitting)[
  names(liver_disease_demography_for_fitting)=="dp_period_assign_years"] <- "time"

]write.csv(liver_disease_demography_for_fitting, file = here(outpath_hbvdata, "liver_disease_demography.csv"), row.names = FALSE)

