###########################################################
### Imperial HBV model:                                  ###
### Clean input natural history data for The Gambia     ###
### Source: mapping rev                                 ###
###########################################################
# Load packages and set directories
require(tidyr)  # for data processing
require(dplyr)  # for data processing
require(here)  # for setting working directory
require(ggplot2)
inpath_hbvdata <- "data-raw"

## HBsAg prevalence in The Gambia dataset
# Age- and sex-specific datasets
input_hbsag_prev <- read.csv(here(inpath_hbvdata,
                                  "hbsag_prevalence.csv"),
                                  header = TRUE, check.names = FALSE,
                                  stringsAsFactors = FALSE)


subset_hbsag_prev <- select(input_hbsag_prev,
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
                            sample_size) %>%
  filter(hbsag_positive_prop != "DUPLICATE" & hbsag_positive_prop != "NR")

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

# 2) Assign to the pre- or post-vaccination period (1991)
# Split datapoint collection period column into minimum and maximum year
subset_hbsag_prev <- subset_hbsag_prev %>%
  separate(col = dp_period, into = c("dp_period_min", "dp_period_max"),
           sep = "-", remove = FALSE) %>%
  mutate(dp_period_max = coalesce(dp_period_max, dp_period_min))
# Assign post-vaccination status to data points collected in 1991 or after
# and pre-vaccination status to data points collected before 1991
subset_hbsag_prev$dp_period_vacc <- "post-vacc"
subset_hbsag_prev$dp_period_vacc[subset_hbsag_prev$dp_period_max < 1991] <- "pre-vacc"

# 3) Assign data points to a series based on the study link (if available) or the IDs
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
  geom_text(aes(label = dp_period), vjust=1.5) +
  theme_bw() + ylim(0, 0.4) + xlim(0,80) + ggtitle("Post-1991 data points") +
  ylab("HBsAg prevalence (proportion)") + xlab("Age (years)")  +
  scale_color_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA",
                              "#3F4921", "#D1A33D", "#8569D5", "#5F7FC7",
                              "#673770", "#38333E"))

## Natural history progression rates in West Africa
input_progression_rates <- read.csv(here(inpath_hbvdata,
                                         "natural_history_progression_rates.csv"),
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)

subset_progression_rates <- input_progression_rates %>%
  select(                   id_paper,
                            id_group,
                            id_proc,
                            pop_group_clinical,
                            pop_group_demographic,
                            recruitment_period,
                            dp_period,
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
                            modelling_use,
                            modelling_notes) %>%
  mutate(rate_py = rate_100py/100,  # Convert rate per 100 person-years to per person-year
         rate_py_ci_lower = as.numeric(rate_100py_ci_lower)/100,
         rate_py_ci_upper = as.numeric(rate_100py_ci_upper)/100)

# Split into 2 datasets based on use as input or output within model
prog_rates_for_input <- filter(subset_progression_rates, modelling_use == "input")
prog_rates_for_output <- filter(subset_progression_rates, modelling_use == "output")

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

# Subset for use in model
prog_rates_for_fitting <- prog_rates_for_output %>%
  select(id_paper,
         id_group,
         id_proc,
         start_period_assign_years,
         age_assign_years,
         fu_assign_years,
         starts_with("bl_age"),
         sex,
         pop_group_clinical,
         model_prog_from,
         model_prog_to,
         rate_py,
         rate_py_ci_lower,
         rate_py_ci_upper,
         modelling_notes)

prog_rates_for_fitting$numerator <- c("cum. incident transitions to IC and ENCHB",
                                            "cum. incident HCC cases",
                                            "cum. incident HCC cases",
                                            "cum. incident HCC cases",
                                            "cum. incident DCC cases - cum. transitions from DCC to HCC",
                                            "cum. incident deaths from CC, DCC, HCC and background",
                                            "cum. incident deaths from CC, DCC, HCC and background",
                                            "cum. incident deaths from CC, DCC, HCC and background",
                                            "cum. incident transitions from IC to R",
                                            "cum. incident transitions from S to IT and S to R",
                                            "cum. incident transitions from S to IT and S to R",
                                            "cum. incident transitions from S to IT")
prog_rates_for_fitting$denominator <- c("personyears in IT and IR",
                                             "personyears in chronic compartments except HCC",
                                             "personyears in chronic compartments except HCC",
                                             "personyears in chronic compartments except HCC",
                                             "personyears in chronic compartments except DCC and HCC",
                                             "personyears in chronic compartments",
                                             "personyears in chronic compartments",
                                             "personyears in CC, DCC and HCC",
                                             "personyears in IC",
                                             "personyears in S",
                                             "personyears in S",
                                             "personyears in S")

