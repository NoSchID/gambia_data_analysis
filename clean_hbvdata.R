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

## HBsAg and anti-HBc prevalence in The Gambia dataset ----
# Age- and sex-specific datasets
input_hbsag_antihbc_prev <- read.csv(here(inpath_hbvdata,
                                  "hbsag_prevalence.csv"),
                                  header = TRUE, check.names = FALSE,
                                  stringsAsFactors = FALSE)


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

hbsag_dataset_for_fitting <- cbind(prevalence_outcome = "HBsAg_prevalence",
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
names(hbsag_dataset_for_fitting)[names(hbsag_dataset_for_fitting)=="hbsag_positive_prop"] <- "value"
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

antihbc_dataset_for_fitting <- cbind(prevalence_outcome = "Anti_HBc_prevalence",
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
names(antihbc_dataset_for_fitting)[names(antihbc_dataset_for_fitting)=="antihbc_positive_prop"] <- "value"
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
# Age- and sex-specific dataset
# This is only for HBeAg prevalence in carriers overall
# (there is a separate dataset for liver disease patients)
input_hbeag_prev <- read.csv(here(inpath_hbvdata,
                                  "hbeag_prevalence.csv"),
                             header = TRUE, check.names = FALSE,
                             stringsAsFactors = FALSE)


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

hbeag_dataset_for_fitting <- cbind(prevalence_outcome = "HBeAg_prevalence",
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
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="hbeag_positive_prop"] <- "value"
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="hbeag_positive_prop_ci_lower"] <- "ci_lower"
names(hbeag_dataset_for_fitting)[names(hbeag_dataset_for_fitting)=="hbeag_positive_prop_ci_upper"] <- "ci_upper"

#write.csv(hbeag_dataset_for_fitting, file = here(outpath_hbvdata, "hbeag_prevalence.csv"), row.names = FALSE)



## Natural history prevalence dataset ----
# Here the age and timepoint assignment is different from sAg/anti-HBc/eAg datasets!
input_natural_history_prev <- read.csv(here(inpath_hbvdata,
                                          "natural_history_prevalence.csv"),
                                     header = TRUE, check.names = FALSE,
                                     stringsAsFactors = FALSE)

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
                                    sample_size) %>%
  filter(is.na(model_denominator) == FALSE)

# Turn number columns into numeric format
natural_history_prev_for_fitting[,8:13] <- apply(natural_history_prev_for_fitting[,8:13], 2, 
                                                 function(x) as.numeric(x))

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

# Replace spaces with underscores
natural_history_prev_for_fitting$prevalence_outcome <- gsub("[[:space:]]", "_", 
                                                         natural_history_prev_for_fitting$prevalence_outcome)

# Add output IDs based on numerator to ensure correct mapping:
# Only separate compartment names by underscores and convert to lowercase
natural_history_prev_for_fitting$id_output <- gsub("_prevalence.*", "", natural_history_prev_for_fitting$prevalence_outcome)
natural_history_prev_for_fitting$id_output <- gsub("and_", "", natural_history_prev_for_fitting$id_output )
natural_history_prev_for_fitting$id_output <- gsub(",", "", natural_history_prev_for_fitting$id_output)
natural_history_prev_for_fitting$id_output <- tolower(natural_history_prev_for_fitting$id_output)

# Create a unique identifier for each row from paper ID, group ID, timepoint and output ID
natural_history_prev_for_fitting$id_unique <- paste0(natural_history_prev_for_fitting$id_paper, "_",
                                                     natural_history_prev_for_fitting$id_group, "_",
                                                     natural_history_prev_for_fitting$dp_period_assign_years, "_",
                                                     natural_history_prev_for_fitting$id_output)
natural_history_prev_for_fitting$id_output <- NULL

# Round ages down to the nearest 0.5 to match model output
natural_history_prev_for_fitting$age_assign_min_years <- floor(natural_history_prev_for_fitting$age_assign_min_years/0.5)*0.5
natural_history_prev_for_fitting$age_assign_max_years <- floor(natural_history_prev_for_fitting$age_assign_max_years/0.5)*0.5

# Rename columns to match model output
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="dp_period_assign_years"] <- "time"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="age_assign_min_years"] <- "age_min"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="age_assign_max_years"] <- "age_max"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="prevalence"] <- "value"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="prevalence_ci_lower"] <- "ci_lower"
names(natural_history_prev_for_fitting)[names(natural_history_prev_for_fitting)=="prevalence_ci_upper"] <- "ci_upper"

#write.csv(natural_history_prev_for_fitting, file = here(outpath_hbvdata, "natural_history_prevalence.csv"), row.names = FALSE)


## Natural history progression rates in West Africa ----
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

