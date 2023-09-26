# Chrystelle Kiang
# Aim 3. Hormone therapy adherence pre-/post-pandemic 
# September 2023
# this file is for dataset creation and prelim data exploration/cleaning
library(here) 
library(tidyverse)
library(ggplot2)
library(arsenal)
# after installing new package, update lock file
# renv::snapshot()

# Imported datasets and saved as rds for future use
# faster loading... and in case something happens to csv 
# GCRcancerdata <- read.csv("/Volumes/Kiang/GCR_CancerData.csv", header = TRUE, na.strings = c("","NA"))
# saveRDS(GCRcancerdata, file = "GCRcancerdata.rds")
# GCRpharmdata <- read.csv("/Volumes/Kiang/GCR_PhamacyData.csv", header = TRUE, na.strings = c("","NA"))
# saveRDS(GCRpharmdata, file = "GCRpharmdata.rds")
# GCRstagedata <- read.csv("/Volumes/Kiang/GCR_Stage.csv", header = TRUE, na.strings = c("","NA"))
# saveRDS(GCRstagedata, file = "GCRstagedata.rds")


# Reading datasets in. These are R version of csv file; don't want to edit.
GCRcancerdata <- readRDS("GCRcancerdata.rds")
GCRpharmdata <- readRDS("GCRpharmdata.rds")
stagedata <- readRDS("GCRstagedata.rds")

table(stagedata$CTCDERIVED_AJCC_6_STAGE_GRP, useNA = "always")
table(stagedata$CTCDERIVED_SUMMARY_STAGE_2018, useNA = "always")

table(cancerdata$CTCDERIVED_AJCC_6_STAGE_GRP, cancerdata$CTCDATE_OF_DIAGNOSIS_YYYY, useNA = "always")
table(cancerdata$CTCDERIVED_AJCC_7_STAGE_GRP, cancerdata$CTCDATE_OF_DIAGNOSIS_YYYY, useNA = "always")

table(cancerdata$CTCDERIVED_AJCC_7_STAGE_GRP, cancerdata$CTCDERIVED_AJCC_6_STAGE_GRP, useNA = "always")

##### Creating analysis dataset: 
cancerdata <- GCRcancerdata %>%
  left_join(stagedata, by = "StudyID", relationship = "many-to-many")
# looks like there are some repeats?? 

# Creating a new ER variable based on year of diagnosis 
# following CTCESTROGEN_RECEPTOR_SUMMARY convention
cancerdata <- cancerdata %>%
  mutate(ERstatus = case_when(
    # dx year 2015-2017
    CTCCS_SITE_SPECIFIC_FACTOR1 == 10 ~ 1, #positive 
    CTCCS_SITE_SPECIFIC_FACTOR1 == 20 ~ 0, #negative
    CTCCS_SITE_SPECIFIC_FACTOR1 == 30 ~ 9, #borderline aka indeterminate
    CTCCS_SITE_SPECIFIC_FACTOR1 == 997 ~ 7, #test ordered, not in chart 
    CTCCS_SITE_SPECIFIC_FACTOR1 %in% c(996, 998, 999) ~ 9, #not documented/indeterminate/unknown/etc
    # dx year 2018-2019
    CTCESTROGEN_RECEPTOR_SUMMARY == 0 ~ 0, #negative
    CTCESTROGEN_RECEPTOR_SUMMARY == 1 ~ 1, #positive
    CTCESTROGEN_RECEPTOR_SUMMARY == 7 ~ 7, #test ordered, not in chart
    CTCESTROGEN_RECEPTOR_SUMMARY == 9 ~ 9 #not documented/indeterminate/unknown/etc
  ))

# subset to population of interest:stages I - III, ER positive
# TODO update once dataset is updated

# inner join only retains observations that exist in both datasets 
# will need to take note of number included vs. not
cancer_pharm_data <- GCRpharmdata %>%
  inner_join(cancer_eligible, by = "StudyID", relationship = "many-to-many") %>%
  mutate(dispense_dt = ymd(dispense_dt),
         rx_written_dt = ymd(rx_written_dt))
# note that I am first joining then filtering to ensure consistency in the dataset
# (vs. filtering pharm and then doing different types of joins.)

# setting a start date, this is arbitrary
start_date <- ymd("2018-01-01")

# subset and cleaning
# subsetting to only people whose dispense date starts after 2018 Jan 1
cancer_pharm_data <- cancer_pharm_data %>%
  filter(dispense_dt > start_date) %>%
  mutate(drug_class_collapsed = case_when(
    startsWith(therapeutic_drug_class, "5-HT3") ~ "5-HT3 RECEPTOR ANTAGONISTS",
    startsWith(therapeutic_drug_class, "ANTINEOPLASTIC") ~ "ANTINEOPLASTIC AGENTS",
    startsWith(therapeutic_drug_class, "BONE") ~ "BONE RESORPTION INHIBITORS",
    therapeutic_drug_class %in% c("ESTROGEN","ESTROGENS") ~ "ESTROGEN",
    startsWith(therapeutic_drug_class, "GI") ~ "GI DRUGS, MISC",
    startsWith(therapeutic_drug_class, "HEMA") ~ "HEMATOPOIETIC AGENTS",
    startsWith(therapeutic_drug_class, "PROKINETIC") ~ "PROKINETIC AGENTS",
    startsWith(therapeutic_drug_class, "SKIN") ~ "SKIN AND MUCOUS MEMBRANE AGENTS, MISC",
    TRUE ~ therapeutic_drug_class)
    )

cancer_eligible %>%
  distinct(StudyID) %>%
  n_distinct()

dupes <- cancer_pharm_data %>%
  

# Save final analytic dataset that will actually be used.
saveRDS(cancer_pharm_data, file = "studypopdata.rds")
############# end of dataset creation code

# Below: data cleaning code

# TODO filter by drug category or by drug class?
# filter after saving dataset
cancer_pharm_data <- cancer_pharm_data %>%
  filter(drug_category == "HORMONE")

# TODO: how to code the actual adherence?
# will be a separate R code file
# MPR: days supply/days interval 

############# Repeat StudyIDs
# TODO investigate these people
cancer_eligible <- cancerdata %>% 
  filter(ERstatus == 1, CTCDERIVED_SUMMARY_STAGE_2018 %in% c(1,2,3))

repeat_cancer_IDs <- cancer_eligible %>% 
  group_by(StudyID) %>%
  filter(n() > 1) %>%
  summarise(n = n())

repeatedIDs <- repeat_cancer_IDs %>%
  select(StudyID)

repeated_cancer <- cancer_eligible %>%
  semi_join(repeatedIDs, by = "StudyID")

repeated_cancer %>%
  distinct(StudyID) %>%
  n_distinct()
# looks like they might just differ in tumor #?

################ archive code
# table of brand_name and hormone
cancer_pharm_data %>% 
  distinct(therapeutic_drug_class)

table(cancer_pharm_data$therapeutic_drug_class, useNA = "always")
# figuring out what to use for missing drug class
missing_drug_class <- cancer_pharm_data %>%
  filter(is.na(therapeutic_drug_class))
missing_drug_class %>%
  distinct(StudyID) %>%
  n_distinct()
table(missing_drug_class$brand_name, useNA = "always")
table(missing_drug_class$generic_name, useNA = "always")
# have generic for most 

# creating indicator for number of drug classes
# this should go in the results eventually
temp <- cancer_pharm_data %>%
  group_by(StudyID) %>%
  summarise(n_drug_class = n_distinct(therapeutic_drug_class),
            n_drug_category = n_distinct(drug_category))

table(temp$n_drug_class)
table(temp$n_drug_category)

# adding n_drug_class back to dataset then subset to unique combinations
drug_data <- cancer_pharm_data %>%
  inner_join(temp, by = "StudyID") %>%
  group_by(StudyID, drug_category, therapeutic_drug_class) %>%
  slice(1) %>%
  ungroup()

# combining repeated drug classes
# not sure this is the best way, but setting them equal then overwriting repeats
drug_data <- drug_data %>%
  mutate(drug_class_collapsed = case_when(
    startsWith(therapeutic_drug_class, "5-HT3") ~ "5-HT3 RECEPTOR ANTAGONISTS",
    startsWith(therapeutic_drug_class, "ANTINEOPLASTIC") ~ "ANTINEOPLASTIC AGENTS",
    startsWith(therapeutic_drug_class, "BONE") ~ "BONE RESORPTION INHIBITORS",
    therapeutic_drug_class %in% c("ESTROGEN","ESTROGENS") ~ "ESTROGEN",
    startsWith(therapeutic_drug_class, "GI") ~ "GI DRUGS, MISC",
    startsWith(therapeutic_drug_class, "HEMA") ~ "HEMATOPOIETIC AGENTS",
    startsWith(therapeutic_drug_class, "PROKINETIC") ~ "PROKINETIC AGENTS",
    startsWith(therapeutic_drug_class, "SKIN") ~ "SKIN AND MUCOUS MEMBRANE AGENTS, MISC",
    TRUE ~ therapeutic_drug_class
  )
  )

# tables of drug categories and etc.
drug_table1 <- tableby(drug_category ~ drug_class_collapsed, data = drug_data)
cat_by_class <- as.data.frame(summary(drug_table1, text = TRUE))
write.csv(cat_by_class, "category by class.csv")

temp_drug <- filter(drug_data, is.na(therapeutic_drug_class))
drug_table2 <- tableby(drug_category ~ generic_name, data = temp_drug)
cat_by_generic <- as.data.frame(summary(drug_table2, text = TRUE))
write.csv(cat_by_generic, "category by generics for NA class.csv")

drug_class_name <- cancer_pharm_data %>%
  group_by(drug_category, drug_class_collapsed) %>%
  distinct(generic_name)
write.csv(drug_class_name, file = "category class and name.csv")


# code used to get info, not create datasets
combined_data %>%
  group_by(drug_category) %>%
  summarise(n_distinct(StudyID))
# lol about 7,000 have 'hormone' 

table(combined_data$drug_category, combined_data$CTCRX_SUMM_HORMONE)
# not all who are prescribed have SUMM_HORMONE = 1
cancer_subset %>% 
  distinct(StudyID) %>%
  n_distinct()
