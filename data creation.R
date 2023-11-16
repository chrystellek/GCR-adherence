# Chrystelle Kiang
# Aim 3. Hormone therapy adherence pre-/post-pandemic 
# September 2023
# this file is for dataset creation and prelim data exploration/cleaning
library(here) 
library(tidyverse)
library(arsenal) # used for frequency tables 
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

##### Creating analysis dataset: 
cancerdata <- GCRcancerdata %>%
  filter(CTCTUMOR_RECORD_NUMBER == 1) %>%
  select(-CTCTUMOR_RECORD_NUMBER) %>%
  # need to verify
  left_join(stagedata, by = "StudyID", relationship = "many-to-many") %>%
  # Creating a new ER variable based on year of diagnosis 
  # following CTCESTROGEN_RECEPTOR_SUMMARY convention
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
    CTCESTROGEN_RECEPTOR_SUMMARY == 9 ~ 9) #not documented/indeterminate/unknown/etc
  ) %>%
    # filtering to only one line per patient
  group_by(StudyID) %>%
  top_n(1, CTCTUMOR_RECORD_NUMBER)
#    filter(CTCTUMOR_RECORD_NUMBER == max(CTCTUMOR_RECORD_NUMBER)) 

# subset to population of interest:stages I - III, ER positive
cancer_eligible <- cancerdata %>% 
  filter(ERstatus == 1, CTCDERIVED_SUMMARY_STAGE_2018 %in% c(1,2,3))

# inner join only retains observations that exist in both datasets 
# will need to take note of number included vs. not
cancer_pharm_data <- GCRpharmdata %>%
  inner_join(cancer_eligible, by = "StudyID", relationship = "many-to-many") %>%
  mutate(dispense_dt = ymd(dispense_dt),
         rx_written_dt = ymd(rx_written_dt))
# note that I am first joining then filtering to ensure consistency in the dataset
# (vs. filtering pharm and then doing different types of joins.)

# setting a start date, modifiable
start_date <- ymd("2018-01-01")

# subset and cleaning
# subsetting to only people whose dispense date starts after 2018 Jan 1
cancer_pharm_data <- cancer_pharm_data %>%
  filter(dispense_dt > start_date) %>%
  # not using drug class, but I wrote this so why not keep i guess.
  mutate(drug_class_collapsed = case_when(
    startsWith(therapeutic_drug_class, "5-HT3") ~ "5-HT3 RECEPTOR ANTAGONISTS",
    startsWith(therapeutic_drug_class, "ANTINEOPLASTIC") ~ "ANTINEOPLASTIC AGENTS",
    startsWith(therapeutic_drug_class, "BONE") ~ "BONE RESORPTION INHIBITORS",
    therapeutic_drug_class %in% c("ESTROGEN","ESTROGENS") ~ "ESTROGEN",
    startsWith(therapeutic_drug_class, "GI") ~ "GI DRUGS, MISC",
    startsWith(therapeutic_drug_class, "HEMA") ~ "HEMATOPOIETIC AGENTS",
    startsWith(therapeutic_drug_class, "PROKINETIC") ~ "PROKINETIC AGENTS",
    startsWith(therapeutic_drug_class, "SKIN") ~ "SKIN AND MUCOUS MEMBRANE AGENTS, MISC",
    TRUE ~ therapeutic_drug_class)) %>%
  mutate(drug_name = case_when(
    startsWith(generic_name, "ANASTROZOLE") ~ "ANASTROZOLE",
    startsWith(generic_name, "EXEMESTANE") ~ "EXEMESTANE",
    startsWith(generic_name, "LETROZOLE") ~ "LETROZOLE",
    startsWith(generic_name, "TAMOXIFEN") ~ "TAMOXIFEN",
    TRUE ~ "OTHER")) %>%
  mutate(drug_group = case_when(
    startsWith(generic_name, "ANASTROZOLE") ~ "AI",
    startsWith(generic_name, "EXEMESTANE") ~ "AI",
    startsWith(generic_name, "LETROZOLE") ~ "AI",
    startsWith(generic_name, "TAMOXIFEN") ~ "TAMOXIFEN",
    TRUE ~ "OTHER")) %>%
  # trying to create a sequential ID that isn't linked to StudyID
  group_by(StudyID) %>% 
  mutate(PersonID = row_number())

# TODO update variables to meaningful values
# studypopdata <- readRDS("studypopdata.rds")


# Save final analytic dataset that will actually be used.
saveRDS(cancer_pharm_data, file = "studypopdata.rds")
############# end of dataset creation code

# Some data exploration
cancer_eligible %>%
  distinct(StudyID) %>%
  n_distinct()

# TODO how to subset further? 
# Should perhaps do that in quarto doc so we can export tables and such 

cancer_pharm_data %>%
  distinct(StudyID) %>%
  n_distinct()

# TODO: how to code the actual adherence?
# will be a separate R code file
# MPR: days supply/days interval 

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


# tables of drug categories and etc.
# plain table
category_class <- as.data.frame(table(drug_data$drug_category, drug_data$drug_class_collapsed))
write.csv(category_class, "category and class.csv")

temp_drug <- filter(drug_data, is.na(therapeutic_drug_class))
names_no_class <- table(temp_drug$drug_category, temp_drug$generic_name)
write.csv(names_no_class, "category and name for no class.csv")


# using arsenal
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