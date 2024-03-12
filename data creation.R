# Chrystelle Kiang
# Aim 3. Hormone therapy adherence pre-/post-pandemic 
# Last updated 2024 March 6 
# this file is for dataset creation and data cleaning
library(here) 
library(tidyverse)

# Imported datasets and saved as rds for future use
# faster loading... and to prevent accidental edits to csv
# GCRcancerdata <- read.csv("/Volumes/Kiang/GCR_CancerData.csv", header = TRUE, na.strings = c("","NA"))
# saveRDS(GCRcancerdata, file = "GCRcancerdata.rds")
# GCRpharmdata <- read.csv("/Volumes/Kiang/GCR_PhamacyData.csv", header = TRUE, na.strings = c("","NA"))
# saveRDS(GCRpharmdata, file = "GCRpharmdata.rds")
# GCRstagedata <- read.csv("/Volumes/Kiang/GCR_Stage.csv", header = TRUE, na.strings = c("","NA"))
# saveRDS(GCRstagedata, file = "GCRstagedata.rds")
# GCRadditional <- read.csv("/Volumes/Kiang/GCR_AdditionalVars.csv", header = TRUE, na.strings = c("","NA"))
# saveRDS(GCRadditional, file = "./data/GCRadditional.rds")

# Reading datasets in. These are R version of csv file; don't want to edit.
GCRcancerdata <- readRDS("./data/GCRcancerdata.rds")
GCRpharmdata <- readRDS("./data/GCRpharmdata.rds")
stagedata <- readRDS("./data/GCRstagedata.rds")
GCRadditional <- readRDS("./data/GCRadditional.rds") %>%
  rename(CTCTUMOR_RECORD_NUMBER = TRN)

##### Creating analysis dataset
cancerdata <- GCRcancerdata %>%
  group_by(StudyID) %>%
  left_join(stagedata, by = c("StudyID","CTCTUMOR_RECORD_NUMBER"), relationship = "many-to-many") %>%
  left_join(GCRadditional, by = c("StudyID","CTCTUMOR_RECORD_NUMBER"), relationship = "many-to-many") %>%
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
    CTCESTROGEN_RECEPTOR_SUMMARY == 9 ~ 9), #not documented/indeterminate/unknown/etc
  # creating combined stage variable 
  stage_summary = case_when(
    # dx year 2015-2017
    CTCDATE_OF_DIAGNOSIS_YYYY %in% c(2015,2016,2017) ~ CTCDERIVED_SS2000, 
    # dx year 2018-2019
    CTCDATE_OF_DIAGNOSIS_YYYY %in% c(2018,2019) ~ CTCDERIVED_SUMMARY_STAGE_2018),
   sex = factor(PATIENTSEX, 
                        levels = c(1, 2), 
                        labels = c('Male', 'Female')),
           marital_status = factor(CTCMARITAL_STATUS_AT_DX,
                                   levels = c(1,2,3,4,5,6,9),
                                   labels = c('Single (never married)', 'Married (including common law)', 'Separated', 'Divorced', 'Widowed', 'Unmarried or Domestic Partner', 'Unknown')),
           race = factor(PATIENTRACE_1,
                         levels = c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,20,21,22,25,26,27,28,30,31,32,96,97,98,99),
                         labels = c('White','Black or African American','American Indian or Alaska Native','Chinese','Japanese','Filipino','Native Hawaiian','Korean','Vietnamese','Laotian','Hmong','Cambodian','Thai','Asian Indian or Pakistani','Asian Indian','Pakistani','Micronesian','Chamorro','Guamanian','Polynesian','Tahitian','Samoan','Tongan','Melanesian','Fiji Islander','Papua New Guinean','Other Asian','Pacific Islander','Other race','Unknown by patient')),
           hispanic = factor(PATIENTNHIA, 
                             levels = c(0,1,2,3,4,5,6,7,8),
                             labels = c('Non-Hispanic','Mexican','Puerto Rican','Cuban','South or Central American (except Brazil)','Other specified Spanish/Hispanic origin','Spanish, Hispanic, or Latino','NHIA surname match','Dominican Republic')),
           race_collapsed = case_when(
             race == "Black or African American" ~ "Black or African American",
             race == "American Indian or Alaska Native" ~ "AI, AN, or other race",
             race == "White" ~ "White",
             race == "Other race" ~ "AI, AN, or other race",
             race == "Unknown by patient" ~ "Unknown by patient",
             .default = "Asian"),
           race_ethnicity = case_when(
             (hispanic == "Non-Hispanic" & race_collapsed == "Asian") ~ "NH Asian",
             (hispanic == "Non-Hispanic" & race_collapsed == "Black or African American") ~ "NH Black or AA",
             (hispanic == "Non-Hispanic" & race_collapsed == "AI, AN, or other race") ~ "NH AI, AN, or other race",
             (hispanic == "Non-Hispanic" & race_collapsed == "White") ~ "NH White",
             (hispanic != "Non-Hispanic") ~ "Hispanic",
             TRUE ~ race_collapsed),
           laterality = factor(CTCLATERALITY,
                               levels = c(0,1,2,3,4,5,9),
                               labels = c('Not a paired site','Right: origin of primary','Left: origin of primary','Only one side involved, right or left origin unspecified','Bilateral involvement at time of diagnosis','Paired site: midline tumor','Paired site, but no information concerning laterality')),
           grade = factor(CTCGRADE,
                          levels = c(1,2,3,4,5,6,7,8,9), 
                          labels = c('Grade I','Grade II','Grade III','Grade IV','T-cell','B-cell','Null cell','Natural killer (NK) cetabll', 'Unknown, not stated, or NA')),
           reporting_source = factor(CTCTYPE_OF_REPORTING_SOURCE,
                                     levels = c(1,2,3,4,5,6,7,8),
                                     labels = c("Hospital inpatient","Radiation Treatment Centers or Medical Oncology Centers","Laboratory","Physician's office/private medical practitioner","Nursing/convalescent home/hospice","Autopsy","Death certificate","Other hospital outpatient units/surgery centers")),
           estrogen_receptor = factor(ERstatus,
                                      levels = c(0,1,7,9),
                                      labels = c('ER negative','ER positive','Results not in chart','Not documented'))
      )

# subset to population of interest: stages I - III, ER positive
cancer_eligible <- cancerdata %>% 
  group_by(StudyID) %>%
  filter(ERstatus == 1, 
         stage_summary %in% c(1,2,3))
saveRDS(cancer_eligible, file = "./data/cancerstudypop.rds")
# inefficient but doing this so I can refer back to above more easily
# in case I need to look more into sequence vs. tumor number filtering
cancer_eligible <- cancer_eligible %>%
  filter(SEQNo == min(SEQNo)) 

GCRpharmdata <- GCRpharmdata %>%
  mutate(dispense_dt = ymd(dispense_dt),
         rx_written_dt = ymd(rx_written_dt),
    drug_class_collapsed = case_when(
    # not using drug class, but wrote this and don't have reason to delete 
    startsWith(therapeutic_drug_class, "5-HT3") ~ "5-HT3 RECEPTOR ANTAGONISTS",
    startsWith(therapeutic_drug_class, "ANTINEOPLASTIC") ~ "ANTINEOPLASTIC AGENTS",
    startsWith(therapeutic_drug_class, "BONE") ~ "BONE RESORPTION INHIBITORS",
    therapeutic_drug_class %in% c("ESTROGEN","ESTROGENS") ~ "ESTROGEN",
    startsWith(therapeutic_drug_class, "GI") ~ "GI DRUGS, MISC",
    startsWith(therapeutic_drug_class, "HEMA") ~ "HEMATOPOIETIC AGENTS",
    startsWith(therapeutic_drug_class, "PROKINETIC") ~ "PROKINETIC AGENTS",
    startsWith(therapeutic_drug_class, "SKIN") ~ "SKIN AND MUCOUS MEMBRANE AGENTS, MISC",
    TRUE ~ therapeutic_drug_class),
  drug_name = case_when(
    startsWith(generic_name, "ANASTROZOLE") ~ "ANASTROZOLE",
    startsWith(generic_name, "EXEMESTANE") ~ "EXEMESTANE",
    startsWith(generic_name, "LETROZOLE") ~ "LETROZOLE",
    startsWith(generic_name, "TAMOXIFEN") ~ "TAMOXIFEN",
    TRUE ~ "OTHER"),
  drug_group = case_when(
    startsWith(generic_name, "ANASTROZOLE") ~ "AI",
    startsWith(generic_name, "EXEMESTANE") ~ "AI",
    startsWith(generic_name, "LETROZOLE") ~ "AI",
    startsWith(generic_name, "TAMOXIFEN") ~ "TAMOXIFEN",
    TRUE ~ "OTHER"))

# subsetting GCR pharm data to relevant meds before join
start_date <- as.Date("2015-01-01", format = "%Y-%m-%d")
GCRpharmdata <- GCRpharmdata %>%
  mutate(elig_rx = case_when(
    as.Date(dispense_dt) >= start_date ~ 1,
    TRUE ~ 0))

AETpharmdata <- GCRpharmdata %>%
   filter(drug_group %in% c("AI","TAMOXIFEN")) 
# not efficient but I was having trouble with this 

# left join to keep rows in AETpharm that have StudyIDs that ARE in cancer_eligible
# this keeps those eligible who do not have pharm data so we can report the total vs. who have data on
cancer_pharm_data <- cancer_eligible %>%
  left_join(AETpharmdata, by = "StudyID", relationship = "many-to-many") 

# updating labels based on data dictionary
cancer_pharm_data <- cancer_pharm_data %>%
  mutate(AETdrug = case_when(drug_group %in% c("AI","TAMOXIFEN") ~ "ET"),
         dx_month = case_when(CTCDATE_OF_DIAGNOSIS_MM == 99 ~ 1,
                              TRUE ~ CTCDATE_OF_DIAGNOSIS_MM),
         dx_day = case_when(CTCDATE_OF_DIAGNOSIS_DD == 99 ~ 1,
                            TRUE ~ CTCDATE_OF_DIAGNOSIS_DD),
         dt_diagnosis = as.Date(
           paste(CTCDATE_OF_DIAGNOSIS_YYYY,
                 sprintf("%02d", as.numeric(dx_month)),
                 sprintf("%02d", as.numeric(dx_day)),
                 sep = "-"),
           format = "%Y-%m-%d")
  ) %>%
    select(-dx_month, -dx_day)
# note that I used excel and saved these as csv then copy/pasted from txt file 

# Save final analytic dataset that will be used
saveRDS(cancer_pharm_data, file = "studypopdata.rds")
