no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
masstools::setwd_project()

source("code/tools.R")

load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("data_analysis/combined_omics/DE_SWAN")
setwd("data_analysis/combined_omics/DE_SWAN")

library("DEswan")

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id))

expression_data <-
  object_cross_section_loess@expression_data

sample_info <-
  object_cross_section_loess@sample_info

variable_info <-
  object_cross_section_loess@variable_info

load("temp_data")

temp_data <-
  temp_data %>%
  dplyr::left_join(variable_info,
                   by = c("variable_id"))

###Pathway enrichment
temp_data_clinical_test <-
  temp_data %>%
  dplyr::filter(class == "clinical_test" & p_value_adjust < 0.05)

####clinical_test
dir.create("clinical_test_pathway_crest1")

####crest 1
temp_data_clinical_test_crest1 <-
  temp_data_clinical_test %>%
  dplyr::filter(center == 51)

dim(temp_data_clinical_test_crest1)

length(unique(temp_data_clinical_test$variable_id))


library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
addWorksheet(wb, sheetName = "clinical_tests in crest1", gridLines = TRUE)
freezePane(wb,
           sheet = 1,
           firstRow = TRUE,
           firstCol = TRUE)
writeDataTable(
  wb,
  sheet = 1,
  x = temp_data_clinical_test_crest1 %>% dplyr::select(-Class),
  colNames = TRUE,
  rowNames = FALSE
)
saveWorkbook(
  wb,
  "lipidomics_pathway_crest1/lipid_minion_result/lipid_enrichment.xlsx",
  overwrite = TRUE
)
