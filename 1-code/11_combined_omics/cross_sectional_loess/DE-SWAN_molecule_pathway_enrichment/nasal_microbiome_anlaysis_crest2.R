no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
setwd(masstools::get_project_wd())

source("1-code/100-tools.R")

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("3-data_analysis/combined_omics/DE_SWAN")
setwd("3-data_analysis/combined_omics/DE_SWAN")

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
temp_data_nasal_microbiome <-
  temp_data %>%
  dplyr::filter(class == "nasal_microbiome" & p_value_adjust < 0.05)

####nasal_microbiome
dir.create("nasal_microbiome_pathway_crest2")

####crest 2
temp_data_nasal_microbiome_crest2 <-
  temp_data_nasal_microbiome %>%
  dplyr::filter(center == 58)

dim(temp_data_nasal_microbiome_crest2)

length(unique(temp_data_nasal_microbiome$variable_id))

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
addWorksheet(wb, sheetName = "nasal_microbiomes in crest2", gridLines = TRUE)
freezePane(wb,
           sheet = 1,
           firstRow = TRUE,
           firstCol = TRUE)
writeDataTable(
  wb,
  sheet = 1,
  x = temp_data_nasal_microbiome_crest2 %>% dplyr::select(-Class),
  colNames = TRUE,
  rowNames = FALSE
)
saveWorkbook(
  wb,
  "nasal_microbiome_pathway_crest2/nasal_microbiome_makrers.xlsx",
  overwrite = TRUE
)
