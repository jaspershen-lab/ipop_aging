no_function()

###
setwd(r4projects::get_project_wd())
rm(list = ls())

# source("1-code/100-tools.R")
library(tidyverse)
library(tidymass)

setwd("3-data_analysis/oral_microbiome/data_preparation")

load("object")

dim(object)

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(adjusted_age))

####only remain the health vist
object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(CL4 == "Healthy")

###for each participant, calculate the median value
object_cross_section <-
  massdataset::summarise_samples(object = object,
                                 what = "mean_intensity",
                                 group_by = "subject_id")

dim(object_cross_section)

save(object_cross_section, file = "object_cross_section")
