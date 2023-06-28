no_function()

###
setwd(masstools::get_project_wd())
rm(list = ls())
# source("code/tools.R")
library(tidyverse)
library(tidymass)

load("data_analysis/plasma_proteomics/data_preparation/object")

setwd("data_analysis/plasma_proteomics/data_preparation")

object

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
