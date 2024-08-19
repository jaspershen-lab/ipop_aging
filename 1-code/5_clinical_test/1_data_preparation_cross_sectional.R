no_function()

###
setwd(r4projects::get_project_wd())
rm(list = ls())

# source("1-code/100-tools.R")
library(tidyverse)
library(tidymass)

setwd("3-data_analysis/clinical_test/data_preparation")

load("object")


massdataset::export_mass_dataset(object = object, 
                                 file_type = "xlsx")

sum(object < 0, na.rm = TRUE)

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

# massdataset::export_mass_dataset(object = object, file_type = "xlsx")





###for each participant, calculate the median value
object_cross_section <-
  massdataset::summarise_samples(object = object,
                                 what = "mean_intensity",
                                 group_by = "subject_id")

dim(object_cross_section)



object_cross_section

dim(object_cross_section)

sum(is.na(object_cross_section))

show_variable_missing_values(object_cross_section, percentage = TRUE)

####remove variables has lots of NA
object_cross_section <-
  object_cross_section %>%
  mutate_variable_na_freq()

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(na_freq < 0.5)

show_sample_missing_values(object_cross_section, percentage = TRUE)

####remove samples has lots of NA
object_cross_section <-
  object_cross_section %>%
  mutate_sample_na_freq()

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(na_freq < 0.5)

object_cross_section <-
  object_cross_section %>%
  impute_mv(method = "knn")

sum(is.na(object_cross_section))

save(object_cross_section, file = "object_cross_section")
