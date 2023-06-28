# no_source()
# rm(list = ls())
# library(tidyverse)
# library(tidymass)
# setwd(masstools::get_project_wd())
# 
# source("code/tools.R")
# 
# load(
#   "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
# )
# 
# dir.create("data_analysis/combined_omics/DE_SWAN")
# setwd("data_analysis/combined_omics/DE_SWAN")
# 
# library("DEswan")
# 
# object_cross_section_loess <-
#   object_cross_section_loess %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::mutate(age = as.numeric(sample_id))
# 
# expression_data <-
#   object_cross_section_loess@expression_data
# 
# sample_info <-
#   object_cross_section_loess@sample_info
# 
# variable_info <-
#   object_cross_section_loess@variable_info
# 
# 
# ###crest 1
# object_cross_section_loess@sample_info$age
# sum(
#   object_cross_section_loess@sample_info$age >= 40 &
#     object_cross_section_loess@sample_info$age <= 47
# )
# 
# sum(
#   object_cross_section_loess@sample_info$age >= 57 &
#     object_cross_section_loess@sample_info$age <= 65
# )
# 
# hist(object_cross_section_loess@sample_info$age)
# 
# object <-
#   object_cross_section_loess %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(age >= 40 & age <= 47)
# 
# cor_data <-
#   massstat::cor_mass_dataset(
#     x = object,
#     margin = "variable",
#     method = "spearman",
#     data_type = "longer",
#     p_adjust_method = "BH"
#   )
# 
# 
# cor_data %>% 
#   dplyr::filter(p_adjust < 0.001) %>% 
#   head()
# 
# plot(unlist(object["transcriptome_A1CF",,drop = TRUE]),
#      unlist(object["transcriptome_A1BG",,drop = TRUE]))
