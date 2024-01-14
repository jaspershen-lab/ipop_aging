no_source()

rm(list = ls())

setwd(masstools::get_project_wd())
library(tidyverse)

load("3-data_analysis/plasma_proteomics/data_preparation/object_cross_section_loess")

variable_info <-
  object_cross_section_loess@variable_info

grep("MMP12", variable_info$variable_id)


object_cross_section_loess %>% 
  activa
