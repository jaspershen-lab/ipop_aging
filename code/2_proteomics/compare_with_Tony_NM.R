no_source()

rm(list = ls())

masstools::setwd_project()
library(tidyverse)

load("data_analysis/plasma_proteomics/data_preparation/object_cross_section_loess")

variable_info <-
  object_cross_section_loess@variable_info

grep("MMP12", variable_info$variable_id)


object_cross_section_loess %>% 
  activa
