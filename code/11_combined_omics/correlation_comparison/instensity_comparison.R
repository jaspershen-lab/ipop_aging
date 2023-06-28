no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("code/tools.R")

library(tidyverse)
library(tidymass)

###load("data)
###transcriptome
load(
  "data_analysis/plasma_transcriptome/data_preparation/same_samples/object_cross_section_loess"
)

load("data_analysis/plasma_transcriptome/data_preparation/object_cross_section")

dim(object_cross_section_loess)
dim(object_cross_section)

rownames(object_cross_section_loess) == rownames(object_cross_section)

transcriptome_cor <-
  purrr::map(1:nrow(object_cross_section_loess), function(i) {
    cat(i, " ")
    temp <-
      cor.test(unlist(object_cross_section_loess[i, , drop = TRUE]),
               unlist(object_cross_section[i, , drop = TRUE]))
    data.frame(cor = unname(temp$estimate),
               p = unname(temp$p.value))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



###proteomics
load(
  "data_analysis/plasma_proteomics/data_preparation/same_samples/object_cross_section_loess"
)

load("data_analysis/plasma_proteomics/data_preparation/object_cross_section")

dim(object_cross_section_loess)
dim(object_cross_section)

rownames(object_cross_section_loess) == rownames(object_cross_section)

proteomics_cor <-
  purrr::map(1:nrow(object_cross_section_loess), function(i) {
    cat(i, " ")
    temp <-
      cor.test(unlist(object_cross_section_loess[i, , drop = TRUE]),
               unlist(object_cross_section[i, , drop = TRUE]))
    data.frame(cor = unname(temp$estimate),
               p = unname(temp$p.value))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()





###proteomics
load(
  "data_analysis/plasma_metabolomics/data_preparation/metabolite/same_samples/object_cross_section_loess"
)

load(
  "data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

dim(object_cross_section_loess)
dim(object_cross_section)

rownames(object_cross_section_loess) == rownames(object_cross_section)

metabolomics_cor <-
  purrr::map(1:nrow(object_cross_section_loess), function(i) {
    cat(i, " ")
    temp <-
      cor.test(unlist(object_cross_section_loess[i, , drop = TRUE]),
               unlist(object_cross_section[i, , drop = TRUE]))
    data.frame(cor = unname(temp$estimate),
               p = unname(temp$p.value))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()
