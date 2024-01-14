no_source()

library(masstools)
setwd_project()
rm(list = ls())
source("1-code/100-tools.R")

###transcriptome
load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("3-data_analysis/combined_omics/correlation_data")
setwd("3-data_analysis/combined_omics/correlation_data")

cor_data <-
  cor_mass_dataset(
    x = object_cross_section_loess,
    margin = "variable",
    method = "spearman",
    data_type = "longer",
    p_adjust_method = "BH"
  )

cor_data <-
  cor_data %>%
  dplyr::filter(p_adjust < 0.05)

cor_data <-
  cor_data %>%
  dplyr::filter(abs(correlation) > 0.5)

cor_data <-
  cor_data %>%
  dplyr::mutate(
    from_class =
      stringr::str_split(from, "_", n = 2) %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist,
    to_class =
      stringr::str_split(to, "_", n = 2) %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist
  )

cor_data$from_class[cor_data$from_class %in% c("gut", "skin", "oral", "nasal")] <- 
  paste0(cor_data$from_class[cor_data$from_class %in% c("gut", "skin", "oral", "nasal")],
        "_microbiome")

cor_data$to_class[cor_data$to_class %in% c("gut", "skin", "oral", "nasal")] <- 
  paste0(cor_data$to_class[cor_data$to_class %in% c("gut", "skin", "oral", "nasal")],
         "_microbiome")

save(cor_data, file = "cor_data")
