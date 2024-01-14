no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "3-data_analysis/plasma_cytokine/data_preparation/object_cross_section"
)
load(
  "3-data_analysis/plasma_cytokine/data_preparation/object_cross_section_loess"
)

dir.create(
  "3-data_analysis/plasma_cytokine/correlation_with_age/cross_section/",
  recursive = TRUE
)

dir.create(
  "3-data_analysis/plasma_cytokine/correlation_with_age/cross_section_loess/",
  recursive = TRUE
)

setwd("3-data_analysis/plasma_cytokine/correlation_with_age")

object_cross_section
object_cross_section_loess

dim(object_cross_section)

###calculate the correlation between age and molecules
cor_data <-
  seq_len(nrow(object_cross_section)) %>%
  purrr::map(function(i) {
    value <-
      object_cross_section[i, , drop = TRUE] %>%
      unlist()
    age <-
      object_cross_section@sample_info$adjusted_age
    test <-
      cor.test(age, value, method = "spearman")
    data.frame(
      variable_id = rownames(object_cross_section)[i],
      correlation = unname(test$estimate),
      p_value = test$p.value
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

save(cor_data, file = "cross_section/cor_data")

###calculate the correlation between age and molecules

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id))

cor_data_loess <-
  seq_len(nrow(object_cross_section_loess)) %>%
  purrr::map(function(i) {
    value <-
      object_cross_section_loess[i, , drop = TRUE] %>%
      unlist()
    age <-
      object_cross_section_loess@sample_info$age
    test <-
      cor.test(age, value, method = "spearman")
    data.frame(
      variable_id = rownames(object_cross_section_loess)[i],
      correlation = unname(test$estimate),
      p_value = test$p.value
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

save(cor_data_loess, file = "cross_section_loess/cor_data_loess")

plot(cor_data$correlation, cor_data_loess$correlation)
cor(cor_data$correlation, cor_data_loess$correlation)
