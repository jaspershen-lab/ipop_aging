no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("code/tools.R")

library(tidyverse)
library(tidymass)

load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)


dir.create(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/correlation_distributation",
  recursive = TRUE
)

setwd(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/correlation_distributation"
)

load("../cor_data")
load("../final_cluster_info")

variable_info <-
  object_cross_section_loess@variable_info

temp_data <-
  cor_data %>%
  dplyr::left_join(final_cluster_info[, c("variable_id", "cluster")]) %>%
  dplyr::left_join(variable_info[, c("variable_id", "class")])

plot <-
  temp_data %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  ggplot(aes(correlation, cluster)) +
  geom_jitter(aes(alpha = 0.5, color = class)) +
  ggplot2::geom_vline(xintercept = 0, color = "red") +
  geom_boxplot(outlier.shape = NA, fill = "transparent") +
  facet_grid(rows = vars(cluster),
             scales = "free") +
  scale_color_manual(values = omics_color) +
  theme_base +
  theme(panel.grid = element_blank())

plot

ggsave(plot,
       filename = "correlation_distributation.pdf",
       width = 7,
       height = 5)



