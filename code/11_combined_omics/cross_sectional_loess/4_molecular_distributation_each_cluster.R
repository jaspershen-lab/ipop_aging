no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("code/tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/",
  recursive = TRUE
)

setwd("data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess")

object_cross_section_loess

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

load("cor_data")

load("final_cluster_info")

head(final_cluster_info)

dim(final_cluster_info)

dim(object_cross_section_loess)

object_cross_section_loess@variable_info

####mosaic plot
library(ggmosaic)

variable_info <-
  object_cross_section_loess@variable_info[, c("variable_id", "class")]

final_cluster_info <-
  final_cluster_info %>%
  dplyr::left_join(variable_info, by = "variable_id")

plot <-
  final_cluster_info %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::mutate(cluster = factor(cluster, levels = stringr::str_sort(unique(cluster), numeric = TRUE))) %>%
  dplyr::mutate(class = factor(class, levels = names(omics_color))) %>%
  ggplot() +
  geom_mosaic(aes(x = product(class, cluster), fill = class),
              offset = 0.01) +
  theme_mosaic() +
  scale_fill_manual(values = omics_color) +
  labs(x = "", y = "") +
  theme(
    panel.border = element_rect(color = "black",
                                fill = "transparent"),
    legend.position = "bottom"
  )

plot

# ggsave(plot, filename = "cluster_omics_distributation.pdf", width = 7, height = 7)
