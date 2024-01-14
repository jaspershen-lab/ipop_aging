no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

load("3-data_analysis/clinical_test/data_preparation/object_corss_section")

setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_2/clinical_test_pathway/"
)

object_cross_section_loess

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

load("../../final_cluster_info")

head(final_cluster_info)

dim(final_cluster_info)

dim(object_cross_section_loess)

object_cross_section_loess@variable_info

variable_info <-
  object_cross_section_loess@variable_info %>%
  dplyr::select(-cluster)

final_cluster_info <-
  final_cluster_info %>%
  dplyr::left_join(variable_info, by = "variable_id")

table(final_cluster_info$cluster, final_cluster_info$class)

library(org.Hs.eg.db)
library(clusterProfiler)

###cluster 2

cluster2_clinical_test <-
  final_cluster_info %>%
  dplyr::filter(cluster == 2 & class == "clinical_test")

dim(cluster2_clinical_test)

cluster2_clinical_test$mol_name
cluster2_clinical_test$variable_id
cluster2_clinical_test$test_name

plot <-
  object_cross_section_loess %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "clinical_test_BUN") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

plot

ggsave(plot,
       file = "clinical_test_BUN.pdf",
       width = 8,
       height = 6)

plot <-
  object_cross_section_loess %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "clinical_test_GLU") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

plot

ggsave(plot,
       file = "clinical_test_GLU.pdf",
       width = 8,
       height = 6)
