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

load(
  "3-data_analysis/plasma_cytokine/data_preparation/object"
)

setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_2/cytokine_pathway/"
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
cluster2_cytokine <-
  final_cluster_info %>%
  dplyr::filter(cluster == 2 & class == "cytokine")
dim(cluster2_cytokine)

cluster2_cytokine$mol_name
cluster2_cytokine$variable_id

object_cross_section_loess %>% 
  ggplot_mass_dataset(direction = "variable", variable_id = "cytokine_ENA78") +
  geom_point(aes(x = age))

object %>% 
  ggplot_mass_dataset(direction = "variable", variable_id = "ENA78") +
  geom_point(aes(x = adjusted_age))







