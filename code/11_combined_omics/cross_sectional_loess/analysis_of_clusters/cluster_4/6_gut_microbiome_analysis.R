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
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_4/gut_microbiome_pathway/"
)

setwd(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_4/gut_microbiome_pathway/"
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

###cluster 4
cluster4_gut_microbiome <-
  final_cluster_info %>%
  dplyr::filter(cluster == 4 & class == "gut_microbiome")
dim(cluster4_gut_microbiome)

cluster4_gut_microbiome <- 
cluster4_gut_microbiome[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

write.csv(cluster4_gut_microbiome, file = "cluster4_gut_microbiome.csv", row.names = FALSE)

cluster4_gut_microbiome[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]$Genus %>% 
  unique()





