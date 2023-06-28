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

setwd(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_5/metabolomics_pathway/"
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

###cluster 5
cluster5_metabolomics <-
  final_cluster_info %>%
  dplyr::filter(cluster == 5 & class == "metabolomics")

save(cluster5_metabolomics, file = "cluster5_metabolomics")

write.csv(cluster5_metabolomics,
          "cluster5_metabolomics.csv",
          row.names = FALSE)

cluster5_metabolomics$KEGG.ID

load("KEGG_result/metabolomics_kegg")
load("HMDB_result/metabolomics_hmdb")

metabolomics_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05)

metabolomics_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05)

cluster5_pathway <-
  rbind(
    data.frame(metabolomics_kegg@result, class = "KEGG"),
    data.frame(metabolomics_hmdb@result, class = "HMDB")
  ) %>%
  dplyr::filter(p_value_adjust < 0.05)

save(cluster5_pathway, file = "cluster5_pathway")

write.csv(cluster5_pathway, file = "cluster5_pathway.csv", row.names = FALSE)
