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

setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_10/metabolomics_pathway/"
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

###cluster 10
cluster10_metabolomics <-
  final_cluster_info %>%
  dplyr::filter(cluster == 10 & class == "metabolomics")

write.csv(cluster10_metabolomics,
          "cluster10_metabolomics.csv",
          row.names = FALSE)

save(cluster10_metabolomics, file = "cluster10_metabolomics")

cluster10_metabolomics$KEGG.ID

load("KEGG_result/metabolomics_kegg")
load("HMDB_result/metabolomics_hmdb")

metabolomics_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05)

metabolomics_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05)

cluster10_pathway <-
  rbind(
    data.frame(metabolomics_kegg@result, class = "KEGG"),
    data.frame(metabolomics_hmdb@result, class = "HMDB")
  ) %>%
  dplyr::filter(p_value_adjust < 0.05)

save(cluster10_pathway, file = "cluster10_pathway")

write.csv(cluster10_pathway, file = "cluster10_pathway.csv", row.names = FALSE)

stringr::str_split(cluster10_pathway$mapped_id, ";")[[1]] %>% 
  purrr::map(function(x){
    masstools::trans_ID(query = x, 
                        from = "KEGG", to = "Chemical Name")    
  })


