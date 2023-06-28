no_soure()

setwd(masstools::get_project_wd())
rm(list = ls())

library(ggraph)
library(igraph)
library(tidygraph)
library(tidyverse)
library(tidymass)

source("code/tools.R")

load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)
load(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

# load("data_analysis/combined_omics/correlation_data/cor_data")

dir.create(
  "data_analysis/combined_omics/correlation_network/cross_section_loess",
  recursive = TRUE
)

setwd("data_analysis/combined_omics/correlation_network/cross_section_loess")

table(final_cluster_info$cluster)

variable_info <-
  object_cross_section_loess@variable_info

#network for each cluster
for (cluster_idx in unique(final_cluster_info$cluster)) {
  cat(cluster_idx, " ")
  
  cluster_info <-
    final_cluster_info %>%
    dplyr::filter(cluster == cluster_idx)
  
  temp_data <-
    object_cross_section_loess@expression_data[cluster_info$variable_id, ]
 
  # Distance matrix
  d <- dist(temp_data, method = "euclidean")
  
  # Hierarchical clustering
  hc <- hclust(d, method = "complete")
  hc_cluster <- cutree(hc, k = 20)
  table(hc_cluster)
  
}



for (i in unique(final_cluster_info$cluster)) {
  cat(i, " ")
  temp <-
    final_cluster_info %>%
    dplyr::filter(cluster == i)
  
  temp_object <-
    object_cross_section_loess[temp$variable_id,]
  
  dir.create(paste0("cluster_", i))
  
  
  
}
