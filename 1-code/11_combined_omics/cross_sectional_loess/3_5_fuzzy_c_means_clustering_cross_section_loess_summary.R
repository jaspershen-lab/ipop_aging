no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

dir.create(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess_summary/",
  recursive = TRUE
)

setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess_summary"
)

load("../cross_section_loess/final_cluster_info")
cluster_info_all <-
  final_cluster_info

load("../cross_section_loess_only_transcriptomics/final_cluster_info")
cluster_info_transcriptomics <-
  final_cluster_info

load("../cross_section_loess_except_transcriptomics/final_cluster_info")
cluster_info_not_transcriptomics <-
  final_cluster_info

dim(cluster_info_all)
dim(cluster_info_transcriptomics)
dim(cluster_info_not_transcriptomics)

######Only transcriptomics
###cluster 2
cluster2_id <-
  cluster_info_all %>%
  dplyr::filter(cluster == 2)

library(plyr)
cluster_info_transcriptomics %>%
  plyr::dlply(.variables = .(cluster)) %>%
  purrr::map(function(x) {
    intersect_id <-
      intersect(x$variable_id,
                cluster2_id$variable_id)
    data.frame(
      cluster = x$cluster[1],
      length = length(intersect_id),
      percentage = length(intersect_id) / sum(stringr::str_detect(cluster2_id$variable_id, "transcriptome"))
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####2 vs 11



###cluster 4
cluster4_id <-
  cluster_info_all %>%
  dplyr::filter(cluster == 4)

library(plyr)
cluster_info_transcriptomics %>%
  plyr::dlply(.variables = .(cluster)) %>%
  purrr::map(function(x) {
    intersect_id <-
      intersect(x$variable_id,
                cluster4_id$variable_id)
    data.frame(
      cluster = x$cluster[1],
      length = length(intersect_id),
      percentage = length(intersect_id) / sum(stringr::str_detect(cluster4_id$variable_id, "transcriptome"))
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####4 vs 6




###cluster 5
cluster5_id <-
  cluster_info_all %>%
  dplyr::filter(cluster == 5)

library(plyr)
cluster_info_transcriptomics %>%
  plyr::dlply(.variables = .(cluster)) %>%
  purrr::map(function(x) {
    intersect_id <-
      intersect(x$variable_id,
                cluster5_id$variable_id)
    data.frame(
      cluster = x$cluster[1],
      length = length(intersect_id),
      percentage = length(intersect_id) / sum(stringr::str_detect(cluster5_id$variable_id, "transcriptome"))
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####5 vs 4 and 5
















######except transcriptomics
###cluster 2
cluster2_id <-
  cluster_info_all %>%
  dplyr::filter(cluster == 2)

library(plyr)
cluster_info_not_transcriptomics %>%
  plyr::dlply(.variables = .(cluster)) %>%
  purrr::map(function(x) {
    intersect_id <-
      intersect(x$variable_id,
                cluster2_id$variable_id)
    data.frame(
      cluster = x$cluster[1],
      length = length(intersect_id),
      percentage = length(intersect_id) / sum(!stringr::str_detect(cluster2_id$variable_id, "transcriptome"))
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####2 vs 2 and 11



###cluster 4
cluster4_id <-
  cluster_info_all %>%
  dplyr::filter(cluster == 4)

library(plyr)
cluster_info_not_transcriptomics %>%
  plyr::dlply(.variables = .(cluster)) %>%
  purrr::map(function(x) {
    intersect_id <-
      intersect(x$variable_id,
                cluster4_id$variable_id)
    data.frame(
      cluster = x$cluster[1],
      length = length(intersect_id),
      percentage = length(intersect_id) / sum(!stringr::str_detect(cluster4_id$variable_id, "transcriptome"))
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####4 vs 3 and 4




###cluster 5
cluster5_id <-
  cluster_info_all %>%
  dplyr::filter(cluster == 5)

library(plyr)
cluster_info_not_transcriptomics %>%
  plyr::dlply(.variables = .(cluster)) %>%
  purrr::map(function(x) {
    intersect_id <-
      intersect(x$variable_id,
                cluster5_id$variable_id)
    data.frame(
      cluster = x$cluster[1],
      length = length(intersect_id),
      percentage = length(intersect_id) / sum(!stringr::str_detect(cluster5_id$variable_id, "transcriptome"))
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####5 vs 2
