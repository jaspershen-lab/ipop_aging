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

setwd("3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess")

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

variable_info <-
  object_cross_section_loess@variable_info %>%
  dplyr::select(-cluster)

variable_info <-
  variable_info %>%
  dplyr::left_join(cor_data, by = "variable_id")

sample_info <-
  object_cross_section_loess@sample_info

expression_data <-
  object_cross_section_loess@expression_data

final_cluster_info <-
  final_cluster_info %>%
  dplyr::left_join(variable_info, by = "variable_id")

library(org.Hs.eg.db)
library(clusterProfiler)

###Cluster 1
result_all_cluster1 <-
  tryCatch(
    expr = {
      load('cluster_1/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )

load('cluster_1/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster1 <- proteomics_go
load('cluster_1/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster1 <- proteomics_kegg
load('cluster_1/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster1 <- proteomics_reactome

result_all_cluster1$pathway_id

###Cluster 2
result_all_cluster2 <-
  tryCatch(
    expr = {
      load('cluster_2/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )

load('cluster_2/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster2 <- proteomics_go
load('cluster_2/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster2 <- proteomics_kegg
load('cluster_2/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster2 <- proteomics_reactome

result_all_cluster2$pathway_id

###Cluster 3
result_all_cluster3 <-
  tryCatch(
    expr = {
      load('cluster_3/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_3/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster3 <- proteomics_go
load('cluster_3/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster3 <- proteomics_kegg
load('cluster_3/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster3 <- proteomics_reactome

result_all_cluster3$pathway_id

###Cluster 4
result_all_cluster4 <-
  tryCatch(
    expr = {
      load('cluster_4/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_4/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster4 <- proteomics_go
load('cluster_4/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster4 <- proteomics_kegg
load('cluster_4/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster4 <- proteomics_reactome

result_all_cluster4$pathway_id

###Cluster 5
result_all_cluster5 <-
  tryCatch(
    expr = {
      load('cluster_5/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_5/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster5 <- proteomics_go
load('cluster_5/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster5 <- proteomics_kegg
load('cluster_5/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster5 <- proteomics_reactome

result_all_cluster5$pathway_id

###Cluster 6
result_all_cluster6 <-
  tryCatch(
    expr = {
      load('cluster_6/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_6/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster6 <- proteomics_go
load('cluster_6/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster6 <- proteomics_kegg
load('cluster_6/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster6 <- proteomics_reactome

result_all_cluster6$pathway_id

###Cluster 7
result_all_cluster7 <-
  tryCatch(
    expr = {
      load('cluster_7/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_7/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster7 <- proteomics_go
load('cluster_7/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster7 <- proteomics_kegg
load('cluster_7/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster7 <- proteomics_reactome

result_all_cluster7$pathway_id

###Cluster 8
result_all_cluster8 <-
  tryCatch(
    expr = {
      load('cluster_8/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_8/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster8 <- proteomics_go
load('cluster_8/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster8 <- proteomics_kegg
load('cluster_8/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster8 <- proteomics_reactome

result_all_cluster8$pathway_id

###Cluster 9
result_all_cluster9 <-
  tryCatch(
    expr = {
      load('cluster_9/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_9/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster9 <- proteomics_go
load('cluster_9/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster9 <- proteomics_kegg
load('cluster_9/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster9 <- proteomics_reactome

result_all_cluster9$pathway_id

###Cluster 10
result_all_cluster10 <-
  tryCatch(
    expr = {
      load('cluster_10/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_10/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster10 <- proteomics_go
load('cluster_10/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster10 <- proteomics_kegg
load('cluster_10/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster10 <- proteomics_reactome

result_all_cluster10$pathway_id

###Cluster 11
result_all_cluster11 <-
  tryCatch(
    expr = {
      load('cluster_11/proteomics_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_11/proteomics_pathway/GO_result/proteomics_go')
proteomics_go_cluster11 <- proteomics_go
load('cluster_11/proteomics_pathway/KEGG_result/proteomics_kegg')
proteomics_kegg_cluster11 <- proteomics_kegg
load('cluster_11/proteomics_pathway/Reactome_result/proteomics_reactome')
proteomics_reactome_cluster11 <- proteomics_reactome

result_all_cluster11$pathway_id

###only Cluter 8 and 11 have enrichmed pathways
result_all_cluster1$Description
result_all_cluster2$Description
result_all_cluster3$Description
result_all_cluster4$Description
result_all_cluster5$Description
result_all_cluster6$Description
result_all_cluster7$Description
result_all_cluster8$Description
result_all_cluster9$Description
result_all_cluster10$Description
result_all_cluster11$Description

result_all_list <-
  list(
    cluster11 = result_all_cluster11 %>% dplyr::arrange(p.adjust),
    cluster8 = result_all_cluster8 %>% dplyr::arrange(p.adjust)
  )

word_cloud_list <-
  list(cluster11 = result_all_cluster11,
       cluster8 = result_all_cluster8) %>%
  purrr::map(function(x) {
    y <-
      stringr::str_split(x$Description, ";") %>%
      unlist() %>%
      stringr::str_split(pattern = " ") %>%
      unlist() %>%
      # stringr::str_to_lower() %>%
      stringr::str_replace_all("\\(", "") %>%
      stringr::str_replace_all("\\)", "") %>%
      stringr::str_replace_all("\\,", "") %>%
      stringr::str_replace_all("\\.", "") %>%
      stringr::str_replace_all("\\:", "")
    y <-
      combine_workds(y = y)
    y <- y[!y %in% remove_words] %>%
      sort()
    combine_workds(y = y)
  })

plot <-
  pahtway_heatmap_all_cluster(
    result_all_list = result_all_list,
    word_cloud_list = word_cloud_list,
    gene_marker = variable_info,
    expression_data = expression_data,
    top_pathway = 10,
    font_size = c(5, 10),
    word_count_breaks = c(1, 30),
    add_text = TRUE,
    show_column_names = FALSE
  )

plot

ggsave(plot,
       filename = "proteomics_pathway_heatmap.pdf",
       width = 7,
       height = 4)
