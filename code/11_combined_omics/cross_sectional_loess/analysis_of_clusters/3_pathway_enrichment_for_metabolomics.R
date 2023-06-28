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

1:11 %>%
  purrr::walk(function(i) {
    message(i)
    name <-
      paste0(
        "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_",
        i
      )
    dir.create(name,
               recursive = TRUE)
  })

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

variable_info <-
  object_cross_section_loess@variable_info %>%
  dplyr::select(-cluster)

final_cluster_info <-
  final_cluster_info %>%
  dplyr::left_join(variable_info, by = "variable_id")

library(org.Hs.eg.db)
library(clusterProfiler)

# ###metabolomics
# for (i in 1:11) {
#   path <- paste0("cluster_", i)
#   dir.create(file.path(path, "metabolomics_pathway"))
#   dir.create(file.path(path, "metabolomics_pathway/KEGG_result/"))
#   dir.create(file.path(path, "metabolomics_pathway/HMDB_result/"))
#   message(i)
#   cluster_info <-
#     final_cluster_info %>%
#     dplyr::filter(cluster == i)
# 
#   ####metabolomics
#   temp_data_metabolomics <-
#     cluster_info %>%
#     dplyr::filter(class == "metabolomics")
# 
#   if (nrow(temp_data_metabolomics) == 0) {
#     next()
#   }
# 
#   ####HMDB pathway
#   data("hmdb_pathway", package = "metpath")
#   hmdb_pathway
#   get_pathway_class(hmdb_pathway)
#   #get the class of pathways
#   pathway_class =
#     metpath::pathway_class(hmdb_pathway)
# 
#   remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")
# 
#   hmdb_pathway =
#     hmdb_pathway[remain_idx]
# 
#   hmdb_pathway
# 
#   ###HMDB enrichment
#   query_id <-
#     temp_data_metabolomics$HMDB.ID[!is.na(temp_data_metabolomics$HMDB.ID)] %>%
#     stringr::str_split("\\{\\}") %>%
#     unlist() %>%
#     unique()
# 
#   if (all(dir(file.path(
#     path,
#     "metabolomics_pathway/HMDB_result/"
#   )) != "metabolomics_hmdb")) {
#     metabolomics_hmdb <-
#       enrich_hmdb(
#         query_id = query_id,
#         query_type = "compound",
#         id_type = "HMDB",
#         pathway_database = hmdb_pathway,
#         only_primary_pathway = TRUE,
#         p_cutoff = 0.05,
#         p_adjust_method = "fdr",
#         threads = 3
#       )
# 
#     save(
#       metabolomics_hmdb,
#       file = file.path(
#         path,
#         "metabolomics_pathway/HMDB_result/metabolomics_hmdb"
#       )
#     )
#   } else{
#     load(file.path(
#       path,
#       "metabolomics_pathway/HMDB_result/metabolomics_hmdb"
#     ))
#   }
# 
#   if (!is.null(metabolomics_hmdb)) {
#     metabolomics_hmdb <-
#       metabolomics_hmdb %>%
#       dplyr::filter(p_value < 0.05 & mapped_number >= 3) %>%
#       dplyr::arrange(p_value)
#   }
# 
#   ###KEGG enrichment
#   data("kegg_hsa_pathway", package = "metpath")
#   kegg_hsa_pathway
#   get_pathway_class(kegg_hsa_pathway)
# 
#   #get the class of pathways
#   pathway_class =
#     metpath::pathway_class(kegg_hsa_pathway)
# 
#   remain_idx =
#     pathway_class %>%
#     unlist() %>%
#     stringr::str_detect("Disease") %>%
#     `!`() %>%
#     which()
# 
#   kegg_hsa_pathway =
#     kegg_hsa_pathway[remain_idx]
# 
#   kegg_hsa_pathway
# 
#   query_id <-
#     temp_data_metabolomics$KEGG.ID[!is.na(temp_data_metabolomics$KEGG.ID)] %>%
#     stringr::str_split("\\{\\}") %>%
#     unlist() %>%
#     unique()
# 
#   if (all(dir(file.path(
#     path,
#     "metabolomics_pathway/KEGG_result"
#   )) != "metabolomics_kegg")) {
#     metabolomics_kegg <-
#       enrich_kegg(
#         query_id = query_id,
#         query_type = "compound",
#         id_type = "KEGG",
#         pathway_database = kegg_hsa_pathway,
#         p_cutoff = 0.2,
#         p_adjust_method = "fdr",
#         threads = 3
#       )
# 
#     save(
#       metabolomics_kegg,
#       file = file.path(
#         path,
#         "metabolomics_pathway/kegg_result/metabolomics_kegg"
#       )
#     )
# 
#   } else{
#     load(file.path(
#       path,
#       "metabolomics_pathway/kegg_result/metabolomics_kegg"
#     ))
#   }
# 
#   if (!is.null(metabolomics_kegg)) {
#     metabolomics_kegg <-
#       metabolomics_kegg %>%
#       dplyr::filter(p_value < 0.2 & mapped_number >= 3) %>%
#       dplyr::arrange(p_value)
#   }
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "KEGG")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   if (!is.null(metabolomics_kegg)) {
#     x <- metabolomics_kegg@result
#   } else{
#     x <-
#       data.frame(
#         "pathway_id" = character(),
#         "pathway_name" = character(),
#         "describtion" = character(),
#         "pathway_class" = character(),
#         "p_value" = numeric(),
#         "p_value_adjust" = numeric(),
#         "all_id" = character(),
#         "all_number" = numeric(),
#         "mapped_id" = character(),
#         "mapped_number" = numeric(),
#         "mapped_percentage" = numeric()
#       )
#   }
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = x,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(
#     wb,
#     file = file.path(
#       path,
#       "metabolomics_pathway/KEGG_result/metabolomics_kegg.xlsx"
#     ),
#     overwrite = TRUE
#   )
# 
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "HMDB")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   if (!is.null(metabolomics_hmdb)) {
#     x <- metabolomics_hmdb@result
#   } else{
#     x <-
#       data.frame(
#         "pathway_id" = character(),
#         "pathway_name" = character(),
#         "describtion" = character(),
#         "pathway_class" = character(),
#         "p_value" = numeric(),
#         "p_value_adjust" = numeric(),
#         "all_id" = character(),
#         "all_number" = numeric(),
#         "mapped_id" = character(),
#         "mapped_number" = numeric(),
#         "mapped_percentage" = numeric()
#       )
#   }
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = x,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(
#     wb,
#     file.path(
#       path,
#       "metabolomics_pathway/HMDB_result/metabolomics_hmdb.xlsx"
#     ),
#     overwrite = TRUE
#   )
# 
#   if (!is.null(metabolomics_kegg)) {
#     plot <-
#       metpath::enrich_bar_plot(
#         metabolomics_kegg,
#         x_axis = "p_value_adjust",
#         cutoff = 0.05,
#         top = 10
#       )
# 
#     if (is.null(plot)) {
#       plot <-
#         patchwork::plot_spacer()
#     }
# 
#   } else{
#     plot <-
#       patchwork::plot_spacer()
#   }
# 
# 
#   ggsave(
#     plot,
#     filename = file.path(
#       path,
#       "metabolomics_pathway/KEGG_result/enrichment_pathway.pdf"
#     ),
#     width = 7,
#     height = 7
#   )
# 
#   if (!is.null(metabolomics_kegg)) {
#     plot <-
#       metpath::enrich_scatter_plot(metabolomics_kegg,
#                                    y_axis = "p_value_adjust",
#                                    y_axis_cutoff = 0.05)
# 
#     if (is.null(plot)) {
#       plot <-
#         patchwork::plot_spacer()
#     }
# 
#   } else{
#     plot <-
#       patchwork::plot_spacer()
#   }
# 
# 
#   ggsave(
#     plot,
#     filename = file.path(
#       path,
#       "metabolomics_pathway/KEGG_result/enrichment_pathway2.pdf"
#     ),
#     width = 7,
#     height = 7
#   )
# }
