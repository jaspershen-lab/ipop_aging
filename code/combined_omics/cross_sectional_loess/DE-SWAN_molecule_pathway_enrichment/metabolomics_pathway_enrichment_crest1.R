no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
masstools::setwd_project()

source("code/tools.R")

load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("data_analysis/combined_omics/DE_SWAN")
setwd("data_analysis/combined_omics/DE_SWAN")

library("DEswan")

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id))

expression_data <-
  object_cross_section_loess@expression_data

sample_info <-
  object_cross_section_loess@sample_info

variable_info <-
  object_cross_section_loess@variable_info

load("temp_data")

temp_data <-
  temp_data %>%
  dplyr::left_join(variable_info,
                   by = c("variable_id"))

###Pathway enrichment
temp_data_metabolomics <-
  temp_data %>%
  dplyr::filter(class == "metabolomics" & p_value_adjust < 0.05)

####metabolomics
dir.create("metabolomics_pathway_crest1")
dir.create("metabolomics_pathway_crest1/HMDB_result/")
dir.create("metabolomics_pathway_crest1/KEGG_result/")

####crest 1
temp_data_metabolomics_crest1 <-
  temp_data_metabolomics %>%
  dplyr::filter(center == 47)

library(metpath)

####HMDB pathway
data("hmdb_pathway", package = "metpath")
hmdb_pathway
get_pathway_class(hmdb_pathway)
#get the class of pathways
pathway_class =
  metpath::pathway_class(hmdb_pathway)

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

hmdb_pathway =
  hmdb_pathway[remain_idx]

hmdb_pathway

###HMDB enrichment
query_id <-
  temp_data_metabolomics_crest1$HMDB.ID[!is.na(temp_data_metabolomics_crest1$HMDB.ID)] %>%
  stringr::str_split("\\{\\}") %>%
  unlist() %>%
  unique()

# metabolomics_crest1_hmdb <-
#   enrich_hmdb(
#     query_id = query_id,
#     query_type = "compound",
#     id_type = "HMDB",
#     pathway_database = hmdb_pathway,
#     only_primary_pathway = TRUE,
#     p_cutoff = 0.05,
#     p_adjust_method = "fdr",
#     threads = 3
#   )
#
# metabolomics_crest1_hmdb <-
#   metabolomics_crest1_hmdb %>%
#   dplyr::filter(p_value < 0.05)
#
# save(metabolomics_crest1_hmdb, file = "metabolomics_pathway_crest1/HMDB_result/metabolomics_crest1_hmdb")
load("metabolomics_pathway_crest1/HMDB_result/metabolomics_crest1_hmdb")


###KEGG enrichment
data("kegg_hsa_pathway", package = "metpath")
kegg_hsa_pathway
get_pathway_class(kegg_hsa_pathway)

#get the class of pathways
pathway_class =
  metpath::pathway_class(kegg_hsa_pathway)

remain_idx =
  pathway_class %>%
  unlist() %>%
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

kegg_hsa_pathway =
  kegg_hsa_pathway[remain_idx]

kegg_hsa_pathway

query_id <-
  temp_data_metabolomics_crest1$KEGG.ID[!is.na(temp_data_metabolomics_crest1$KEGG.ID)] %>%
  stringr::str_split("\\{\\}") %>%
  unlist() %>%
  unique()

# metabolomics_crest1_kegg <-
#   enrich_kegg(
#     query_id = query_id,
#     query_type = "compound",
#     id_type = "KEGG",
#     pathway_database = kegg_hsa_pathway,
#     p_cutoff = 0.05,
#     p_adjust_method = "fdr",
#     threads = 3
#   )
#
# metabolomics_crest1_kegg <-
#   metabolomics_crest1_kegg %>%
#   dplyr::filter(p_value < 0.05)
#
# save(metabolomics_crest1_kegg, file = "metabolomics_pathway_crest1/kegg_result/metabolomics_crest1_kegg")

load("metabolomics_pathway_crest1/kegg_result/metabolomics_crest1_kegg")
load("metabolomics_pathway_crest1/HMDB_result/metabolomics_crest1_hmdb")

metabolomics_crest1_kegg <-
  metabolomics_crest1_kegg %>%
  dplyr::filter(p_value < 0.05 & mapped_number >= 3) %>%
  dplyr::arrange(p_value)

###KEGG
load("metabolomics_pathway_crest1/KEGG_result/metabolomics_crest1_kegg")
metabolomics_crest1_hmdb <-
  metabolomics_crest1_hmdb %>%
  dplyr::filter(p_value < 0.05 & mapped_number >= 3) %>%
  dplyr::arrange(p_value)

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "KEGG")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = metabolomics_crest1_kegg@result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(
  wb,
  "metabolomics_pathway_crest1/KEGG_result/metabolomics_crest1_kegg.xlsx",
  overwrite = TRUE
)


library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "HMDB")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = metabolomics_crest1_hmdb@result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(
  wb,
  "metabolomics_pathway_crest1/HMDB_result/metabolomics_crest1_hmdb.xlsx",
  overwrite = TRUE
)

plot <- 
metpath::enrich_bar_plot(
  metabolomics_crest1_kegg,
  x_axis = "p_value_adjust",
  cutoff = 0.05,
  top = 10
)

ggsave(plot, filename = "metabolomics_pathway_crest1/KEGG_result/enrichment_pathway.pdf",
       width = 7, height = 7)

plot <-
metpath::enrich_scatter_plot(
  metabolomics_crest1_kegg,
  y_axis = "p_value_adjust",
  y_axis_cutoff = 0.05
)

ggsave(plot, filename = "metabolomics_pathway_crest1/KEGG_result/enrichment_pathway2.pdf",
       width = 7, height = 7)



