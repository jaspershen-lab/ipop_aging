no_source()
rm(list = ls())
gc()
library(tidyverse)
library(tidymass)
library(plyr)
setwd(masstools::get_project_wd())

source("code/tools.R")

dir.create("data_analysis/combined_omics/DE_SWAN")
dir.create(
  "data_analysis/combined_omics/DE_SWAN/transcriptome_proteome_metabolome_summary"
)
setwd("data_analysis/combined_omics/DE_SWAN")

load("transcriptome_pathway/result_all")

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

result_all_crest1 <- result_all

####crest2
load("transcriptome_pathway_crest2/result_all")

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

result_all_crest2 <- result_all

result_all_crest1 <-
  result_all_crest1 %>%
  dplyr::arrange(p.adjust)

result_all_crest2 <-
  result_all_crest2 %>%
  dplyr::arrange(p.adjust)

temp_crest1 <-
  result_all_crest1 %>%
  head(20)

temp_crest2 <-
  result_all_crest2 %>%
  head(20)

temp_crest1_transcriptome <-
  temp_crest1

temp_crest2_transcriptome <-
  temp_crest2

load("proteomics_pathway/result_all")

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

result_all_crest1 <- result_all

####crest2
load("proteomics_pathway_crest2/result_all")

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

result_all_crest2 <- result_all

result_all_crest1 <-
  result_all_crest1 %>%
  dplyr::arrange(p.adjust)

result_all_crest2 <-
  result_all_crest2 %>%
  dplyr::arrange(p.adjust)

temp_crest1 <-
  result_all_crest1 %>%
  head(20)

temp_crest2 <-
  result_all_crest2 %>%
  head(20)

temp_crest1_proteome <-
  temp_crest1

temp_crest2_proteome <-
  temp_crest2

load("metabolomics_pathway_crest1/KEGG_result/metabolomics_crest1_kegg")
load("metabolomics_pathway_crest1/HMDB_result/metabolomics_crest1_hmdb")

metabolomics_crest1_kegg <-
  metabolomics_crest1_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::mutate(database = "KEGG")

metabolomics_crest1_hmdb <-
  metabolomics_crest1_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::mutate(database = "HMDB")

result_all_crest1 <-
  rbind(metabolomics_crest1_kegg,
        metabolomics_crest1_hmdb)

load("metabolomics_pathway_crest2/KEGG_result/metabolomics_crest2_kegg")
load("metabolomics_pathway_crest2/HMDB_result/metabolomics_crest2_hmdb")

metabolomics_crest2_kegg <-
  metabolomics_crest2_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::mutate(database = "KEGG")

metabolomics_crest2_hmdb <-
  metabolomics_crest2_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::mutate(database = "HMDB")

result_all_crest2 <-
  rbind(metabolomics_crest2_kegg,
        metabolomics_crest2_hmdb)

temp_crest1 <-
  result_all_crest1 %>%
  head(20)

temp_crest2 <-
  result_all_crest2 %>%
  head(20)

temp_crest1_metabolome <-
  temp_crest1

temp_crest2_metabolome <-
  temp_crest2

######network to show the function of omics data
dim(temp_crest1_transcriptome)
dim(temp_crest1_proteome)
dim(temp_crest1_metabolome)

dim(temp_crest2_transcriptome)
dim(temp_crest2_proteome)
dim(temp_crest2_metabolome)

library(ggraph)
library(tidygraph)

colnames(temp_crest1_transcriptome)

###crest 1
temp_transcriptome_crest1 <-
  temp_crest1_transcriptome %>%
  dplyr::select(module_annotation, p.adjust, database, geneID) %>%
  dplyr::mutate(node = paste(module_annotation, "transcriptome", sep = "_"),
                class = "transcriptome")

temp_proteome_crest1 <-
  temp_crest1_proteome %>%
  dplyr::select(module_annotation, p.adjust, database, geneID) %>%
  dplyr::mutate(node = paste(module_annotation, "proteome", sep = "_"),
                class = "proteomics")

temp_metabolome_crest1 <-
  temp_crest1_metabolome %>%
  dplyr::select(pathway_name, p_value_adjust, database, mapped_id) %>%
  dplyr::rename(module_annotation = pathway_name,
                p.adjust = p_value_adjust,
                geneID = mapped_id) %>%
  dplyr::mutate(node = paste(module_annotation, "metabolome", sep = "_"),
                class = "metabolomics")

node_data <-
  rbind(temp_transcriptome_crest1,
        temp_proteome_crest1,
        temp_metabolome_crest1) %>% 
  dplyr::select(node, everything()) 

###calculate the jaccard index between nodes
edge_data <-
  seq_len(nrow(node_data) - 1) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <-
      node_data$geneID[i] %>% stringr::str_split(pattern = "/") %>% `[[`(1)
    
    (i + 1):nrow(node_data) %>%
      purrr::map(function(j) {
        y <-
          node_data$geneID[j] %>% stringr::str_split(pattern = "/") %>% `[[`(1)
        data.frame(
          from = node_data$node[i],
          to = node_data$node[j],
          jaccard_index = calculate_jaccard_index(x, y)
        )
      }) %>%
      dplyr::bind_rows()
  }) %>%
  dplyr::bind_rows()

edge_data <-
  edge_data %>%
  dplyr::filter(jaccard_index > 0.1)

netwrok <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data)

plot =
  ggraph(netwrok,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(color = "grey",
                 aes(width = jaccard_index),
                 show.legend = TRUE) +
  geom_node_point(aes(fill = class,
                      size = -log(p.adjust, 10)),
                  shape = 21,
                  show.legend = TRUE) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = module_annotation
    ),
    size = 3
  ) +
  # shadowtext::geom_shadowtext(
  #   aes(
  #     x = x,
  #     y = y,
  #     label = module_annotation
  #   ),
  #   color = "black",
  #   bg.color = "white",
  #   size = 5.5,
  #   show.legend = FALSE
  # ) +
  scale_fill_manual(values = c("transcriptome" = unname(omics_color["transcriptome"]),
                                "proteomics" = unname(omics_color["proteomics"]),
                                "metabolomics" = unname(omics_color["metabolomics"]))) +
  scale_size_continuous(range = c(2,8)) +
  scale_edge_width_continuous(range = c(0.5,4)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# library(extrafont)
# extrafont::loadfonts()
# ggsave(plot,
#        filename = "transcriptome_proteome_metabolome_summary/crest1_network.pdf",
#        width = 7,
#        height = 6)


###crest 2
temp_transcriptome_crest2 <-
  temp_crest2_transcriptome %>%
  dplyr::select(module_annotation, p.adjust, database, geneID) %>%
  dplyr::mutate(node = paste(module_annotation, "transcriptome", sep = "_"),
                class = "transcriptome")

temp_proteome_crest2 <-
  temp_crest2_proteome %>%
  dplyr::select(module_annotation, p.adjust, database, geneID) %>%
  dplyr::mutate(node = paste(module_annotation, "proteome", sep = "_"),
                class = "proteomics")

temp_metabolome_crest2 <-
  temp_crest2_metabolome %>%
  dplyr::select(pathway_name, p_value_adjust, database, mapped_id) %>%
  dplyr::rename(module_annotation = pathway_name,
                p.adjust = p_value_adjust,
                geneID = mapped_id) %>%
  dplyr::mutate(node = paste(module_annotation, "metabolome", sep = "_"),
                class = "metabolomics")

node_data <-
  rbind(temp_transcriptome_crest2,
        temp_proteome_crest2,
        temp_metabolome_crest2) %>% 
  dplyr::select(node, everything()) 

###calculate the jaccard index between nodes
edge_data <-
  seq_len(nrow(node_data) - 1) %>%
  purrr::map(function(i) {
    cat(i, " ")
    x <-
      node_data$geneID[i] %>% stringr::str_split(pattern = "/") %>% `[[`(1)
    
    (i + 1):nrow(node_data) %>%
      purrr::map(function(j) {
        y <-
          node_data$geneID[j] %>% stringr::str_split(pattern = "/") %>% `[[`(1)
        data.frame(
          from = node_data$node[i],
          to = node_data$node[j],
          jaccard_index = calculate_jaccard_index(x, y)
        )
      }) %>%
      dplyr::bind_rows()
  }) %>%
  dplyr::bind_rows()

edge_data <-
  edge_data %>%
  dplyr::filter(jaccard_index > 0.1)

netwrok <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data)

plot =
  ggraph(netwrok,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(color = "grey",
                 aes(width = jaccard_index),
                 show.legend = TRUE) +
  geom_node_point(aes(fill = class,
                      size = -log(p.adjust, 10)),
                  shape = 21,
                  show.legend = TRUE) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = module_annotation
    ),
    size = 3
  ) +
  # shadowtext::geom_shadowtext(
  #   aes(
  #     x = x,
  #     y = y,
  #     label = module_annotation
  #   ),
  #   color = "black",
  #   bg.color = "white",
  #   size = 5.5,
  #   show.legend = FALSE
  # ) +
scale_fill_manual(values = c("transcriptome" = unname(omics_color["transcriptome"]),
                             "proteomics" = unname(omics_color["proteomics"]),
                             "metabolomics" = unname(omics_color["metabolomics"]))) +
  scale_size_continuous(range = c(2,8)) +
  scale_edge_width_continuous(range = c(0.5,4)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# library(extrafont)
# extrafont::loadfonts()
# ggsave(plot,
#        filename = "transcriptome_proteome_metabolome_summary/crest2_network.pdf",
#        width = 7,
#        height = 6)
