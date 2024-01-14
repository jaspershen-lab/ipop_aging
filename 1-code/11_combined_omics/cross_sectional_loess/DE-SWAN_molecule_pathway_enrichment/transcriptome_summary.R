no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
library(plyr)
setwd(masstools::get_project_wd())

source("1-code/100-tools.R")

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("3-data_analysis/combined_omics/DE_SWAN")
dir.create("3-data_analysis/combined_omics/DE_SWAN/transcriptome_summary")
setwd("3-data_analysis/combined_omics/DE_SWAN")

# library("DEswan")

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
temp_data_transcriptome <-
  temp_data %>%
  dplyr::filter(class == "transcriptome" & p_value_adjust < 0.05)

####crest 1
temp_data_transcriptome_crest1 <-
  temp_data_transcriptome %>%
  dplyr::filter(center == 44)

####crest 2
temp_data_transcriptome_crest2 <-
  temp_data_transcriptome %>%
  dplyr::filter(center == 61)

library(org.Hs.eg.db)
library(clusterProfiler)

load("transcriptome_pathway/GO_result/transcriptome_crest1_go")
load("transcriptome_pathway/KEGG_result/transcriptome_crest1_kegg")
load("transcriptome_pathway/Reactome_result/transcriptome_crest1_reactome")

transcriptome_crest1_go <-
  transcriptome_crest1_go@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(Count >= 3) %>%
  dplyr::arrange(p.adjust)

transcriptome_crest1_kegg <-
  transcriptome_crest1_kegg@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::filter(Count >= 3)

transcriptome_crest1_reactome <-
  transcriptome_crest1_reactome@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  # dplyr::filter(!disease) %>%
  dplyr::filter(Count >= 3) %>%
  dplyr::arrange(p.adjust)

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
load("transcriptome_pathway_crest2/GO_result/transcriptome_crest2_go")
load("transcriptome_pathway_crest2/KEGG_result/transcriptome_crest2_kegg")
load("transcriptome_pathway_crest2/Reactome_result/transcriptome_crest2_reactome")

transcriptome_crest2_go <-
  transcriptome_crest2_go@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(Count >= 3) %>%
  dplyr::arrange(p.adjust)

transcriptome_crest2_kegg <-
  transcriptome_crest2_kegg@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::filter(Count >= 3)

transcriptome_crest2_reactome <-
  transcriptome_crest2_reactome@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  # dplyr::filter(!disease) %>%
  dplyr::filter(Count >= 3) %>%
  dplyr::arrange(p.adjust)

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

dim(result_all_crest1)
dim(result_all_crest2)

result_all_crest1 <-
  result_all_crest1 %>%
  dplyr::arrange(p.adjust)

result_all_crest2 <-
  result_all_crest2 %>%
  dplyr::arrange(p.adjust)

head(result_all_crest1$module_annotation)
head(result_all_crest2$module_annotation)

result_all_crest1$Description[1] %>%
  stringr::str_split(pattern = ";") %>%
  `[[`(1) %>%
  unique()
result_all_crest2$Description[1] %>%
  stringr::str_split(pattern = ";") %>%
  `[[`(1) %>%
  unique()

intersect(result_all_crest1$module_annotation,
          result_all_crest2$module_annotation)

temp_crest1 <-
  result_all_crest1 %>%
  head(20)

temp_crest2 <-
  result_all_crest2 %>%
  head(20)

intersect(temp_crest1$module_annotation,
          temp_crest2$module_annotation)

temp_data <-
  seq_len(nrow(temp_crest1)) %>%
  purrr::map(function(i) {
    x <-
      temp_crest1$geneID[i] %>% stringr::str_split(pattern = "/") %>% `[[`(1)
    
    seq_len(nrow(temp_crest2)) %>%
      purrr::map(function(j) {
        y <-
          temp_crest2$geneID[j] %>% stringr::str_split(pattern = "/") %>% `[[`(1)
        data.frame(
          annotation1 = temp_crest1$module_annotation[i],
          annotation2 = temp_crest2$module_annotation[j],
          jaccard_index = calculate_jaccard_index(x, y)
        )
      }) %>%
      dplyr::bind_rows()
  }) %>%
  dplyr::bind_rows()


edge_data <-
  temp_data %>%
  dplyr::filter(jaccard_index >= 0.5) %>%
  dplyr::mutate(
    from = paste(annotation1, "1", sep = "_"),
    to = paste(annotation2, "2", sep = "_")
  )

node_data1 <-
  temp_crest1 %>%
  dplyr::mutate(module_annotation = paste(module_annotation, "1", sep = "_")) %>%
  dplyr::rename(node = module_annotation) %>%
  dplyr::mutate(database = stringr::str_split(database, ';') %>%
                  lapply(function(x) {
                    x[1]
                  }) %>%
                  unlist()) %>%
  dplyr::select(node, p.adjust, Count, database, module)

node_data2 <-
  temp_crest2 %>%
  dplyr::mutate(module_annotation = paste(module_annotation, "2", sep = "_")) %>%
  dplyr::rename(node = module_annotation) %>%
  dplyr::mutate(database = stringr::str_split(database, ';') %>%
                  lapply(function(x) {
                    x[1]
                  }) %>%
                  unlist()) %>%
  dplyr::select(node, p.adjust, Count, database, module)

node_data <-
  rbind(node_data1,
        node_data2) %>%
  dplyr::mutate(Count = as.numeric(Count))

library(ggraph)
library(tidygraph)
library(igraph)

graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

g <- graph_data

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite")
# dplyr::select(x, y)

coords$x[stringr::str_detect(coords$node, "_1")] <- 0
coords$x[stringr::str_detect(coords$node, "_2")] <- 1

coords$y[coords$x == 0] <-
  1:sum(coords$x == 0)

coords$y[coords$x == 1] <-
  1:sum(coords$x == 1)

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
  )

plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_link(color = "black",
                 aes(width = jaccard_index),
                 show.legend = TRUE) +
  geom_node_point(aes(fill = database,
                      size = Count),
                  shape = 21,
                  show.legend = TRUE) +
  scale_fill_manual(values = database_color) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = stringr::str_replace_all(node, "_[0-2]{1}", "")
    ),
    # size = 3,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(
    edge_width = guide_legend(title = "Jaccard index",
                              override.aes = list(shape = NA)),
    fill = guide_legend(
      title = "Database class",
      override.aes = list(size = 7, linetype = "blank")
    ),
    size = guide_legend(title = "Gene counts",
                        override.aes = list(linetype = 0))
  ) +
  ggraph::scale_edge_color_gradientn(colours = pal) +
  ggraph::scale_edge_width(range = c(0.1, 1.5)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom"
  )

plot

ggsave(plot,
       filename = "transcriptome_summary/top20_pathways.pdf",
       width = 2,
       height = 7)
