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

dir.create(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/metabolomics_summary"
)
setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/metabolomics_summary"
)

object_cross_section_loess

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

###cluster 2, 4, 5

load("../cluster_2/metabolomics_pathway/cluster2_pathway")

load("../cluster_4/metabolomics_pathway/cluster4_pathway")

load("../cluster_5/metabolomics_pathway/cluster5_pathway")

<<<<<<< HEAD
temp_data <-
  rbind(data.frame(cluster2_pathway, cluster = 2),
        data.frame(cluster4_pathway, cluster = 4))

edge_data <-
  seq_len(nrow(temp_data)) %>%
  purrr::map(function(i) {
    data.frame(
      from = temp_data$pathway_id[i],
      to = stringr::str_split(temp_data$mapped_id[i], ";")[[1]]
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

node_data1 <-
  data.frame(
    node = temp_data$pathway_id,
    node_name = temp_data$pathway_name,
    p_value_adjust = -log(temp_data$p_value_adjust, 10),
    mapped_percentage = temp_data$mapped_percentage,
    class = "pathway"
  )

node_data2 <-
  unique(unlist(stringr::str_split(temp_data$mapped_id, ";"))) %>%
  data.frame(node = .) %>%
  left_join(object_cross_section_loess@variable_info[, c("KEGG.ID", "Compound.name")],
            by = c("node" = "KEGG.ID")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  left_join(object_cross_section_loess@variable_info[, c("HMDB.ID", "Compound.name")],
            by = c("node" = "HMDB.ID")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::mutate(Compound.name = case_when(
    !is.na(Compound.name.x) ~ Compound.name.x,
    !is.na(Compound.name.y) ~ Compound.name.y
  )) %>%
  dplyr::select(-c(Compound.name.x, Compound.name.y)) %>%
  dplyr::rename(node_name = Compound.name) %>%
  dplyr::mutate(
    p_value_adjust = 1,
    mapped_percentage = NA,
    class = "metabolite"
  )

node_data <-
  rbind(node_data1,
        node_data2)

library(tidygraph)
library(ggraph)
library(igraph)

graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE)

plot <-
  ggraph(graph_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "grey"
  ) +
  geom_node_point(
    aes(
      color = class,
      size = p_value_adjust,
      shape = class
    ),
    shape = 16,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      angle = node_angle(x, y),
      label = node_name
    ),
    hjust = 'outward',
    size = 3
  ) +
  ggraph::theme_graph()

plot

library(extrafont)
extrafont::loadfonts()
ggsave(
  plot = plot,
  filename = "pathway.pdf",
  width = 6.6,
  height = 5.38
)











###cluster 1, 3, 6, 7, 8, 9, 10, 11

load("../cluster_1/metabolomics_pathway/cluster1_pathway")
load("../cluster_3/metabolomics_pathway/cluster3_pathway")
load("../cluster_6/metabolomics_pathway/cluster6_pathway")
load("../cluster_7/metabolomics_pathway/cluster7_pathway")
load("../cluster_8/metabolomics_pathway/cluster8_pathway")
load("../cluster_9/metabolomics_pathway/cluster9_pathway")
load("../cluster_10/metabolomics_pathway/cluster10_pathway")
load("../cluster_11/metabolomics_pathway/cluster11_pathway")

=======
>>>>>>> abb848b5143c5dd49832cc7a6dee08d71aa25e89
temp_data <-
  rbind(data.frame(cluster2_pathway, cluster = 2),
        data.frame(cluster4_pathway, cluster = 4))

edge_data <-
  seq_len(nrow(temp_data)) %>%
  purrr::map(function(i) {
    data.frame(
      from = temp_data$pathway_id[i],
      to = stringr::str_split(temp_data$mapped_id[i], ";")[[1]]
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

node_data1 <-
  data.frame(
    node = temp_data$pathway_id,
    node_name = temp_data$pathway_name,
    p_value_adjust = -log(temp_data$p_value_adjust, 10),
    mapped_percentage = temp_data$mapped_percentage,
    class = "pathway"
  )

node_data2 <-
  unique(unlist(stringr::str_split(temp_data$mapped_id, ";"))) %>%
  data.frame(node = .) %>%
  left_join(object_cross_section_loess@variable_info[, c("KEGG.ID", "Compound.name")],
            by = c("node" = "KEGG.ID")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  left_join(object_cross_section_loess@variable_info[, c("HMDB.ID", "Compound.name")],
            by = c("node" = "HMDB.ID")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::mutate(Compound.name = case_when(
    !is.na(Compound.name.x) ~ Compound.name.x,
    !is.na(Compound.name.y) ~ Compound.name.y
  )) %>%
  dplyr::select(-c(Compound.name.x, Compound.name.y)) %>%
  dplyr::rename(node_name = Compound.name) %>%
  dplyr::mutate(
    p_value_adjust = 1,
    mapped_percentage = NA,
    class = "metabolite"
  )

node_data <-
  rbind(node_data1,
        node_data2)

library(tidygraph)
library(ggraph)
library(igraph)

graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE)

plot <-
  ggraph(graph_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "grey"
  ) +
  geom_node_point(
    aes(
      color = class,
      size = p_value_adjust,
      shape = class
    ),
    shape = 16,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      angle = node_angle(x, y),
      label = node_name
    ),
    hjust = 'outward',
    size = 3
  ) +
  ggraph::theme_graph()

plot

library(extrafont)
extrafont::loadfonts()
ggsave(
  plot = plot,
  filename = "pathway.pdf",
  width = 6.6,
  height = 5.38
)











###cluster 1, 3, 6, 7, 8, 9, 10, 11
load("../cluster_1/metabolomics_pathway/cluster1_pathway")
load("../cluster_3/metabolomics_pathway/cluster3_pathway")
load("../cluster_6/metabolomics_pathway/cluster6_pathway")
load("../cluster_7/metabolomics_pathway/cluster7_pathway")
load("../cluster_8/metabolomics_pathway/cluster8_pathway")
load("../cluster_9/metabolomics_pathway/cluster9_pathway")
load("../cluster_10/metabolomics_pathway/cluster10_pathway")
load("../cluster_11/metabolomics_pathway/cluster11_pathway")

temp_data <-
  rbind(data.frame(cluster3_pathway, cluster = 3),
        data.frame(cluster7_pathway, cluster = 7),
        data.frame(cluster8_pathway, cluster = 8),
        data.frame(cluster9_pathway, cluster = 9),
        data.frame(cluster10_pathway, cluster = 10))

temp_data$pathway_id <-
  paste(temp_data$pathway_id, temp_data$cluster, sep = "_")

temp_data$pathway_name <-
  paste(temp_data$pathway_name, temp_data$cluster, sep = "_")

edge_data <-
  seq_len(nrow(temp_data)) %>%
  purrr::map(function(i) {
    data.frame(
      from = temp_data$pathway_id[i],
      to = stringr::str_split(temp_data$mapped_id[i], ";")[[1]]
    )
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

node_data1 <-
  data.frame(
    node = temp_data$pathway_id,
    node_name = temp_data$pathway_name,
    p_value_adjust = -log(temp_data$p_value_adjust, 10),
    mapped_percentage = temp_data$mapped_percentage,
    class = "pathway"
  )

node_data2 <-
  unique(unlist(stringr::str_split(temp_data$mapped_id, ";"))) %>%
  data.frame(node = .) %>%
  left_join(object_cross_section_loess@variable_info[, c("KEGG.ID", "Compound.name")],
            by = c("node" = "KEGG.ID")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  left_join(object_cross_section_loess@variable_info[, c("HMDB.ID", "Compound.name")],
            by = c("node" = "HMDB.ID")) %>%
  dplyr::distinct(node, .keep_all = TRUE) %>%
  dplyr::mutate(Compound.name = case_when(
    !is.na(Compound.name.x) ~ Compound.name.x,
    !is.na(Compound.name.y) ~ Compound.name.y
  )) %>%
  dplyr::select(-c(Compound.name.x, Compound.name.y)) %>%
  dplyr::rename(node_name = Compound.name) %>%
  dplyr::mutate(
    p_value_adjust = 1,
    mapped_percentage = NA,
    class = "metabolite"
  )

node_data <-
  rbind(node_data1,
        node_data2)

library(tidygraph)
library(ggraph)
library(igraph)

graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE)

plot <-
  ggraph(graph_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "grey"
  ) +
  geom_node_point(
    aes(
      color = class,
      size = p_value_adjust,
      shape = class
    ),
    shape = 16,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      angle = node_angle(x, y),
      label = node_name
    ),
    hjust = 'outward',
    size = 3
  ) +
  ggraph::theme_graph()

plot

library(extrafont)
extrafont::loadfonts()
ggsave(
  plot = plot,
  filename = "pathway2.pdf",
  width = 6.6,
  height = 5.38
)
