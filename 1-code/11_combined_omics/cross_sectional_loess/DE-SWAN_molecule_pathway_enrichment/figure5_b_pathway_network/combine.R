no_source()
rm(list = ls())
gc()
library(tidyverse)
library(tidymass)
library(plyr)
setwd(r4projects::get_project_wd())

source("1-code/100-tools.R")

dir.create("3-data_analysis/combined_omics/DE_SWAN")
dir.create(
  "3-data_analysis/combined_omics/DE_SWAN/transcriptome_proteome_metabolome_summary"
)
setwd("3-data_analysis/combined_omics/DE_SWAN/transcriptome_proteome_metabolome_summary/")

load("transcriptomics_node")
load("transcriptomics_edge")

temp = 
transcriptomics_edge %>% 
  dplyr::filter(stringr::str_detect(from, "crest") & stringr::str_detect(to, "crest")) %>% 
  dplyr::mutate(from = stringr::str_extract(from, "crest[0-2]{1,2}"),
                to = stringr::str_extract(to, "crest[0-2]{1,2}"))

dim(temp)
sum(temp$from == temp$to)

load("proteomics_node")
load("proteomics_edge")

temp = 
  proteomics_edge %>% 
  dplyr::filter(stringr::str_detect(from, "crest") & stringr::str_detect(to, "crest")) %>% 
  dplyr::mutate(from = stringr::str_extract(from, "crest[0-2]{1,2}"),
                to = stringr::str_extract(to, "crest[0-2]{1,2}"))

dim(temp)
sum(temp$from == temp$to)

load("metabolomics_node")
load("metabolomics_edge")

temp = 
  metabolomics_edge %>% 
  dplyr::filter(stringr::str_detect(from, "crest") & stringr::str_detect(to, "crest")) %>% 
  dplyr::mutate(from = stringr::str_extract(from, "crest[0-2]{1,2}"),
                to = stringr::str_extract(to, "crest[0-2]{1,2}"))

dim(temp)
sum(temp$from == temp$to)


colnames(transcriptomics_edge)
colnames(metabolomics_edge)

transcriptomics_edge <-
transcriptomics_edge %>% 
  dplyr::left_join(transcriptomics_node[,c("node_id", "class")],
                   by = c("from" = "node_id")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::left_join(transcriptomics_node[,c("node_id", "class")],
                   by = c("to" = "node_id")) %>% 
  dplyr::rename(to_class = class) %>% 
  dplyr::mutate(class = paste(from_class, to_class, sep = "_")) %>% 
  dplyr::select(-c(from_class, to_class))

transcriptomics_edge$class[transcriptomics_edge$class !="function_module"] <-
  "transcriptome"


proteomics_edge <-
  proteomics_edge %>% 
  dplyr::left_join(proteomics_node[,c("node_id", "class")],
                   by = c("from" = "node_id")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::left_join(proteomics_node[,c("node_id", "class")],
                   by = c("to" = "node_id")) %>% 
  dplyr::rename(to_class = class) %>% 
  dplyr::mutate(class = paste(from_class, to_class, sep = "_")) %>% 
  dplyr::select(-c(from_class, to_class))

proteomics_edge$class[proteomics_edge$class !="function_module"] <-
  "transcriptome"

proteomics_edge <-
  proteomics_edge %>% 
  dplyr::left_join(proteomics_node[,c("node_id", "class")],
                   by = c("from" = "node_id")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::left_join(proteomics_node[,c("node_id", "class")],
                   by = c("to" = "node_id")) %>% 
  dplyr::rename(to_class = class) %>% 
  dplyr::mutate(class = paste(from_class, to_class, sep = "_")) %>% 
  dplyr::select(-c(from_class, to_class))

proteomics_edge$class[proteomics_edge$class !="function_module"] <-
  "proteomics"

metabolomics_edge <-
  metabolomics_edge %>% 
  dplyr::left_join(metabolomics_node[,c("node_id", "class")],
                   by = c("from" = "node_id")) %>% 
  dplyr::rename(from_class = class) %>% 
  dplyr::left_join(metabolomics_node[,c("node_id", "class")],
                   by = c("to" = "node_id")) %>% 
  dplyr::rename(to_class = class) %>% 
  dplyr::mutate(class = paste(from_class, to_class, sep = "_")) %>% 
  dplyr::select(-c(from_class, to_class))

metabolomics_edge$class[metabolomics_edge$class !="function_module"] <-
  "metabolomics"


edge_data <-
  rbind(
    transcriptomics_edge,
    proteomics_edge,
    metabolomics_edge
  )

colnames(transcriptomics_node)
colnames(proteomics_node)
colnames(metabolomics_node)

node_data <-
rbind(
  transcriptomics_node %>% 
    dplyr::select(colnames(metabolomics_node)) %>% 
    dplyr::mutate(omics_class = "transcriptome"),
  proteomics_node %>% 
    dplyr::select(colnames(metabolomics_node)) %>% 
    dplyr::mutate(omics_class = "proteomics"),
  metabolomics_node %>% 
    dplyr::mutate(omics_class = "metabolomics")
) %>% 
  dplyr::distinct(node_id, .keep_all = TRUE)

node_data$omics_class[node_data$class == "function"] <-"function"

node_data$label <-
  node_data$label %>% 
  stringr::str_replace("_crest[1-2]{1,2}", "")

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

coords$y[coords$crest == 1 & coords$class == "function"] <- 6
coords$y[coords$crest == 1 & coords$class == "module"] <- 5
coords$y[coords$crest == 1 & coords$class == "pathway"] <- 4
coords$y[coords$class == "gene"] <- 3
coords$y[coords$class == "protein"] <- 3
coords$y[coords$class == "metabolite"] <- 3
coords$y[coords$crest == 2 & coords$class == "pathway"] <- 2
coords$y[coords$crest == 2 & coords$class == "module"] <- 1
coords$y[coords$crest == 2 & coords$class == "function"] <- 0

range(coords$x[coords$y == 3])

##function
coords$x[coords$y == 0] <-
  new_coords(range_left = 500, range_right = 1500, old_x = coords$x[coords$y == 0])
coords$x[coords$y == 6] <-
  new_coords(range_left = 500, range_right = 1500, old_x = coords$x[coords$y == 6])

##module
coords$x[coords$y == 1] <-
  new_coords(range_left = 200, range_right = 1800, old_x = coords$x[coords$y == 1])
coords$x[coords$y == 5] <-
  new_coords(range_left = 200, range_right = 1800, old_x = coords$x[coords$y == 5])

##pathway
##crest 1
##transcriptomics
range(coords$x[coords$y == 2 & coords$omics_class == "transcriptome"])
coords$x[coords$y == 2 & coords$omics_class == "transcriptome"] <-
  new_coords(range_left = 50,
             range_right = 1000,
             old_x = coords$x[coords$y == 2 & coords$omics_class == "transcriptome"],
             scale = TRUE)
###proteomics
range(coords$x[coords$y == 2 & coords$omics_class == "proteomics"])
coords$x[coords$y == 2 & coords$omics_class == "proteomics"] <-
  new_coords(range_left = 1010,
             range_right = 1700,
             old_x = coords$x[coords$y == 2 & coords$omics_class == "proteomics"],
             scale = TRUE)
###metabolomics
range(coords$x[coords$y == 2 & coords$omics_class == "metabolomics"])
coords$x[coords$y == 2 & coords$omics_class == "metabolomics"] <-
  new_coords(range_left = 1710,
             range_right = 1950,
             old_x = coords$x[coords$y == 2 & coords$omics_class == "metabolomics"],
             scale = FALSE)

# coords$x[coords$y == 2] <-
#   new_coords(range_left = 50,
#              range_right = 1950,
#              old_x = coords$x[coords$y == 2])

##crest 2
###transcriptomics
range(coords$x[coords$y == 4 & coords$omics_class == "transcriptome"])
coords$x[coords$y == 4 & coords$omics_class == "transcriptome"] <-
  new_coords(range_left = 50,
             range_right = 900,
             old_x = coords$x[coords$y == 4 & coords$omics_class == "transcriptome"],
             scale = TRUE)
###proteomics
range(coords$x[coords$y == 4 & coords$omics_class == "proteomics"])
coords$x[coords$y == 4 & coords$omics_class == "proteomics"] <-
  new_coords(range_left = 910,
             range_right = 1800,
             old_x = coords$x[coords$y == 4 & coords$omics_class == "proteomics"],
             scale = TRUE)
###metabolomics
range(coords$x[coords$y == 4 & coords$omics_class == "metabolomics"])
coords$x[coords$y == 4 & coords$omics_class == "metabolomics"] <-
  new_coords(range_left = 1810,
             range_right = 1950,
             old_x = coords$x[coords$y == 4 & coords$omics_class == "metabolomics"],
             scale = FALSE)


# coords$x[coords$y == 4] <-
#   new_coords(range_left = 50,
#              range_right = 1950,
#              old_x = coords$x[coords$y == 4])

# 
###gene
coords$x[coords$y == 3] <-
  new_coords(range_left = 0, 
             range_right = 2000, 
             old_x = coords$x[coords$y == 3], 
             scale = TRUE)

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
  geom_edge_diagonal(show.legend = FALSE,
                 aes(color = class)) +
  geom_node_point(aes(fill = omics_class,
                      size = Count,
                      shape = class),
                  show.legend = TRUE) +
  scale_edge_color_manual(
    values = c("function_module" = unname(network_class_color["function"]),
    omics_color["transcriptome"],
    omics_color["proteomics"],
    omics_color["metabolomics"])) +
  scale_fill_manual(values = c(network_class_color["function"],
                               omics_color["transcriptome"],
                               omics_color["proteomics"],
                               omics_color["metabolomics"])) +
  scale_shape_manual(values = c("function" = 23,
                                "module" = 22,
                                "pathway" = 21,
                                "gene" = 24,
                                "metabolite" = 24,
                                "protein" = 24
                                )) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = label
    ),
    # size = 3,
    alpha = 1,
    show.legend = FALSE,
    angle = 45, hjust = 0, vjust = 0.5
  ) +
  scale_size_continuous(range = c(1, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom"
  )

plot

# ggsave(plot, filename = "network.pdf", width = 10, height = 6)

