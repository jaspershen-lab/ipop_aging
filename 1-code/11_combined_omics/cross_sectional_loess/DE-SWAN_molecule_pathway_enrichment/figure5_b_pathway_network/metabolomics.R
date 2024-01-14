# no_source()
rm(list = ls())
gc()
library(tidyverse)
library(tidymass)
library(plyr)
setwd(masstools::get_project_wd())

source("1-code/100-tools.R")

dir.create("3-data_analysis/combined_omics/DE_SWAN")
dir.create(
  "3-data_analysis/combined_omics/DE_SWAN/transcriptome_proteome_metabolome_summary"
)
setwd("3-data_analysis/combined_omics/DE_SWAN")

####Metabolomics
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

###crest 1
metabolomics_node_crest1 <-
  temp_crest1_metabolome %>%
  dplyr::rename(
    node_id = pathway_id,
    p.adjust = p_value_adjust,
    module_annotation = pathway_name,
    Count = mapped_number
  ) %>%
  dplyr::mutate(crest = 1,
                label = module_annotation) %>%
  dplyr::select(node_id,
                p.adjust,
                module_annotation,
                Count,
                crest,
                label,
                mapped_id) %>%
  dplyr::mutate(node_id =
                  paste("metabolomics_crest1", node_id, sep = "_"))

metabolomics_edge_crest1 <-
  purrr::map2(metabolomics_node_crest1$node_id,
              metabolomics_node_crest1$mapped_id,
              function(x, y) {
                y <- stringr::str_split(y, ";")[[1]]
                data.frame(from = x,
                           to = y)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

metabolomics_node_crest1 <-
  rbind(
    metabolomics_node_crest1 %>% dplyr::select(-mapped_id) %>%
      dplyr::mutate(class = "pathway"),
    data.frame(
      node_id = stringr::str_split(metabolomics_node_crest1$mapped_id, ";") %>%
        unlist(),
      p.adjust = NA,
      module_annotation = NA,
      Count = NA,
      crest = 1,
      label = NA
    ) %>%
      dplyr::mutate(class = "metabolite")
  )

metabolomics_edge_crest1_temp <-
  rbind(
    data.frame(from = "lipid_metabolism_crest1",
               to = c("ABC transporters")),
    data.frame(from = "Caffeine_metabolism_crest1",
               to = c("Caffeine metabolism")),
    data.frame(from = "CVD_crest1",
               to = c("Phenylalanine metabolism"))
  ) %>%
  dplyr::left_join(metabolomics_node_crest1[, c("node_id", "module_annotation")],
                   by = c("to" = "module_annotation")) %>%
  dplyr::select(-to) %>%
  dplyr::rename(to = node_id)

match("Phenylalanine metabolism",
      metabolomics_node_crest1$module_annotation)

metabolomics_node_crest1_temp <-
  data.frame(
    node_id = unique(metabolomics_edge_crest1_temp$from),
    p.adjust = NA,
    module_annotation = unique(metabolomics_edge_crest1_temp$from),
    Count = NA,
    crest = 1,
    label = unique(metabolomics_edge_crest1_temp$from),
    class = "function"
  )

metabolomics_edge_crest1 <-
  rbind(metabolomics_edge_crest1,
        metabolomics_edge_crest1_temp)

metabolomics_node_crest1 <-
  metabolomics_node_crest1 %>%
  dplyr::filter((class != "pathway") |
                  (class == "pathway" &
                     node_id %in% metabolomics_edge_crest1_temp$to)
  )

metabolomics_node_crest1 <-
  rbind(metabolomics_node_crest1,
        metabolomics_node_crest1_temp)

metabolomics_edge_crest1 <-
  metabolomics_edge_crest1 %>% 
  dplyr::filter(from %in% metabolomics_node_crest1$node_id &
                  to %in% metabolomics_node_crest1$node_id)

save(metabolomics_edge_crest1, file = "transcriptome_proteome_metabolome_summary/metabolomics_edge_crest1")
save(metabolomics_node_crest1, file = "transcriptome_proteome_metabolome_summary/metabolomics_node_crest1")

###crest 2
metabolomics_node_crest2 <-
  temp_crest2_metabolome %>%
  dplyr::rename(
    node_id = pathway_id,
    p.adjust = p_value_adjust,
    module_annotation = pathway_name,
    Count = mapped_number
  ) %>%
  dplyr::mutate(crest = 2,
                label = module_annotation) %>%
  dplyr::select(node_id,
                p.adjust,
                module_annotation,
                Count,
                crest,
                label,
                mapped_id) %>%
  dplyr::mutate(node_id =
                  paste("metabolomics_crest2", node_id, sep = "_"))

metabolomics_edge_crest2 <-
  purrr::map2(metabolomics_node_crest2$node_id,
              metabolomics_node_crest2$mapped_id,
              function(x, y) {
                y <- stringr::str_split(y, ";")[[1]]
                data.frame(from = x,
                           to = y)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

metabolomics_node_crest2 <-
  rbind(
    metabolomics_node_crest2 %>% dplyr::select(-mapped_id) %>%
      dplyr::mutate(class = "pathway"),
    data.frame(
      node_id = stringr::str_split(metabolomics_node_crest2$mapped_id, ";") %>%
        unlist(),
      p.adjust = NA,
      module_annotation = NA,
      Count = NA,
      crest = 2,
      label = NA
    ) %>%
      dplyr::mutate(class = "metabolite")
  )

metabolomics_edge_crest2_temp <-
  rbind(
    data.frame(from = "Caffeine_metabolism_crest2",
               to = c("Caffeine metabolism")),
    data.frame(
      from = "CVD_crest2",
      to = c(
        "Phenylalanine metabolism",
        "Valine, leucine and isoleucine biosynthesis",
        "Alpha Linolenic Acid and Linoleic Acid Metabolism",
        "Alanine, aspartate and glutamate metabolism"
      )
    ),
    data.frame(
      from = "Skin-Muscle_crest2",
      to = c("Glycine, serine and threonine metabolism")
    )
  ) %>%
  dplyr::left_join(metabolomics_node_crest2[, c("node_id", "module_annotation")],
                   by = c("to" = "module_annotation")) %>%
  dplyr::select(-to) %>%
  dplyr::rename(to = node_id)

match(
  "Glycine, serine and threonine metabolism",
  metabolomics_node_crest2$module_annotation
)

metabolomics_node_crest2_temp <-
  data.frame(
    node_id = unique(metabolomics_edge_crest2_temp$from),
    p.adjust = NA,
    module_annotation = unique(metabolomics_edge_crest2_temp$from),
    Count = NA,
    crest = 2,
    label = unique(metabolomics_edge_crest2_temp$from),
    class = "function"
  )

metabolomics_node_crest2 <-
  metabolomics_node_crest2 %>%
  dplyr::filter((class != "pathway") |
                  (class == "pathway" &
                     node_id %in% metabolomics_edge_crest2_temp$to)
  )

metabolomics_node_crest2 <-
  rbind(metabolomics_node_crest2,
        metabolomics_node_crest2_temp)

metabolomics_edge_crest2 <-
  rbind(metabolomics_edge_crest2,
        metabolomics_edge_crest2_temp)

metabolomics_edge_crest2 <-
  metabolomics_edge_crest2 %>% 
  dplyr::filter(from %in% metabolomics_node_crest2$node_id &
                  to %in% metabolomics_node_crest2$node_id)

save(metabolomics_edge_crest2, file = "transcriptome_proteome_metabolome_summary/metabolomics_edge_crest2")
save(metabolomics_node_crest2, file = "transcriptome_proteome_metabolome_summary/metabolomics_node_crest2")

###combine
#####combine crest1 and crest2
metabolomics_node <-
  rbind(metabolomics_node_crest1,
        metabolomics_node_crest2) %>%
  dplyr::distinct(node_id, .keep_all = TRUE)

metabolomics_node$label[metabolomics_node$class == "function"] <-
  metabolomics_node$node_id[metabolomics_node$class == "function"] %>%
  stringr::str_replace("_crest[0-1]{1,2}", "")

# metabolomics_node$label[metabolomics_node$class == "pathway"] <-
#   NA

metabolomics_node$Count[metabolomics_node$class == "metabolite"] <-
  1

metabolomics_node$Count[metabolomics_node$class == "function"] <-
  300

metabolomics_edge <-
  rbind(metabolomics_edge_crest1,
        metabolomics_edge_crest2)

metabolomics_node <-
metabolomics_node %>% 
  dplyr::filter(node_id %in% c(metabolomics_edge$from, metabolomics_edge$to))

save(metabolomics_node, file = "transcriptome_proteome_metabolome_summary/metabolomics_node")
save(metabolomics_edge, file = "transcriptome_proteome_metabolome_summary/metabolomics_edge")

library(tidygraph)
library(ggraph)
library(igraph)

graph_data <-
  tidygraph::tbl_graph(nodes = metabolomics_node,
                       edges = metabolomics_edge,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

g <- graph_data

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite")
# dplyr::select(x, y)

coords$y[coords$crest == 1 & coords$class == "function"] <- 0
coords$y[coords$crest == 1 & coords$class == "module"] <- 1
coords$y[coords$crest == 1 & coords$class == "pathway"] <- 2
coords$y[coords$class == "metabolite"] <- 3
coords$y[coords$crest == 2 & coords$class == "pathway"] <- 4
coords$y[coords$crest == 2 & coords$class == "module"] <- 5
coords$y[coords$crest == 2 & coords$class == "function"] <- 6

range(coords$x[coords$y == 3])

##function
coords$x[coords$y == 0] <-
  new_coords(range_left = 250,
             range_right = 750,
             old_x = coords$x[coords$y == 0])
coords$x[coords$y == 6] <-
  new_coords(range_left = 250,
             range_right = 750,
             old_x = coords$x[coords$y == 6])

##module
coords$x[coords$y == 1] <-
  new_coords(range_left = 100,
             range_right = 900,
             old_x = coords$x[coords$y == 1])
coords$x[coords$y == 5] <-
  new_coords(range_left = 100,
             range_right = 900,
             old_x = coords$x[coords$y == 5])

##pathway
coords$x[coords$y == 2] <-
  new_coords(range_left = 50,
             range_right = 950,
             old_x = coords$x[coords$y == 2])
coords$x[coords$y == 4] <-
  new_coords(range_left = 50,
             range_right = 950,
             old_x = coords$x[coords$y == 4])

###metabolite
coords$x[coords$y == 3] <-
  new_coords(range_left = 0,
             range_right = 1000,
             old_x = coords$x[coords$y == 3])

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
  )

# plot <-
ggraph(my_graph,
       layout = 'bipartite') +
  geom_edge_diagonal(color = "grey",
                     show.legend = TRUE) +
  geom_node_point(shape = 21,
                  aes(fill = class,
                      size = Count),
                  show.legend = TRUE) +
  scale_fill_manual(values = network_class_color) +
  geom_node_text(
    aes(x = x,
        y = y,
        label = label),
    # size = 3,
    alpha = 1,
    show.legend = FALSE,
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ) +
  scale_size_continuous(range = c(1, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom"
  )

