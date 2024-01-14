no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
setwd(masstools::get_project_wd())

source("1-code/100-tools.R")

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("3-data_analysis/combined_omics/DE_SWAN")
setwd("3-data_analysis/combined_omics/DE_SWAN")

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
temp_data_lipidomics <-
  temp_data %>%
  dplyr::filter(class == "lipidomics" & p_value_adjust < 0.05)

temp_data_lipidomics$Lipid_Name =
  temp_data_lipidomics$Lipid_Name %>%
  stringr::str_replace_all("\\_", "\\/")

temp_data_lipidomics$Lipid_Name[grep("TAG", temp_data_lipidomics$Lipid_Name)] =
  temp_data_lipidomics$Lipid_Name[grep("TAG", temp_data_lipidomics$Lipid_Name)] %>%
  purrr::map(function(x) {
    # main =
    stringr::str_extract(x, "TAG[0-9]{1,3}\\.[0-9]{1,2}") %>%
      stringr::str_replace("TAG", "") %>%
      stringr::str_split("\\.") %>%
      `[[`(1) %>%
      paste(collapse = ":") %>%
      paste("TAG(", ., ")", sep = "")
    
    # chain =
    #   stringr::str_extract(x, "FA[0-9]{1,3}\\.[0-9]{1,2}") %>%
    #   stringr::str_replace("FA", "") %>%
    #   stringr::str_split("\\.") %>%
    #   `[[`(1) %>%
    #   paste(collapse = ":") %>%
    #   paste("(FA ",.,")", sep = "")
    # paste(main, chain, sep = " ")
  }) %>%
  unlist()

###TAG to TG and DAG to TG
temp_data_lipidomics$Lipid_Name =
  temp_data_lipidomics$Lipid_Name %>%
  stringr::str_replace_all("TAG", "TG") %>%
  stringr::str_replace_all("DAG", "DG")

####lipidomics
dir.create("lipidomics_pathway_crest1")
dir.create("lipidomics_pathway_crest1/lipid_minion_result/")

####crest 1
temp_data_lipidomics_crest1 <-
  temp_data_lipidomics %>%
  dplyr::filter(center == 45)

####lipid_minion
temp_data_lipidomics_crest1$Lipid_Name

# write.table(
#   temp_data_lipidomics_crest1 %>%
#     dplyr::select(Lipid_Name) %>%
#     dplyr::filter(!is.na(Lipid_Name)) %>%
#     dplyr::rename(lipid = Lipid_Name) %>%
#     dplyr::distinct(lipid),
#   file = "lipidomics_pathway_crest1/lipid_minion_result/lipid_list.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )

universe_lipid =
  variable_info %>%
  dplyr::select(Lipid_Name) %>%
  dplyr::filter(!is.na(Lipid_Name)) %>%
  dplyr::rename(lipid = Lipid_Name) %>%
  dplyr::distinct(lipid)

# write.table(
#   temp_data_lipidomics %>%
#     dplyr::select(Lipid_Name) %>%
#     dplyr::filter(!is.na(Lipid_Name)) %>%
#     dplyr::rename(lipid = Lipid_Name) %>%
#     dplyr::distinct(lipid),
#   file = "lipidomics_pathway_crest1/lipid_minion_result/universal_lipid_list.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )

#####read the lipid minion results
enriched_result <-
  readr::read_delim(
    "lipidomics_pathway_crest1/lipid_minion_result/Fisher output table (0 pvals _ 0.05).txt",
    delim = " "
  )

network_node <-
  readr::read_delim("lipidomics_pathway_crest1/lipid_minion_result/Network_nodes.txt",
                    delim = "\t")
network_edge <-
  readr::read_delim("lipidomics_pathway_crest1/lipid_minion_result/Network_edges.txt",
                    delim = "\t")
network_edge_attr <-
  readr::read_delim(
    "lipidomics_pathway_crest1/lipid_minion_result/Network_edge_attributes.txt",
    delim = "\t"
  )

unique(network_node$title)

# enriched_result$Classifier[1] = "Free facty acid"
# enriched_result$Classifier[3] = "FFA"
# enriched_result$Classifier[5] = "FFA("

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
# addWorksheet(wb, sheetName = "Lipid enrichment result", gridLines = TRUE)
# freezePane(wb,
#            sheet = 1,
#            firstRow = TRUE,
#            firstCol = TRUE)
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = enriched_result,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# saveWorkbook(wb, "lipidomics_pathway_crest1/lipid_minion_result/lipid_enrichment.xlsx", overwrite = TRUE)

###network
network_node

library(igraph)
library(ggraph)
library(tidygraph)

network_node

node =
  network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class"))

# node$title[node$title == "Uncategorized"] = "Free facty acid"

edge =
  network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

node <-
  node %>%
  dplyr::filter(!stringr::str_detect(title, "\\($")) %>%
  dplyr::filter(!stringr::str_detect(title, "^[0-9]{1,3}:[0-9]{1,3}$"))

edge <-
  edge %>%
  dplyr::filter(from %in% node$node & to %in% node$node)

netwrok =
  tidygraph::tbl_graph(nodes = node,
                       edges = edge)

plot =
  ggraph(netwrok,
         layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "grey"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class == "class", title, NA)
    ),
    color = "black",
    bg.color = "white",
    size = 5.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("lipid" = unname(omics_color["lipidomics"]),
                                "class" = "red")) +
  scale_size_manual(values = c("lipid" = 2,
                               "class" = 5)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "lipidomics_pathway_crest1/lipid_minion_result/network.pdf",
       width = 12,
       height = 7)
