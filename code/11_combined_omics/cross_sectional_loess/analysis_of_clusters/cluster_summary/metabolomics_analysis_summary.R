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

dir.create(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_summary/metabolomics",
  recursive = TRUE
)

setwd(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_summary/metabolomics"
)

object_cross_section_loess

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

load("../../final_cluster_info")

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

table(final_cluster_info$cluster, final_cluster_info$class)

###cluster 2
load("../../cluster_2/metabolomics_pathway/cluster2_metabolomics")
load("../../cluster_2/metabolomics_pathway/cluster2_pathway")

load("../../cluster_4/metabolomics_pathway/cluster4_metabolomics")
load("../../cluster_4/metabolomics_pathway/cluster4_pathway")

load("../../cluster_5/metabolomics_pathway/cluster5_metabolomics")
load("../../cluster_5/metabolomics_pathway/cluster5_pathway")

temp_data <-
  rbind(
    data.frame(cluster2_pathway, cluster = "2"),
    data.frame(cluster4_pathway, cluster = "4"),
    data.frame(cluster5_pathway, cluster = "5")
  )

library(ggh4x)

plot <-
  temp_data %>%
  dplyr::mutate(p_value_adjust = -log(p_value_adjust, 10)) %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::mutate(pathway_name = factor(pathway_name, levels = pathway_name)) %>%
  ggplot(aes(p_value_adjust, pathway_name)) +
  geom_segment(aes(
    x = 0,
    xend = p_value_adjust,
    y = pathway_name,
    yend = pathway_name
  )) +
  geom_point(aes(size = mapped_number,
                 color = class)) +
  scale_size_continuous(range = c(2, 5)) +
  facet_wrap(
    facets = vars(cluster),
    strip.position = "left",
    scales = "free_y",
    ncol = 1,
    shrink = TRUE,
    as.table = TRUE
  ) +
  force_panelsizes(rows = c(1, 4, 3)) +
  scale_color_manual(values = database_color, limits = unique(temp_data$class)) +
  theme_base +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.ticks.y = element_blank()
  ) +
  labs(y = "")

plot


###Caffeine metabolism
caffeine_metabolism_kegg_id <-
  temp_data$mapped_id[3] %>% stringr::str_split(";") %>% `[[`(1)

caffeine_metabolism_metabolites <-
  variable_info %>%
  dplyr::filter(KEGG.ID %in% caffeine_metabolism_kegg_id)


