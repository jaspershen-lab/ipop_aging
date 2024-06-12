no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/class_distributation",
  recursive = TRUE
)

setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/class_distributation"
)

load("../final_cluster_info")

variable_info <-
  object_cross_section_loess@variable_info

temp_data <-
  final_cluster_info[, c("variable_id", "cluster")] %>% 
  dplyr::left_join(variable_info[, c("variable_id", "class")])

library(plyr)

plot <-
  temp_data %>%
  dplyr::select(-variable_id) %>% 
  dplyr::mutate(cluster = as.character(cluster)) %>%
  plyr::dlply(.variables = .(cluster)) %>% 
  purrr::map(function(x){
    x %>% 
      dplyr::group_by(class) %>% 
      dplyr::summarise(n = n()) %>% 
      dplyr::mutate(freq = n*100/sum(n)) %>% 
      dplyr::mutate(cluster = x$cluster[1])
  }) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(class = factor(class, names(omics_color))) %>% 
  ggplot(aes(y = cluster, x = freq)) +
  geom_bar(aes(fill = class), stat = "identity") +
  facet_grid(rows = vars(cluster),
             scales = "free") +
  scale_fill_manual(values = omics_color) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme_base +
  theme(panel.grid = element_blank())

plot

# ggsave(plot,
#        filename = "class_distributation.pdf",
#        width = 7,
#        height = 5)



