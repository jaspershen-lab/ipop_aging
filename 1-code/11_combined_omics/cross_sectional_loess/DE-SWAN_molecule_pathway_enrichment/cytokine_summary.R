no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
library(plyr)
setwd(r4projects::get_project_wd())

source("1-code/100-tools.R")

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

cytokine_class <-
  readxl::read_xlsx("3-data_analysis/plasma_cytokine/data_preparation/cy_cla.xlsx")

dir.create("3-data_analysis/combined_omics/DE_SWAN")
dir.create("3-data_analysis/combined_omics/DE_SWAN/cytokine_summary")
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

variable_info <-
  variable_info %>%
  dplyr::filter(class == "cytokine")

expression_data <-
  expression_data[variable_info$variable_id,]

match(variable_info$variable_id,
      paste("cytokine", cytokine_class$variable_id, sep = "_"))

variable_info <-
  variable_info %>%
  dplyr::left_join(cytokine_class %>%
                     dplyr::mutate(variable_id =
                                     paste("cytokine", variable_id, sep = "_")),
                   by = "variable_id")

load("temp_data")

temp_data <-
  temp_data %>%
  dplyr::left_join(variable_info,
                   by = c("variable_id"))

###Pathway enrichment
temp_data_cytokine <-
  temp_data %>%
  dplyr::filter(class == "cytokine" & p_value_adjust < 0.05)

####crest 1
temp_data_cytokine_crest1 <-
  temp_data_cytokine %>%
  dplyr::filter(center == 51)

####crest 2
temp_data_cytokine_crest2 <-
  temp_data_cytokine %>%
  dplyr::filter(center == 61)

cytokine_crest1 <-
  readxl::read_xlsx("cytokine_pathway_crest1/cytokine_crest1.xlsx") %>%
  dplyr::select(variable_id) %>%
  dplyr::left_join(variable_info[, c("variable_id", "classification")],
                   by = "variable_id") %>%
  # dplyr::mutate(variable_id = paste0(variable_id, "_1")) %>%
  dplyr::mutate(crest = "1") %>%
  dplyr::filter(!stringr::str_detect(variable_id, "CHEX"))

cytokine_crest2 <-
  readxl::read_xlsx("cytokine_pathway_crest2/cytokine_crest2.xlsx") %>%
  dplyr::select(variable_id) %>%
  dplyr::left_join(variable_info[, c("variable_id", "classification")],
                   by = "variable_id") %>%
  # dplyr::mutate(variable_id = paste0(variable_id, "_2")) %>%
  dplyr::mutate(crest = "2") %>%
  dplyr::filter(!stringr::str_detect(variable_id, "CHEX"))

temp_data <-
  cytokine_crest1 %>%
  dplyr::full_join(cytokine_crest2, by = c("variable_id", "classification")) %>% 
  dplyr::rename(crest1 = crest.x,
                crest2 = crest.y)

temp_data$crest1[which(!is.na(temp_data$crest1))] <-
  temp_data$classification[which(!is.na(temp_data$crest1))]

temp_data$crest2[which(!is.na(temp_data$crest2))] <-
  temp_data$classification[which(!is.na(temp_data$crest2))]

temp_data <-
  temp_data %>% 
  make_long(crest1, crest2)

plot <- 
ggplot(
  temp_data,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = node
  )
) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = c(cytokine_class_color)) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) 

ggsave(plot, filename = "cytokine_summary/cytokine_sankey.pdf", width = 7, height = 3)

dim(cytokine_crest1)
dim(cytokine_crest2)

setdiff(cytokine_crest1$variable_id,
          cytokine_crest2$variable_id)

setdiff(cytokine_crest2$variable_id,
        cytokine_crest1$variable_id)

intersect(cytokine_crest2$variable_id,
        cytokine_crest1$variable_id) %>% 
  length()

table(cytokine_crest1$classification)

intersect(cytokine_crest1$variable_id[cytokine_crest1$classification == "Proinflammatory"],
          cytokine_crest2$variable_id[cytokine_crest2$classification == "Proinflammatory"]) %>% 
  length()

intersect(cytokine_crest1$variable_id[cytokine_crest1$classification == "Proinflammatory/Anti-inflammatory"],
          cytokine_crest2$variable_id[cytokine_crest2$classification == "Proinflammatory/Anti-inflammatory"]) %>% 
  length()

