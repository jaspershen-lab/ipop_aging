no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
library(plyr)
setwd(masstools::get_project_wd())

source("code/tools.R")

load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("data_analysis/combined_omics/DE_SWAN")
dir.create("data_analysis/combined_omics/DE_SWAN/gut_microbiome_summary")
setwd("data_analysis/combined_omics/DE_SWAN")

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
  dplyr::filter(class == "gut_microbiome")

expression_data <-
  expression_data[variable_info$variable_id,]

load("temp_data")

temp_data <-
  temp_data %>%
  dplyr::left_join(variable_info,
                   by = c("variable_id"))

###Pathway enrichment
temp_data_gut_microbiome <-
  temp_data %>%
  dplyr::filter(class == "gut_microbiome" & p_value_adjust < 0.05)

####crest 1
temp_data_gut_microbiome_crest1 <-
  temp_data_gut_microbiome %>%
  dplyr::filter(center == 43) %>%
  dplyr::select(variable_id, Phylum)
# dplyr::mutate(crest = "crest1",
#               variable_id_new = paste0(variable_id, "_1"))

####crest 2
temp_data_gut_microbiome_crest2 <-
  temp_data_gut_microbiome %>%
  dplyr::filter(center == 61) %>%
  dplyr::select(variable_id, Phylum)
# dplyr::mutate(crest = "crest2",
#               variable_id_new = paste0(variable_id, "_2"))

temp_data <-
  temp_data_gut_microbiome_crest1 %>%
  dplyr::full_join(temp_data_gut_microbiome_crest2, by = "variable_id") %>%
  dplyr::rename(crest1 = Phylum.x,
                crest2 = Phylum.y)

temp_data$crest1[is.na(temp_data$crest1)] <- "No"
temp_data$crest2[is.na(temp_data$crest2)] <- "No"

temp_data <-
  temp_data %>%
  dplyr::group_by(crest1, crest2) %>%
  dplyr::summarise(Freq = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(class = case_when(crest1 != "No" ~ crest1,
                                  crest2 != "No" ~ crest2))

library(ggalluvial)

is_alluvia_form(as.data.frame(temp_data),
                axes = 1:3,
                silent = TRUE)

unique(temp_data$class)

# temp_data <-
# temp_data %>%
#   dplyr::filter(crest1 != "No" & crest2 != "No")

plot <-
  ggplot(as.data.frame(temp_data),
         aes(y = Freq, axis1 = crest1, axis2 = crest2)) +
  geom_alluvium(aes(fill = class), width = 3 / 12) +
  geom_stratum(width = 3 / 12,
               aes(fill = class),
               color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("crest1", "crest2"),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = genus_color) +
  theme_base +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "bottom", panel.grid = element_blank())

plot

ggsave(plot, file = "gut_microbiome_summary/sankey_plot.pdf", width = 7, height = 5)
