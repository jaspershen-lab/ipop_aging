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
dir.create("data_analysis/combined_omics/DE_SWAN/metabolomics_summary")
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

load("temp_data")


##remove the background
background <-
  temp_data %>%
  group_by(variable_id) %>%
  dplyr::summarise(n = sum(p_value_adjust < 0.05))

sum(background$n > length(unique(temp_data$center)) * 0.8)

background <-
  background %>%
  dplyr::filter(n > length(unique(temp_data$center)) * 0.8)

temp_data <-
  temp_data %>%
  dplyr::filter(!variable_id %in% background$variable_id)

temp_data <-
  temp_data %>%
  dplyr::left_join(variable_info,
                   by = c("variable_id"))

###Pathway enrichment
temp_data_metabolomics <-
  temp_data %>%
  dplyr::filter(class == "metabolomics" & p_value_adjust < 0.05)

####crest 1
temp_data_metabolomics_crest1 <-
  temp_data_metabolomics %>%
  dplyr::filter(center == 47)

####crest 2
temp_data_metabolomics_crest2 <-
  temp_data_metabolomics %>%
  dplyr::filter(center == 57)

load("metabolomics_pathway_crest1/KEGG_result/metabolomics_crest1_kegg")
load("metabolomics_pathway_crest1/HMDB_result/metabolomics_crest1_hmdb")

metabolomics_crest1_kegg <-
  metabolomics_crest1_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust)

metabolomics_crest1_hmdb <-
  metabolomics_crest1_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust)

result_all_crest1 <-
  rbind(metabolomics_crest1_kegg,
        metabolomics_crest1_hmdb)

load("metabolomics_pathway_crest2/KEGG_result/metabolomics_crest2_kegg")
load("metabolomics_pathway_crest2/HMDB_result/metabolomics_crest2_hmdb")

metabolomics_crest2_kegg <-
  metabolomics_crest2_kegg@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust)

metabolomics_crest2_hmdb <-
  metabolomics_crest2_hmdb@result %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  dplyr::filter(mapped_number >= 3) %>%
  dplyr::arrange(p_value_adjust)

result_all_crest2 <-
  rbind(metabolomics_crest2_kegg,
        metabolomics_crest2_hmdb)

dim(result_all_crest1)
dim(result_all_crest2)

intersect(result_all_crest1$pathway_name,
          result_all_crest2$pathway_name)

temp_crest1 <-
  result_all_crest1 %>%
  head(20)

temp_crest2 <-
  result_all_crest2 %>%
  head(20)

intersect(temp_crest1$pathway_name,
          temp_crest2$pathway_name)

setdiff(temp_crest1$pathway_name,
        temp_crest2$pathway_name)

setdiff(temp_crest2$pathway_name,
        temp_crest1$pathway_name)

library(ggVennDiagram)
plot <- 
ggVennDiagram(x = list(temp_crest1$pathway_name,
                       temp_crest2$pathway_name))

plot
dir.create("metabolomics_summary")
# ggsave(plot, filename = "metabolomics_summary/venn_plot.pdf", width = 7, height = 7)

library(openxlsx)

# openxlsx::write.xlsx(x = temp_crest1, file = "metabolomics_summary/temp_crest1.xlsx", asTable = TRUE)
# openxlsx::write.xlsx(x = temp_crest2, file = "metabolomics_summary/temp_crest2.xlsx", asTable = TRUE)
