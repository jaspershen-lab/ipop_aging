no_function()

###
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1-code/100-tools.R")
library(tidyverse)
library(tidymass)

load(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

setwd("3-data_analysis/combined_omics/individual/")

load("data_preparation/transcriptome_object")
load("data_preparation/proteomics_object")
load("data_preparation/metabolomics_object")
load("data_preparation/lipidomics_object")
load("data_preparation/cytokine_object")
load("data_preparation/clinical_test_object")
load("data_preparation/gut_microbiome_object")
load("data_preparation/skin_microbiome_object")
load("data_preparation/oral_microbiome_object")
load("data_preparation/nasal_microbiome_object")

dir.create("scatter_plot")
setwd("scatter_plot/")

dim(object)

dim(final_cluster_info)


###cluster 2, 4, and 5
temp_cluster_info <-
  final_cluster_info %>%
  dplyr::filter(cluster %in% c("2", "4", "5"))

temp_data <-
  rbind(
    transcriptome_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(transcriptome_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "transcriptome"),
    proteomics_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(proteomics_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "proteomics"),
    metabolomics_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(metabolomics_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "metabolomics"),
    lipidomics_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(lipidomics_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "lipidomics"),
    cytokine_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(cytokine_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "cytokine"),
    clinical_test_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(clinical_test_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "clinical_test"),
    skin_microbiome_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(skin_microbiome_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "skin_microbiome"),
    gut_microbiome_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(gut_microbiome_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "gut_microbiome"),
    oral_microbiome_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(oral_microbiome_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "oral_microbiome"),
    nasal_microbiome_object %>%
      activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(variable_id %in% temp_cluster_info$variable_id) %>%
      pivot_longer() %>%
      dplyr::left_join(nasal_microbiome_object@sample_info[, c("sample_id", "real_age")],
                       by = "sample_id") %>%
      dplyr::mutate(class = "nasal_microbiome")
  ) %>%
  dplyr::left_join(final_cluster_info[, c("variable_id", "cluster")],
                   by = "variable_id")

plot <-
temp_data %>%
  ggplot(aes(x = real_age, y = value)) +
  geom_line(aes(group = variable_id, color = class)) +
  facet_wrap(facets = vars(cluster), scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.justification = c(0, 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(
    x = "Age (years)",
    y = "Z-score"
  ) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = omics_color)

ggsave(plot, filename = "omics_cluster_age.pdf", width = 12, height = 6)
