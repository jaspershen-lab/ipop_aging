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

load("3-data_analysis/clinical_test/data_preparation/object_corss_section")

setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_5/clinical_test_pathway/"
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

library(org.Hs.eg.db)
library(clusterProfiler)

###cluster 5
cluster5_clinical_test <-
  final_cluster_info %>%
  dplyr::filter(cluster == 5 & class == "clinical_test")

dim(cluster5_clinical_test)

cluster5_clinical_test$mol_name
cluster5_clinical_test$variable_id
cluster5_clinical_test$test_name

plot <-
  object_cross_section_loess %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "clinical_test_CR") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

plot

ggsave(plot,
       file = "clinical_test_CR.pdf",
       width = 8,
       height = 6)




plot <-
  object_cross_section_loess %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "clinical_test_MCH") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

plot

ggsave(plot,
       file = "clinical_test_MCH.pdf",
       width = 8,
       height = 6)

plot <-
  object_cross_section_loess %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "clinical_test_MONO") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

plot

ggsave(plot,
       file = "clinical_test_MONO.pdf",
       width = 8,
       height = 6)



plot <-
  object_cross_section_loess %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "clinical_test_RDW") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

plot

ggsave(plot,
       file = "clinical_test_RDW.pdf",
       width = 8,
       height = 6)



object_corss_section %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "CR") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = adjusted_age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")


object_corss_section %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "MCH") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = adjusted_age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

object_corss_section %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "MONO") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = adjusted_age, n = 10))) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)")

object_corss_section %>%
  ggplot_mass_dataset(direction = "variable", variable_id = "RDW") +
  geom_boxplot(aes(x = ggplot2::cut_interval(x = adjusted_age, n = 10)),
               outlier.shape = NA) +
  theme_base +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  labs(x = "Age range (years)") +
  scale_y_continuous(limits = c(12.5, 17.5))
