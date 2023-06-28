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

###load marker
load("data_analysis/plasma_transcriptome/DEG/cross_section_loess/all_marker_name")
transcriptome_marker <-
  paste("transcriptome", all_marker_name, sep = "_")

load("data_analysis/plasma_proteomics/DEG/cross_section_loess/all_marker_name")
proteomics_marker <-
  paste("protein", all_marker_name, sep = "_")

load("data_analysis/plasma_metabolomics/DEG/cross_section_loess/all_marker_name")
metabolomics_marker <- all_marker_name

load("data_analysis/plasma_cytokine/DEG/cross_section_loess/all_marker_name")
cytokine_marker <- all_marker_name

load("data_analysis/plasma_lipidomics/DEG/cross_section_loess/all_marker_name")
lipidomics_marker <- all_marker_name

load("data_analysis/clinical_test/DEG/cross_section_loess/all_marker_name")
clinical_test_marker <- all_marker_name

load("data_analysis/gut_microbiome/DEG/cross_section_loess/all_marker_name")
gut_microbiome_marker <-
  paste("gut", all_marker_name, sep = "_")

load("data_analysis/skin_microbiome/DEG/cross_section_loess/all_marker_name")
skin_microbiome_marker <-
  paste("skin", all_marker_name, sep = "_")

load("data_analysis/oral_microbiome/DEG/cross_section_loess/all_marker_name")
oral_microbiome_marker <-
  paste("oral", all_marker_name, sep = "_")

load("data_analysis/nasal_microbiome/DEG/cross_section_loess/all_marker_name")
nasal_microbiome_marker <-
  paste("nasal", all_marker_name, sep = "_")

dir.create(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/",
  recursive = TRUE
)

setwd("data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess")

match(transcriptome_marker,
      rownames(object_cross_section_loess))

match(proteomics_marker,
      rownames(object_cross_section_loess))

match(metabolomics_marker,
      rownames(object_cross_section_loess))

match(cytokine_marker,
      rownames(object_cross_section_loess))

match(lipidomics_marker,
      rownames(object_cross_section_loess))

match(clinical_test_marker,
      rownames(object_cross_section_loess))

match(gut_microbiome_marker,
      rownames(object_cross_section_loess))

match(skin_microbiome_marker,
      rownames(object_cross_section_loess))

match(oral_microbiome_marker,
      rownames(object_cross_section_loess))

match(nasal_microbiome_marker,
      rownames(object_cross_section_loess))


marker <-
  c(
    transcriptome_marker,
    proteomics_marker,
    metabolomics_marker,
    cytokine_marker,
    lipidomics_marker,
    clinical_test_marker,
    gut_microbiome_marker,
    skin_microbiome_marker,
    oral_microbiome_marker,
    nasal_microbiome_marker
  )

object_cross_section_loess <-
  object_cross_section_loess[marker,]

object <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "variable_info") %>%
  massdataset::split_mass_dataset(by = "data_type")

length(object)

####clustering
library(Mfuzz)

######cluster all variables using the age range
object_cross_section_loess@sample_info$adjusted_age %>% range

object_cross_section_loess@sample_info$adjusted_age %>% plot

age_index <-
  data.frame(from = c(25, 40, 45, 50, 55, 60, 65),
             to =   c(40, 45, 50, 55, 60, 65, 75.5))

apply(age_index, 1, function(x) {
  sum(
    object_cross_section_loess@sample_info$adjusted_age > x[1] &
      object_cross_section_loess@sample_info$adjusted_age <= x[2]
  )
})

##transcriptome
temp_data_metabolomics <-
  log(object$metabolomics@expression_data + 1, 2)

temp_data_proteomics <-
  object$proteomics@expression_data

temp_data_transcriptome <-
  object$transcriptome@expression_data

temp_data_cytokine <-
  log(object$cytokine@expression_data + 1, 2)

temp_data_lipidomics <-
  object$lipidomics@expression_data

temp_data_gut_microbiome <-
  object$gut_microbiome@expression_data

temp_data_nasal_microbiome <-
  object$nasal_microbiome@expression_data

temp_data_oral_microbiome <-
  object$oral_microbiome@expression_data

temp_data_skin_microbiome <-
  object$skin_microbiome@expression_data

temp_data_clinical_test <-
  object$clinical_test@expression_data

sample_info <-
  extract_sample_info(object_cross_section_loess) %>%
  dplyr::arrange(adjusted_age)

temp_data <-
  rbind(
    temp_data_metabolomics,
    temp_data_proteomics,
    temp_data_transcriptome,
    temp_data_cytokine,
    temp_data_lipidomics,
    temp_data_gut_microbiome,
    temp_data_nasal_microbiome,
    temp_data_oral_microbiome,
    temp_data_skin_microbiome,
    temp_data_clinical_test
  )[rownames(object_cross_section_loess), sample_info$sample_id]

temp_data[object_cross_section_loess@variable_info$data_type == "metabolomics", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "metabolomics", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "proteomics", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "proteomics", ] %>%
  apply(1, function(x) {
    x - mean(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "transcriptome", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "transcriptome", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "cytokine", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "cytokine", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "lipidomics", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "lipidomics", ] %>%
  apply(1, function(x) {
    x - mean(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "gut_microbiome", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "gut_microbiome", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "nasal_microbiome", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "nasal_microbiome", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "oral_microbiome", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "oral_microbiome", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "skin_microbiome", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "skin_microbiome", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

temp_data[object_cross_section_loess@variable_info$data_type == "clinical_test", ] <-
  temp_data[object_cross_section_loess@variable_info$data_type == "clinical_test", ] %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(circlize)

col_fun = colorRamp2(c(-2, 0, 2),
                     c(
                       viridis::viridis(n = 3)[1],
                       viridis::viridis(n = 3)[2],
                       viridis::viridis(n = 3)[3]
                     ))

plot <-
  Heatmap(
    temp_data,
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    col = col_fun,
    border = TRUE
  )

plot <- ggplotify::as.ggplot(plot)

plot

ggsave(plot,
       file = "heatmap.pdf",
       width = 9,
       height = 7)

ggsave(plot,
       file = "heatmap.png",
       width = 9,
       height = 7)
