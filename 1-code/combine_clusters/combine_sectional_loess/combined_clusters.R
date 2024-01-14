no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
###transcriptome
load(
  "3-data_analysis/plasma_transcriptome/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)
transcriptome_cluster_info <-
  final_cluster_info

transcriptome_data <-
  readr::read_delim(
    "3-data_analysis/plasma_transcriptome/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  transcriptome_data$...1[-1]
transcriptome_data <-
  transcriptome_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

transcriptome_data <-
  transcriptome_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(transcriptome_data) <-
  row_names

load("3-data_analysis/plasma_transcriptome/data_preparation/object_cross_section")

transcriptome_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

###proteomics
load(
  "3-data_analysis/plasma_proteomics/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)
proteomics_cluster_info <-
  final_cluster_info

proteomics_data <-
  readr::read_delim(
    "3-data_analysis/plasma_proteomics/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  proteomics_data$...1[-1]
proteomics_data <-
  proteomics_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

proteomics_data <-
  proteomics_data %>%
  apply(1, function(x) {
    (x - mean(x))
  }) %>%
  t() %>%
  as.data.frame()

rownames(proteomics_data) <-
  row_names

load("3-data_analysis/plasma_proteomics/data_preparation/object_cross_section")

proteomics_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x))
  }) %>%
  t() %>%
  as.data.frame()


###metabolomics
load(
  "3-data_analysis/plasma_metabolomics/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)
metabolomics_cluster_info <-
  final_cluster_info

metabolomics_data <-
  readr::read_delim(
    "3-data_analysis/plasma_metabolomics/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  metabolomics_data$...1[-1]
metabolomics_data <-
  metabolomics_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

metabolomics_data <-
  metabolomics_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(metabolomics_data) <-
  row_names

load(
  "3-data_analysis/plasma_metabolomics/data_preparation/metabolite/object_cross_section"
)

metabolomics_all_data <-
  object_cross_section@expression_data %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


###cytokine
load(
  "3-data_analysis/plasma_cytokine/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

cytokine_cluster_info <-
  final_cluster_info

cytokine_data <-
  readr::read_delim(
    "3-data_analysis/plasma_cytokine/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  cytokine_data$...1[-1]
cytokine_data <-
  cytokine_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

cytokine_data <-
  cytokine_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(cytokine_data) <-
  row_names

load("3-data_analysis/plasma_cytokine/data_preparation/object_cross_section")

cytokine_all_data <-
  object_cross_section@expression_data %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

###clinical_test
load(
  "3-data_analysis/clinical_test/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

clinical_test_cluster_info <-
  final_cluster_info

clinical_test_data <-
  readr::read_delim(
    "3-data_analysis/clinical_test/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  clinical_test_data$...1[-1]
clinical_test_data <-
  clinical_test_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

clinical_test_data <-
  clinical_test_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(clinical_test_data) <-
  row_names


load("3-data_analysis/clinical_test/data_preparation/object_cross_section")

clinical_test_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

###lipidomics
load(
  "3-data_analysis/plasma_lipidomics/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

lipidomics_cluster_info <-
  final_cluster_info

lipidomics_data <-
  readr::read_delim(
    "3-data_analysis/plasma_lipidomics/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  lipidomics_data$...1[-1]
lipidomics_data <-
  lipidomics_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

lipidomics_data <-
  lipidomics_data %>%
  apply(1, function(x) {
    (x - mean(x))
  }) %>%
  t() %>%
  as.data.frame()

rownames(lipidomics_data) <-
  row_names



load("3-data_analysis/plasma_lipidomics/data_preparation/object_cross_section")

lipidomics_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x))
  }) %>%
  t() %>%
  as.data.frame()

###gut_microbiome
load(
  "3-data_analysis/gut_microbiome/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

gut_microbiome_cluster_info <-
  final_cluster_info

gut_microbiome_data <-
  readr::read_delim(
    "3-data_analysis/gut_microbiome/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  gut_microbiome_data$...1[-1]
gut_microbiome_data <-
  gut_microbiome_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

gut_microbiome_data <-
  gut_microbiome_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(gut_microbiome_data) <-
  row_names

load("3-data_analysis/gut_microbiome/data_preparation/object_cross_section")

library(microbiomedataset)

object_cross_section <-
  object_cross_section %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

dim(object_cross_section)

non_zero_per <-
  apply(object_cross_section, 1, function(x) {
    sum(x != 0) / ncol(object_cross_section)
  })

idx <-
  which(non_zero_per > 0.1)

object_cross_section <-
  object_cross_section[idx,]

object_cross_section <-
  transform2relative_intensity(object_cross_section)

gut_microbiome_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


###skin_microbiome
load(
  "3-data_analysis/skin_microbiome/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

skin_microbiome_cluster_info <-
  final_cluster_info

skin_microbiome_data <-
  readr::read_delim(
    "3-data_analysis/skin_microbiome/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  skin_microbiome_data$...1[-1]
skin_microbiome_data <-
  skin_microbiome_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

skin_microbiome_data <-
  skin_microbiome_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(skin_microbiome_data) <-
  row_names

load("3-data_analysis/skin_microbiome/data_preparation/object_cross_section")

object_cross_section <-
  object_cross_section %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

dim(object_cross_section)

non_zero_per <-
  apply(object_cross_section, 1, function(x) {
    sum(x != 0) / ncol(object_cross_section)
  })

idx <-
  which(non_zero_per > 0.1)

object_cross_section <-
  object_cross_section[idx,]

object_cross_section <-
  transform2relative_intensity(object_cross_section)

skin_microbiome_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


###oral_microbiome
load(
  "3-data_analysis/oral_microbiome/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

oral_microbiome_cluster_info <-
  final_cluster_info

oral_microbiome_data <-
  readr::read_delim(
    "3-data_analysis/oral_microbiome/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  oral_microbiome_data$...1[-1]
oral_microbiome_data <-
  oral_microbiome_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

oral_microbiome_data <-
  oral_microbiome_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(oral_microbiome_data) <-
  row_names

load("3-data_analysis/oral_microbiome/data_preparation/object_cross_section")

object_cross_section <-
  object_cross_section %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

dim(object_cross_section)

non_zero_per <-
  apply(object_cross_section, 1, function(x) {
    sum(x != 0) / ncol(object_cross_section)
  })

idx <-
  which(non_zero_per > 0.1)

object_cross_section <-
  object_cross_section[idx,]

object_cross_section <-
  transform2relative_intensity(object_cross_section)

oral_microbiome_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

###nasal_microbiome
load(
  "3-data_analysis/nasal_microbiome/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

nasal_microbiome_cluster_info <-
  final_cluster_info

nasal_microbiome_data <-
  readr::read_delim(
    "3-data_analysis/nasal_microbiome/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt"
  )

row_names <-
  nasal_microbiome_data$...1[-1]
nasal_microbiome_data <-
  nasal_microbiome_data[-1, ] %>%
  tibble::column_to_rownames(var = "...1") %>%
  purrr::map(function(x) {
    as.numeric(x)
  }) %>%
  dplyr::bind_cols()

nasal_microbiome_data <-
  nasal_microbiome_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

rownames(nasal_microbiome_data) <-
  row_names

load("3-data_analysis/nasal_microbiome/data_preparation/object_cross_section")

object_cross_section <-
  object_cross_section %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

dim(object_cross_section)

non_zero_per <-
  apply(object_cross_section, 1, function(x) {
    sum(x != 0) / ncol(object_cross_section)
  })

idx <-
  which(non_zero_per > 0.1)

object_cross_section <-
  object_cross_section[idx,]

object_cross_section <-
  transform2relative_intensity(object_cross_section)

nasal_microbiome_all_data <-
  object_cross_section@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

dir.create("3-data_analysis/combined_clusters/cross_section_loess/",
           recursive = TRUE)

setwd("3-data_analysis/combined_clusters/cross_section_loess/")

temp <-
  transcriptome_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

transcriptome_cluster_info <-
  transcriptome_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

transcriptome_mean_data <-
  unique(transcriptome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(transcriptome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    transcriptome_data[transcriptome_cluster_info$variable_id[transcriptome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(transcriptome_mean_data) <-
  paste("transcriptome",
        unique(transcriptome_cluster_info$cluster),
        sep = "_")

transcriptome_mean_all_data <-
  unique(transcriptome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(transcriptome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    transcriptome_all_data[transcriptome_cluster_info$variable_id[transcriptome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(transcriptome_mean_all_data) <-
  paste("transcriptome",
        unique(transcriptome_cluster_info$cluster),
        sep = "_")

temp <-
  proteomics_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

proteomics_cluster_info <-
  proteomics_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

proteomics_mean_data <-
  unique(proteomics_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(proteomics_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    proteomics_data[proteomics_cluster_info$variable_id[proteomics_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

rownames(proteomics_mean_data) <-
  paste("proteomics", unique(proteomics_cluster_info$cluster), sep = "_")

proteomics_mean_all_data <-
  unique(proteomics_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(proteomics_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    proteomics_all_data[proteomics_cluster_info$variable_id[proteomics_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(proteomics_mean_all_data) <-
  paste("proteomics",
        unique(proteomics_cluster_info$cluster),
        sep = "_")


temp <-
  metabolomics_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

metabolomics_cluster_info <-
  metabolomics_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

metabolomics_mean_data <-
  unique(metabolomics_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(metabolomics_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    metabolomics_data[metabolomics_cluster_info$variable_id[metabolomics_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(metabolomics_mean_data) <-
  paste("metabolomics",
        unique(metabolomics_cluster_info$cluster),
        sep = "_")



metabolomics_mean_all_data <-
  unique(metabolomics_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(metabolomics_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    metabolomics_all_data[metabolomics_cluster_info$variable_id[metabolomics_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(metabolomics_mean_all_data) <-
  paste("metabolomics",
        unique(metabolomics_cluster_info$cluster),
        sep = "_")

temp <-
  cytokine_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

cytokine_cluster_info <-
  cytokine_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

cytokine_mean_data <-
  unique(cytokine_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(cytokine_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    cytokine_data[cytokine_cluster_info$variable_id[cytokine_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(cytokine_mean_data) <-
  paste("cytokine", unique(cytokine_cluster_info$cluster), sep = "_")

cytokine_mean_all_data <-
  unique(cytokine_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(cytokine_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    cytokine_all_data[cytokine_cluster_info$variable_id[cytokine_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(cytokine_mean_all_data) <-
  paste("cytokine",
        unique(cytokine_cluster_info$cluster),
        sep = "_")

temp <-
  clinical_test_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

clinical_test_cluster_info <-
  clinical_test_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

clinical_test_mean_data <-
  unique(clinical_test_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(clinical_test_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    clinical_test_data[clinical_test_cluster_info$variable_id[clinical_test_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

rownames(clinical_test_mean_data) <-
  paste("clinical_test",
        unique(clinical_test_cluster_info$cluster),
        sep = "_")



clinical_test_mean_all_data <-
  unique(clinical_test_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(clinical_test_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    clinical_test_all_data[clinical_test_cluster_info$variable_id[clinical_test_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(clinical_test_mean_all_data) <-
  paste("clinical_test",
        unique(clinical_test_cluster_info$cluster),
        sep = "_")


temp <-
  lipidomics_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

lipidomics_cluster_info <-
  lipidomics_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

lipidomics_mean_data <-
  unique(lipidomics_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(lipidomics_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    lipidomics_data[lipidomics_cluster_info$variable_id[lipidomics_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(lipidomics_mean_data) <-
  paste("lipidomics", unique(lipidomics_cluster_info$cluster), sep = "_")


lipidomics_mean_all_data <-
  unique(lipidomics_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(lipidomics_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    lipidomics_all_data[lipidomics_cluster_info$variable_id[lipidomics_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(lipidomics_mean_all_data) <-
  paste("lipidomics",
        unique(lipidomics_cluster_info$cluster),
        sep = "_")

temp <-
  gut_microbiome_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

gut_microbiome_cluster_info <-
  gut_microbiome_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

gut_microbiome_mean_data <-
  unique(gut_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(gut_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    gut_microbiome_data[gut_microbiome_cluster_info$variable_id[gut_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(gut_microbiome_mean_data) <-
  paste("gut_microbiome",
        unique(gut_microbiome_cluster_info$cluster),
        sep = "_")


gut_microbiome_mean_all_data <-
  unique(gut_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(gut_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    gut_microbiome_all_data[gut_microbiome_cluster_info$variable_id[gut_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(gut_microbiome_mean_all_data) <-
  paste("gut_microbiome",
        unique(gut_microbiome_cluster_info$cluster),
        sep = "_")

temp <-
  skin_microbiome_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

skin_microbiome_cluster_info <-
  skin_microbiome_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)


skin_microbiome_mean_data <-
  unique(skin_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(skin_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    skin_microbiome_data[skin_microbiome_cluster_info$variable_id[skin_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(skin_microbiome_mean_data) <-
  paste("skin_microbiome",
        unique(skin_microbiome_cluster_info$cluster),
        sep = "_")

skin_microbiome_mean_all_data <-
  unique(skin_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(skin_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    skin_microbiome_all_data[skin_microbiome_cluster_info$variable_id[skin_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(skin_microbiome_mean_all_data) <-
  paste("skin_microbiome",
        unique(skin_microbiome_cluster_info$cluster),
        sep = "_")


temp <-
  oral_microbiome_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

oral_microbiome_cluster_info <-
  oral_microbiome_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

oral_microbiome_mean_data <-
  unique(oral_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(oral_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    oral_microbiome_data[oral_microbiome_cluster_info$variable_id[oral_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(oral_microbiome_mean_data) <-
  paste("oral_microbiome",
        unique(oral_microbiome_cluster_info$cluster),
        sep = "_")


oral_microbiome_mean_all_data <-
  unique(oral_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(oral_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    oral_microbiome_all_data[oral_microbiome_cluster_info$variable_id[oral_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(oral_microbiome_mean_all_data) <-
  paste("oral_microbiome",
        unique(oral_microbiome_cluster_info$cluster),
        sep = "_")

temp <-
  nasal_microbiome_cluster_info %>%
  dplyr::count(cluster) %>%
  dplyr::filter(n > 1)

nasal_microbiome_cluster_info <-
  nasal_microbiome_cluster_info %>%
  dplyr::filter(cluster %in% temp$cluster)

nasal_microbiome_mean_data <-
  unique(nasal_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(nasal_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    nasal_microbiome_data[nasal_microbiome_cluster_info$variable_id[nasal_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(nasal_microbiome_mean_data) <-
  paste("nasal_microbiome",
        unique(nasal_microbiome_cluster_info$cluster),
        sep = "_")


nasal_microbiome_mean_all_data <-
  unique(nasal_microbiome_cluster_info$cluster) %>%
  purrr::map(function(idx) {
    if (sum(nasal_microbiome_cluster_info$cluster == idx) == 1) {
      return(NULL)
    }
    nasal_microbiome_all_data[nasal_microbiome_cluster_info$variable_id[nasal_microbiome_cluster_info$cluster == idx],] %>%
      colMeans()
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()
rownames(nasal_microbiome_mean_all_data) <-
  paste("nasal_microbiome",
        unique(nasal_microbiome_cluster_info$cluster),
        sep = "_")

temp_data <-
  rbind(
    transcriptome_mean_data,
    proteomics_mean_data,
    metabolomics_mean_data,
    cytokine_mean_data,
    clinical_test_mean_data,
    lipidomics_mean_data
  )

intersect_sample_id <-
  reduce(
    .x = list(
      colnames(transcriptome_mean_all_data),
      colnames(proteomics_mean_all_data),
      colnames(metabolomics_mean_all_data),
      colnames(cytokine_mean_all_data),
      colnames(clinical_test_mean_all_data),
      colnames(lipidomics_mean_all_data),
      colnames(gut_microbiome_mean_all_data),
      colnames(skin_microbiome_mean_all_data),
      colnames(oral_microbiome_mean_all_data),
      colnames(nasal_microbiome_mean_all_data)
    ),
    .f = intersect
  )

sample_info <-
  object_cross_section@sample_info %>%
  dplyr::arrange(adjusted_age) %>%
  dplyr::filter(sample_id %in% intersect_sample_id)

intersect_sample_id <-
  sample_info$sample_id

temp_data_all <-
  rbind(
    transcriptome_mean_all_data[, intersect_sample_id],
    proteomics_mean_all_data[, intersect_sample_id],
    metabolomics_mean_all_data[, intersect_sample_id],
    cytokine_mean_all_data[, intersect_sample_id],
    clinical_test_mean_all_data[, intersect_sample_id],
    lipidomics_mean_all_data[, intersect_sample_id],
    gut_microbiome_mean_all_data[, intersect_sample_id],
    skin_microbiome_mean_all_data[, intersect_sample_id],
    oral_microbiome_mean_all_data[, intersect_sample_id],
    nasal_microbiome_mean_all_data[, intersect_sample_id]
  )

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2),
                     c(
                       viridis::viridis(n = 3)[1],
                       viridis::viridis(n = 3)[2],
                       viridis::viridis(n = 3)[3]
                     ))


ha <- rowAnnotation(
  class = rownames(temp_data) %>% stringr::str_replace("_[0-9]{1,2}", ""),
  col = list(class = omics_color),
  na_col = "black",
  border = TRUE,
  gp = gpar(col = "black")
)

plot <-
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    border = TRUE,
    column_names_rot = 45,
    rect_gp = gpar(col = "black"),
    col = col_fun,
    name = "Z-score",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    right_annotation = ha
  )

plot <-
  ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "cluster_heatmap.pdf",
       width = 5,
       height = 12)


# ha <- rowAnnotation(
#   class = rownames(temp_data_all) %>% stringr::str_replace("_[0-9]{1,2}", ""),
#   col = list(class = omics_color),
#   na_col = "black",
#   border = TRUE,
#   gp = gpar(col = "black")
# )
# 
# plot <-
#   Heatmap(
#     temp_data_all,
#     cluster_columns = FALSE,
#     border = TRUE,
#     column_names_rot = 45,
#     rect_gp = gpar(col = "black"),
#     col = col_fun,
#     name = "Z-score",
#     clustering_distance_columns = "euclidean",
#     clustering_method_columns = "complete",
#     right_annotation = ha
#   )
# 
# plot <-
#   ggplotify::as.ggplot(plot)
# plot
# ggsave(plot,
#        filename = "cluster_heatmap.pdf",
#        width = 5,
#        height = 12)
