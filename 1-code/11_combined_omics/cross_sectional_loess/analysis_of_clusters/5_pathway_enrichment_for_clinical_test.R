no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

1:11 %>%
  purrr::walk(function(i) {
    message(i)
    name <-
      paste0(
        "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_",
        i
      )
    dir.create(name,
               recursive = TRUE)
  })

setwd("3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess")

object_cross_section_loess

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

load("cor_data")

load("final_cluster_info")

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

library(org.Hs.eg.db)
library(clusterProfiler)

####clinical_test
for (i in 1:11) {
  path <- paste0("cluster_", i)
  dir.create(file.path(path, "clinical_test_pathway"))
  message(i)
  cluster_info <-
    final_cluster_info %>%
    dplyr::filter(cluster == i)
  
  ####clinical_test
  temp_data_clinical_test <-
    cluster_info %>%
    dplyr::filter(class == "clinical_test")
  
  if (nrow(temp_data_clinical_test) == 0) {
    next()
  }
  
  library(openxlsx)
  wb <- createWorkbook()
  modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
  addWorksheet(wb, sheetName = "clinical_test")
  freezePane(
    wb = wb,
    sheet = 1,
    firstRow = TRUE,
    firstCol = TRUE
  )
  
  writeDataTable(
    wb,
    sheet = 1,
    x = temp_data_clinical_test %>%
      dplyr::select(-Class),
    colNames = TRUE,
    rowNames = FALSE
  )
  
  saveWorkbook(
    wb,
    file = file.path(path,
                     "clinical_test_pathway/clinical_test.xlsx"),
    overwrite = TRUE
  )
  
}
