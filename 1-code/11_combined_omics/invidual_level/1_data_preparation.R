no_function()

###
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1-code/100-tools.R")
library(tidyverse)
library(tidymass)

library(microbiomedataset)

###load("data)
load("3-data_analysis/plasma_transcriptome/data_preparation/object")
transcriptome_object <- object

load("3-data_analysis/plasma_proteomics/data_preparation/object")
proteomics_object <- object

load("3-data_analysis/plasma_metabolomics/data_preparation/metabolite/object")
metabolomics_object <- object

load("3-data_analysis/plasma_cytokine/data_preparation/object")
cytokine_object <- object

load("3-data_analysis/clinical_test/data_preparation/object")
clinical_test_object <- object

load("3-data_analysis/plasma_lipidomics/data_preparation/object")
lipidomics_object <- object

load("3-data_analysis/gut_microbiome/data_preparation/object")
gut_microbiome_object <- object

load("3-data_analysis/skin_microbiome/data_preparation/object")
skin_microbiome_object <- object

load("3-data_analysis/oral_microbiome/data_preparation/object")
oral_microbiome_object <- object

load("3-data_analysis/nasal_microbiome/data_preparation/object")
nasal_microbiome_object <- object

dir.create("3-data_analysis/combined_omics/individual/",
           recursive = TRUE)

setwd("3-data_analysis/combined_omics/individual")

dim(transcriptome_object)
dim(proteomics_object)
dim(metabolomics_object)
dim(cytokine_object)
dim(clinical_test_object)
dim(lipidomics_object)
dim(gut_microbiome_object)
dim(oral_microbiome_object)
dim(skin_microbiome_object)
dim(nasal_microbiome_object)

c(
  transcriptome_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  proteomics_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  metabolomics_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  cytokine_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  clinical_test_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  lipidomics_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  gut_microbiome_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  oral_microbiome_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  skin_microbiome_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date),
  nasal_microbiome_object %>%
    activate_microbiome_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == "69-001") %>%
    pull(collection_date)
) %>%
  unique() %>%
  sort()

gut_microbiome_object <-
  gut_microbiome_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("gut_microbiome_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48) 

oral_microbiome_object <-
  oral_microbiome_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("oral_microbiome_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48)

skin_microbiome_object <-
  skin_microbiome_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("skin_microbiome_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48) 

nasal_microbiome_object <-
  nasal_microbiome_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("nasal_microbiome_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48)

proteomics_object <-
  proteomics_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("proteomics_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48)

transcriptome_object <-
  transcriptome_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("transcriptome_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(GENETYPE)) %>%
  dplyr::filter(GENETYPE %in% c("GENETYPE", "protein-coding"))

metabolomics_object <-
  metabolomics_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("metabolomics_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48)

cytokine_object <-
  cytokine_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("cytokine_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48)

lipidomics_object <-
  lipidomics_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("lipidomics_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48)

clinical_test_object <-
  clinical_test_object %>%
  activate_microbiome_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste0("clinical_test_", variable_id)) %>%
  activate_microbiome_dataset(what = "sample_info") %>%
  dplyr::filter(subject_id == "69-001") %>%
  dplyr::filter(CL4 == "Healthy") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(real_age = as.numeric(collection_date - as.Date("04/02/2010", "%m/%d/%Y")) /
                  365 + 59.48)


transcriptome_object@expression_data <-
  transcriptome_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

proteomics_object@expression_data <-
  proteomics_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x))
  }) %>%
  t() %>%
  as.data.frame()

metabolomics_object@expression_data <-
  metabolomics_object@expression_data %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

cytokine_object@expression_data <-
  cytokine_object@expression_data %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

clinical_test_object@expression_data <-
  clinical_test_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

lipidomics_object@expression_data <-
  lipidomics_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x))
  }) %>%
  t() %>%
  as.data.frame()

gut_microbiome_object@expression_data <-
  gut_microbiome_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

skin_microbiome_object@expression_data <-
  skin_microbiome_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

oral_microbiome_object@expression_data <-
  oral_microbiome_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

nasal_microbiome_object@expression_data <-
  nasal_microbiome_object@expression_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()


colnames(proteomics_object)
colnames(transcriptome_object)
colnames(metabolomics_object)
colnames(cytokine_object)
colnames(clinical_test_object)
colnames(lipidomics_object)
colnames(gut_microbiome_object)
colnames(oral_microbiome_object)

dir.create("data_preparation")

save(transcriptome_object, file = "data_preparation/transcriptome_object")
save(proteomics_object, file = "data_preparation/proteomics_object")
save(metabolomics_object, file = "data_preparation/metabolomics_object")
save(cytokine_object, file = "data_preparation/cytokine_object")
save(clinical_test_object, file = "data_preparation/clinical_test_object")
save(lipidomics_object, file = "data_preparation/lipidomics_object")
save(gut_microbiome_object, file = "data_preparation/gut_microbiome_object")
save(oral_microbiome_object, file = "data_preparation/oral_microbiome_object")
save(skin_microbiome_object, file = "data_preparation/skin_microbiome_object")
save(nasal_microbiome_object, file = "data_preparation/nasal_microbiome_object")

# dir.create("3-data_analysis/combined_omics/individual/DE_SWAN")
# setwd("3-data_analysis/combined_omics/individual/DE_SWAN")
#
#
# ##loess fit
# ####fore each molecule, we just use the loess method to
# ####impute the space between the first sample and the last sample
# expression_data <-
#   extract_expression_data(object)
#
# sample_info <-
#   extract_sample_info(object)
#
# variable_info <-
#   extract_variable_info(object)
#
# dim(expression_data)
#
# library(plyr)
#
# new_expression_data <-
#   vector(mode = "list", length = nrow(variable_info))
#
# for (i in seq_along(variable_info$variable_id)) {
#   temp_variable_id <- variable_info$variable_id[i]
#   cat(i, " ")
#   temp_data <-
#     data.frame(value = as.numeric(expression_data[temp_variable_id, ]),
#                sample_info)
#
#   optimize_span <-
#     optimize_loess_span(
#       x = temp_data$real_age,
#       y = temp_data$value,
#       span_range = c(0.3, 0.4, 0.5, 0.6)
#     )
#
#   span <-
#     optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
#
#   value <- temp_data$value
#   real_age <- temp_data$real_age
#
#   ls_reg <-
#     loess(value ~ real_age,
#           span = span)
#
#   prediction_value =
#     predict(ls_reg,
#             newdata = data.frame(real_age = seq(55.5, 62, by = 0.05)))
#   new_expression_data[[i]] <- as.numeric(prediction_value)
# }
#
# new_expression_data <-
#   new_expression_data %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# colnames(new_expression_data) <-
#   as.character(seq(55.5, 62, by = 0.05))
#
# rownames(new_expression_data) <- rownames(expression_data)
# save(new_expression_data, file = "new_expression_data")
#
# load("new_expression_data")
#
# new_expression_data <-
#   new_expression_data %>%
#   dplyr::select(-c(`55.5`, `55.55`))
#
# negative_idx <-
#   new_expression_data %>%
#   apply(1, function(x) {
#     sum(x < 0)
#   }) %>%
#   `>`(0) %>%
#   which()
#
# plot(as.numeric(new_expression_data[negative_idx[1], ]))
#
# new_expression_data[which(new_expression_data < 0, arr.ind = TRUE)] <-
#   0
#
# new_sample_info <-
#   data.frame(
#     sample_id =
#       colnames(new_expression_data),
#     class = "Subject",
#     subject_id = "Subject"
#   )
#
# new_variable_info <-
#   variable_info
#
# save(new_sample_info, file = "new_sample_info")
# save(new_variable_info, file = "new_variable_info")
#
# library(massdataset)
#
# object_loess <-
#   object
#
# object_loess@expression_data <-
#   new_expression_data
#
# object_loess@sample_info <-
#   new_sample_info
#
# new_sample_info_note <-
#   data.frame(
#     name = c("sample_id", "class", "subject_id"),
#     meaning = c("sample_id", "class", "subject_id")
#   )
#
# object_loess@sample_info_note <-
#   new_sample_info_note
#
# colnames(object_loess)
#
# object_loess <-
#   object_loess %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(!sample_id %in% c("55.5", "55.55"))
#
# save(object_loess, file = "object_loess")
#
# sum(is.na(object_loess))
# sum(object_loess < 0)
