no_source()

setwd(r4projects::get_project_wd())
rm(list = ls())
source("1-code/100-tools.R")

library(tidyverse)

load("3-data_analysis/phenotype/data_preparation/phenotype_data")

load("3-data_analysis/plasma_cytokine/data_preparation/object")
cytokine_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

cytokine_object <-
  object

load("3-data_analysis/plasma_lipidomics/data_preparation/object")
lipidomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

library(plyr)
lipidomics_sample_info %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    x %>%
      pull(BMI)
  })

lipidomics_sample_info %>%
  plyr::dlply(.variables = .(subject_id)) %>%
  purrr::map(function(x) {
    x %>%
      pull(SSPG)
  })

lipidomics_object <-
  object

load("3-data_analysis/plasma_metabolomics/data_preparation/metabolite/object")
metabolomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

metabolomics_object <-
  object

load("3-data_analysis/plasma_proteomics/data_preparation/object")
proteomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

proteomics_object <-
  object

load("3-data_analysis/plasma_transcriptome/data_preparation/object")
transcriptome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

transcriptome_object <-
  object

load("3-data_analysis/clinical_test/data_preparation/object")
clinical_test_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

clinical_test_object <-
  object

load("3-data_analysis/gut_microbiome/data_preparation/object")
gut_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

gut_microbiome_object <-
  object

load("3-data_analysis/nasal_microbiome/data_preparation/object")
nasal_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

nasal_microbiome_object <-
  object

load("3-data_analysis/skin_microbiome/data_preparation/object")
skin_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

skin_microbiome_object <-
  object

load("3-data_analysis/oral_microbiome/data_preparation/object")
oral_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

oral_microbiome_object <-
  object

dir.create("3-data_analysis/sample_info/")
setwd("3-data_analysis/sample_info")

dim(transcriptome_object)
dim(proteomics_object)
dim(metabolomics_object)
dim(cytokine_object)
dim(clinical_test_object)
dim(lipidomics_object)
dim(gut_microbiome_object)
unique(gut_microbiome_object@variable_info$Genus) %>%
  length

unique(skin_microbiome_object@variable_info$Genus) %>%
  length

unique(oral_microbiome_object@variable_info$Genus) %>%
  length

unique(nasal_microbiome_object@variable_info$Genus) %>%
  length

all_subject_id <-
  c(
    transcriptome_sample_info$subject_id,
    proteomics_sample_info$subject_id,
    metabolomics_sample_info$subject_id,
    cytokine_sample_info$subject_id,
    clinical_test_sample_info$subject_id,
    lipidomics_sample_info$subject_id,
    gut_microbiome_sample_info$subject_id,
    skin_microbiome_sample_info$subject_id,
    oral_microbiome_sample_info$subject_id,
    nasal_microbiome_sample_info$subject_id
  ) %>%
  unique()

library(plyr)

sample_number <-
  purrr::map2(.x = list(
    cytokine_sample_info,
    lipidomics_sample_info,
    metabolomics_sample_info,
    proteomics_sample_info,
    transcriptome_sample_info,
    clinical_test_sample_info,
    gut_microbiome_sample_info,
    nasal_microbiome_sample_info,
    skin_microbiome_sample_info,
    oral_microbiome_sample_info
  ),
  .y = c(
    "cytokine",
    "lipidomics",
    "metabolomics",
    "proteomics",
    "transcriptome",
    "clinical_test",
    "gut_microbiome",
    "nasal_microbiome",
    "skin_microbiome",
    "oral_microbiome"
  ),
  function(x, y) {
    number = x %>%
      dplyr::count(subject_id,
                   name = paste0(y, "_number"))
    
    days <-
      x %>%
      plyr::dlply(.variables = .(subject_id)) %>%
      purrr::map(function(z) {
        z <-
          z %>%
          arrange(collection_date)
        data.frame(
          subject_id = z$subject_id[1],
          day_range = as.numeric(z$collection_date[nrow(z)] - z$collection_date[1])
        )
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()
    
    colnames(days)[2] <-
      paste0(y, "_day_range")
    
    number %>%
      dplyr::left_join(days, by = "subject_id")
  })

sample_number <-
  sample_number[[1]] %>%
  dplyr::full_join(sample_number[[2]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[3]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[4]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[5]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[6]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[7]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[8]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[9]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[10]], by = c("subject_id"))

transcriptome_subject_info <-
  transcriptome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

proteomics_subject_info <-
  proteomics_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

metabolomics_subject_info <-
  metabolomics_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

cytokine_subject_info <-
  cytokine_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

lipidomics_subject_info <-
  lipidomics_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

clinical_test_subject_info <-
  clinical_test_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

gut_microbiome_subject_info <-
  gut_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

skin_microbiome_subject_info <-
  skin_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

oral_microbiome_subject_info <-
  oral_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

nasal_microbiome_subject_info <-
  nasal_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

subject_info <-
  transcriptome_subject_info %>%
  dplyr::full_join(proteomics_subject_info, by = "subject_id") %>%
  dplyr::full_join(metabolomics_subject_info, by = "subject_id") %>%
  dplyr::full_join(lipidomics_subject_info, by = "subject_id") %>%
  dplyr::full_join(clinical_test_subject_info, by = "subject_id") %>%
  dplyr::full_join(cytokine_subject_info, by = "subject_id") %>%
  dplyr::full_join(gut_microbiome_subject_info, by = "subject_id") %>%
  dplyr::full_join(oral_microbiome_subject_info, by = "subject_id") %>%
  dplyr::full_join(nasal_microbiome_subject_info, by = "subject_id") %>%
  dplyr::full_join(skin_microbiome_subject_info, by = "subject_id")

subject_id_random <-
  subject_info %>%
  dplyr::select(matches("subject_id_random")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

consented <-
  subject_info %>%
  dplyr::select(matches("consented")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

IRIS <-
  subject_info %>%
  dplyr::select(matches("IRIS")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

SSPG <-
  subject_info %>%
  dplyr::select(matches("SSPG")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

FPG <-
  subject_info %>%
  dplyr::select(matches("FPG")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

diabetes_class <-
  subject_info %>%
  dplyr::select(matches("diabetes_class")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

Gender <-
  subject_info %>%
  dplyr::select(matches("Gender")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

Ethnicity <-
  subject_info %>%
  dplyr::select(matches("Ethnicity")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

adjusted_age <-
  subject_info %>%
  dplyr::select(matches("adjusted_age")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

BMI <-
  subject_info %>%
  dplyr::select(matches("BMI")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

subject_info <-
  data.frame(
    subject_id = subject_info$subject_id,
    subject_id_random = subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

subject_info <-
  subject_info %>%
  dplyr::left_join(sample_number, by = "subject_id")


subject_info$subject_id_random <-
  masstools::name_duplicated(subject_info$subject_id_random)


dim(subject_info)

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "subject_info")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = subject_info,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# saveWorkbook(
#   wb,
#   "subject_info.xlsx",
#   overwrite = TRUE
# )

which(clinical_test_object@variable_info$variable_id == "A1C")

temp <-
clinical_test_object %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  dplyr::filter(variable_id == "A1C") %>% 
  pivot_longer() %>% 
  dplyr::left_join(clinical_test_object@sample_info[,c("subject_id_random", "sample_id")],
                   by = "sample_id") %>% 
  plyr::dlply(.variables = .(subject_id_random)) %>% 
  purrr::map(function(x){
    sum(x$value > 6.5, na.rm = TRUE)
  }) %>% 
  unlist()

temp <- 
  temp[temp > 0]


names(temp)
  

clinical_test <- 
clinical_test_object %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  dplyr::filter(variable_id == "A1C") %>% 
  pivot_longer() %>% 
  dplyr::left_join(clinical_test_object@sample_info[,c("subject_id_random", "sample_id", "collection_date")],
                   by = "sample_id")

# 
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "clinical_test")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
# 
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = clinical_test,
#   colNames = TRUE,
#   rowNames = FALSE
# )
# 
# saveWorkbook(
#   wb,
#   "clinical_test.xlsx",
#   overwrite = TRUE
# )
