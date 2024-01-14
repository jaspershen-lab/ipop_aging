no_function()

###
setwd(r4projects::get_project_wd())
rm(list = ls())
# source("1-code/100-tools.R")
library(tidyverse)
library(tidymass)

load("3-data_analysis/plasma_transcriptome/data_preparation/object")

table(object@variable_info$GENETYPE)

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(real_age)) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(GENETYPE)) %>%
  dplyr::filter(GENETYPE %in% c("GENETYPE", "protein-coding"))

####only remain the health vist
object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(CL4 == "Healthy")

object %>%
  activate_mass_dataset(what = "sample_info") %>%
  count(subject_id) %>%
  arrange(desc(n)) %>%
  head()

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(subject_id == "69-001")

object@sample_info$real_age <-
  as.numeric(object@sample_info$collection_date - as.Date("01/01/1955", "%m/%d/%Y")) /
  365

dir.create("3-data_analysis/combined_omics/individual/DE_SWAN")
setwd("3-data_analysis/combined_omics/individual/DE_SWAN")


##loess fit
####fore each molecule, we just use the loess method to
####impute the space between the first sample and the last sample
expression_data <-
  extract_expression_data(object)

sample_info <-
  extract_sample_info(object)

variable_info <-
  extract_variable_info(object)

dim(expression_data)

library(plyr)

new_expression_data <-
  vector(mode = "list", length = nrow(variable_info))

for (i in seq_along(variable_info$variable_id)) {
  temp_variable_id <- variable_info$variable_id[i]
  cat(i, " ")
  temp_data <-
    data.frame(value = as.numeric(expression_data[temp_variable_id,]),
               sample_info)
  
  optimize_span <-
    optimize_loess_span(
      x = temp_data$real_age,
      y = temp_data$value,
      span_range = c(0.3, 0.4, 0.5, 0.6)
    )
  
  span <-
    optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
  
  value <- temp_data$value
  real_age <- temp_data$real_age
  
  ls_reg <-
    loess(value ~ real_age,
          span = span)
  
  prediction_value =
    predict(ls_reg,
            newdata = data.frame(real_age = seq(55.5, 62, by = 0.05)))
  new_expression_data[[i]] <- as.numeric(prediction_value)
}

new_expression_data <-
  new_expression_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(new_expression_data) <-
  as.character(seq(55.5, 62, by = 0.05))

rownames(new_expression_data) <- rownames(expression_data)
save(new_expression_data, file = "new_expression_data")

load("new_expression_data")

new_expression_data <-
  new_expression_data %>% 
  dplyr::select(-c(`55.5`, `55.55`))

negative_idx <-
  new_expression_data %>%
  apply(1, function(x) {
    sum(x < 0)
  }) %>%
  `>`(0) %>%
  which()

plot(as.numeric(new_expression_data[negative_idx[1],]))

new_expression_data[which(new_expression_data < 0, arr.ind = TRUE)] <-
  0

new_sample_info <-
  data.frame(
    sample_id =
      colnames(new_expression_data),
    class = "Subject",
    subject_id = "Subject"
  )

new_variable_info <-
  variable_info

save(new_sample_info, file = "new_sample_info")
save(new_variable_info, file = "new_variable_info")

library(massdataset)

object_loess <-
  object

object_loess@expression_data <-
  new_expression_data

object_loess@sample_info <-
  new_sample_info

new_sample_info_note <-
  data.frame(
    name = c("sample_id", "class", "subject_id"),
    meaning = c("sample_id", "class", "subject_id")
  )

object_loess@sample_info_note <-
  new_sample_info_note

colnames(object_loess)

object_loess <-
  object_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!sample_id %in% c("55.5", "55.55"))

save(object_loess, file = "object_loess")

sum(is.na(object_loess))
sum(object_loess < 0)
