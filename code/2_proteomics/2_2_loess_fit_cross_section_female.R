no_function()

###
setwd(masstools::get_project_wd())
rm(list = ls())

source("code/tools.R")
library(tidyverse)
library(tidymass)

setwd("data_analysis/plasma_proteomics/data_preparation/")
load("object_cross_section")

dir.create("Female")
setwd("Female/")

dim(object_cross_section)

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(adjusted_age)

###only remain male
object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(Gender == "Female")

range(object_cross_section@sample_info$adjusted_age)

##loess fit
####fore each molecule, we just use the loess method to
####impute the space between the first sample and the last sample
expression_data <-
  extract_expression_data(object_cross_section)

sample_info <-
  extract_sample_info(object_cross_section)

variable_info <-
  extract_variable_info(object_cross_section)

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
      x = temp_data$adjusted_age,
      y = temp_data$value,
      span_range = c(0.4, 0.5, 0.6)
    )
  
  span <-
    optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
  
  value <- temp_data$value
  adjusted_age <- temp_data$adjusted_age
  
  ls_reg <-
    loess(value ~ adjusted_age,
          span = span)
  
  prediction_value =
    predict(ls_reg,
            newdata = data.frame(adjusted_age = seq(26, 69, by = 0.5)))
  new_expression_data[[i]] <- as.numeric(prediction_value)
}

new_expression_data <-
  new_expression_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(new_expression_data) <-
  as.character(seq(26, 69, by = 0.5))

rownames(new_expression_data) <- rownames(expression_data)
save(new_expression_data, file = "new_expression_data")

load("new_expression_data")

sum(is.na(new_expression_data))
sum(new_expression_data < 0)
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

object_cross_section_loess <-
  object_cross_section

object_cross_section_loess@expression_data <-
  new_expression_data

object_cross_section_loess@sample_info <-
  new_sample_info

new_sample_info_note <-
  data.frame(
    name = c("sample_id", "class", "subject_id"),
    meaning = c("sample_id", "class", "subject_id")
  )

object_cross_section_loess@sample_info_note <-
  new_sample_info_note

save(object_cross_section_loess, file = "object_cross_section_loess")
