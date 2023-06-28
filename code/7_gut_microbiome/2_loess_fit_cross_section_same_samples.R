no_function()

###
setwd(masstools::get_project_wd())
rm(list = ls())

source("code/tools.R")
library(tidyverse)
library(tidymass)

dir.create("data_analysis/gut_microbiome/data_preparation/same_samples", recursive = TRUE)
setwd("data_analysis/gut_microbiome/data_preparation/same_samples")

load("../object_cross_section")

dim(object_cross_section)

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(adjusted_age)


####only remain the genus level
library(microbiomedataset)

object_cross_section <-
  object_cross_section %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(object_cross_section)

non_zero_per <-
  apply(object_cross_section, 1, function(x) {
    sum(x != 0) / ncol(object_cross_section)
  })

idx <- 
  which(non_zero_per > 0.1)

object_cross_section <-
  object_cross_section[idx,]

dim(object_cross_section)

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
      span_range = c(0.3, 0.4, 0.5, 0.6)
    )
  
  span <-
    optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
  
  value <- temp_data$value
  adjusted_age <- temp_data$adjusted_age
  
  ls_reg <-
    loess(value ~ adjusted_age,
          span = span)
  
  prediction_value2 =
    predict(ls_reg,
            newdata = data.frame(adjusted_age = temp_data$adjusted_age))
  new_expression_data[[i]] <- as.numeric(prediction_value2)
}

new_expression_data <-
  new_expression_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(new_expression_data) <-
  colnames(expression_data)

rownames(new_expression_data) <- rownames(expression_data)

negative_idx <-
  new_expression_data %>%
  apply(1, function(x) {
    sum(x < 0)
  }) %>%
  `>`(0) %>%
  which()

plot(as.numeric(new_expression_data[negative_idx[1], ]))

new_expression_data[which(new_expression_data < 0, arr.ind = TRUE)] <-
  0

library(massdataset)

object_cross_section_loess <-
  object_cross_section

object_cross_section_loess@expression_data <-
  new_expression_data

save(object_cross_section_loess, file = "object_cross_section_loess")

sum(new_expression_data < 0)
