no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "3-data_analysis/combined_omics/data_preparation/cross_section/object_cross_section"
)

dir.create("3-data_analysis/combined_omics/mci/cross_section/",
           recursive = TRUE)

setwd("3-data_analysis/combined_omics/mci/cross_section/")

object_cross_section

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

# ###calculate the correlation between age and molecules
# cor_data <-
#   seq_len(nrow(object_cross_section)) %>%
#   purrr::map(function(i){
#     cat(i, " ")
#     value <-
#       object_cross_section[i,,drop = TRUE] %>%
#       unlist()
#     age <-
#       object_cross_section@sample_info$adjusted_age
#     test <-
#       tryCatch(cor.test(age, value, method = "spearman"), error = function(e) NULL)
#     if(is.null(test)){
#       data.frame(variable_id = rownames(object_cross_section)[i],
#                  correlation = 0,
#                  p_value = 1)
#     }else{
#       data.frame(variable_id = rownames(object_cross_section)[i],
#                  correlation = unname(test$estimate),
#                  p_value = test$p.value)
#     }
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# save(cor_data, file = "cor_data")

load("cor_data")



####calculate the MCI (Maximal Information Coefficient)

# mci_data <-
#   seq_len(nrow(object_cross_section)) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     value <-
#       object_cross_section[i, , drop = TRUE] %>%
#       unlist()
#     age <-
#       object_cross_section@sample_info$adjusted_age
#
#     temp <-
#       tryCatch(mci_test(value, age), error = function(e) NULL)
#     if(is.null(temp)){
#       temp <-
#         data.frame(
#           mic = 0,
#           p_value = 1
#         )
#     }
#     cbind(variable_id = rownames(object_cross_section)[i],
#           temp)
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# save(mci_data, file = "mci_data")
load("mci_data")

temp_data <-
  cor_data %>%
  dplyr::left_join(mci_data, by = "variable_id")

library(viridis)
library(wesanderson)
library(nord)

plot <-
  temp_data %>%
  # dplyr::filter(p_value.y < 0.5) %>%
  ggplot(aes(mic, correlation)) +
  # geom_point() +
  geom_hex(bins = 100) +
  scale_fill_viridis(option = "H") +
  theme_base +
  labs(x = "Maximal information coefficient",
       y = "Spearman correlation coefficient")

ggsave(plot,
       filename = "mic vs cor.pdf",
       width = 8,
       height = 7)
