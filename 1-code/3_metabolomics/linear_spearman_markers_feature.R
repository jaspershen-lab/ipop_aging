no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "3-data_analysis/plasma_metabolomics/data_preparation/peak/object_cross_section"
)

dir.create("3-data_analysis/plasma_metabolomics/linear_spearman_feature/cross_section/",
           recursive = TRUE)

setwd("3-data_analysis/plasma_metabolomics/linear_spearman_feature/cross_section")

object_cross_section

dim(object_cross_section)

###find marker which are change according to aging

###linear regression
library(tidyverse)
library(ggpubr)
library(rstatix)

expression_data <-
  extract_expression_data(object_cross_section) %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

sample_info <-
  object_cross_section@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
expression_data <-
  lm_adjust(expression_data = expression_data,
            sample_info = sample_info,
            threads = 3)

temp_object <- object_cross_section
temp_object@expression_data <- expression_data

# ##correlation
# ##cor
# cor_data <-
#   seq_len(nrow(temp_object)) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     value <-
#       as.numeric(unlist(temp_object[i, , drop = TRUE]))
#     cor_result <-
#       cor.test(value, temp_object@sample_info$adjusted_age, method = "spearman")
#     data.frame(variable_id = rownames(temp_object)[i],
#       cor_p = cor_result$p.value,
#                spearman_cor = cor_result$estimate)
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
# save(cor_data, file = "cor_data")

load("cor_data")

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p = cor_data$cor_p,
                spearman_cor = cor_data$spearman_cor)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p_adjust = p.adjust(cor_p, method = "fdr"))

variable_info <-
  extract_variable_info(temp_object)

volcano_plot <-
  variable_info %>%
  mutate(
    marker = case_when(
      cor_p_adjust < 0.05 & spearman_cor > 0 ~ "Up",
      cor_p_adjust < 0.05 &
        spearman_cor < 0 ~ "Down",
      TRUE ~ "No"
    )
  ) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)")

volcano_plot

ggsave(volcano_plot,
       filename = "volcano_plot.png",
       width = 9,
       height = 7)


dim(variable_info)
sum(variable_info$cor_p < 0.05)

x <-
variable_info %>% 
  dplyr::filter(cor_p < 0.05) %>% 
  pull(Compound.name)

length(x[!is.na(x)])
