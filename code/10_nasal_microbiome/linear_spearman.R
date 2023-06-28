no_source()

rm(list = ls())
masstools::setwd_project()
source("code/tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("data_analysis/nasal_microbiome/data_preparation/object_cross_section")

dir.create("data_analysis/nasal_microbiome/linear_spearman/cross_section/",
           recursive = TRUE)

setwd("data_analysis/nasal_microbiome/linear_spearman/cross_section")

object_cross_section

dim(object_cross_section)

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
  object_cross_section[idx, ]

object_cross_section <-
  object_cross_section %>%
  transform2relative_intensity()

dim(object_cross_section)

###find marker which are change according to aging

###linear regression
library(tidyverse)
library(ggpubr)
library(rstatix)

expression_data <-
  extract_expression_data(object_cross_section) %>%
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

####PCA analysis
pca_object <-
  temp_object %>%
  massstat::run_pca()

plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "adjusted_age",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Reds"))

plot

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
  mutate(cor_p = cor_data$cor_p,
         spearman_cor = cor_data$spearman_cor)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p_adjust = p.adjust(cor_p, method = "fdr"))

######linear mixed model
library(lme4)
# 
# linear_mix_model_data <-
#   seq_len(nrow(temp_object)) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     value <-
#       as.numeric(unlist(temp_object[i, , drop = TRUE]))
#     sample_info <-
#       temp_object@sample_info
#     temp_data <-
#       data.frame(sample_info, value)
# 
# 
#     lm_result <-
#       glm(formula = value ~ adjusted_age, data = temp_data)
# 
#     lm_result <-
#       lm_result %>%
#       broom::tidy()
# 
#     data.frame(cor_p_adjust = lm_result$p.value[2],
#                coefficient = lm_result$estimate[2])
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
# save(linear_mix_model_data, file = "linear_mix_model_data")
load("linear_mix_model_data")

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  mutate(lm_p = linear_mix_model_data$lm_p,
         coefficient = linear_mix_model_data$coefficient)

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(lm_p_adjust = p.adjust(lm_p, method = "fdr"))

plot(temp_object@variable_info$spearman_cor,
     temp_object@variable_info$coefficient)

temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point() +
  theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta")

###volcano plot
variable_info <-
  extract_variable_info(temp_object)

top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) %>%
  dplyr::pull(Genus)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) %>%
  dplyr::pull(Genus)

volcano_plot <-
  variable_info %>%
  mutate(marker = case_when(
    cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
    cor_p_adjust < 0.2 &
      spearman_cor < 0 ~ "Down",
    TRUE ~ "No"
  )) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)") +
  ggrepel::geom_text_repel(aes(label = ifelse(
    Genus %in% top10_up_marker_name,
    Genus, NA
  )), size = 3) +
  ggrepel::geom_text_repel(aes(label = ifelse(
    Genus %in% top10_down_marker_name,
    Genus, NA
  )), size = 3)

volcano_plot

ggsave(volcano_plot,
       filename = "volcano_plot.pdf",
       width = 9,
       height = 7)

####differential expressional metabolites
dim(variable_info)
aging_markers <-
  variable_info %>%
  dplyr::filter(cor_p < 0.05)

pca_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% aging_markers$variable_id) %>%
  massstat::run_pca()

plot <-
  massstat::pca_score_plot(
    object = temp_object,
    pca_object = pca_object,
    color_by = "adjusted_age",
    frame = FALSE
  ) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  scale_y_continuous(limits = c(-0.2, 0.2))

plot
ggsave(plot,
       filename = "pca_age_plot_with_markers.pdf",
       width = 9,
       height = 7)

variable_info <-
  temp_object@variable_info
save(pca_object, file = "pca_object")
save(variable_info, file = "variable_info")
