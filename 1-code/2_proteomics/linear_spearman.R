no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3-data_analysis/plasma_proteomics/data_preparation/object_cross_section")

dir.create("3-data_analysis/plasma_proteomics/linear_spearman/cross_section/",
           recursive = TRUE)

setwd("3-data_analysis/plasma_proteomics/linear_spearman/cross_section")

object_cross_section

dim(object_cross_section)

###find marker which are change according to aging

###linear regression
library(tidyverse)
library(ggpubr)
library(rstatix)

expression_data <-
  extract_expression_data(object_cross_section) %>%
  apply(1, function(x) {
    (x - mean(x))
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
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))

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

####permutation to get the p value
# dir.create("permutaton_cor_data")
# for (idx in 1:100) {
#   cat(idx, " ")
#   permutation_cor_data <-
#     seq_len(nrow(temp_object)) %>%
#     purrr::map(function(i) {
#       if((i %% 100) == 0)
#       cat(i, " ")
#       value <-
#         as.numeric(unlist(temp_object[i, , drop = TRUE]))
#       cor_result <-
#         cor.test(value, sample(temp_object@sample_info$adjusted_age,
#                                replace = FALSE), method = "spearman")
#       data.frame(
#         variable_id = rownames(temp_object)[i],
#         cor_p = cor_result$p.value,
#         spearman_cor = cor_result$estimate
#       )
#     }) %>%
#     dplyr::bind_rows() %>%
#     as.data.frame()
#   save(permutation_cor_data,
#        file = file.path("permutaton_cor_data/",
#                         paste0("permutation_cor_data_", idx)), compress = "xz")
# 
# }

load("cor_data")
load("permutaton_cor_data/permutation_cor_data_1")
permutaton_cor_data_all <-
  permutation_cor_data %>%
  dplyr::select(-cor_p)

for (idx in 2:100) {
  cat(idx, " ")
  load(file.path(
    "permutaton_cor_data/",
    paste0("permutation_cor_data_", idx)
  ))
  permutaton_cor_data_all <-
    cbind(permutaton_cor_data_all, permutation_cor_data[, 3, drop = FALSE])
}

rownames(permutaton_cor_data_all) <- NULL

permutaton_cor_data_all <-
  permutaton_cor_data_all %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "variable_id")

permutated_p_value <-
  1:nrow(cor_data) %>%
  purrr::map(function(idx) {
    cat(idx, " ")
    original_cor <-
      cor_data$spearman_cor[idx]
    permutation_cor <-
      sample(as.numeric(permutaton_cor_data_all[idx,]), 10000, replace = TRUE)
    if (original_cor > 0) {
      sum(permutation_cor > original_cor) / 10000
    } else{
      sum(permutation_cor < original_cor) / 10000
    }
    
  }) %>%
  unlist()

cor_data$permutated_p_value <-
  permutated_p_value

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p = cor_data$cor_p,
    permutated_p_value = cor_data$permutated_p_value,
    spearman_cor = cor_data$spearman_cor
  )

temp_object <-
  temp_object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(
    cor_p_adjust = p.adjust(cor_p, method = "fdr"),
    permutated_p_adjust = p.adjust(permutated_p_value, method = "fdr")
  )

plot(cor_data$cor_p,
     cor_data$permutated_p_value)

sum(temp_object@variable_info$cor_p < 0.05)
sum(temp_object@variable_info$cor_p_adjust < 0.05)
sum(temp_object@variable_info$permutated_p_value < 0.05)
sum(temp_object@variable_info$permutated_p_adjust < 0.05)

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
#     data.frame(lm_p = lm_result$p.value[2],
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
  dplyr::pull(variable_id)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) %>%
  dplyr::pull(variable_id)

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
  ggrepel::geom_text_repel(aes(
    label = ifelse(variable_id %in% top10_up_marker_name,
                   variable_id, NA)
  ), size = 3) +
  ggrepel::geom_text_repel(aes(
    label = ifelse(variable_id %in% top10_down_marker_name,
                   variable_id, NA)
  ), size = 3)

volcano_plot

# ggsave(volcano_plot,
#        filename = "volcano_plot.pdf",
#        width = 9,
#        height = 7)



top10_up_marker_name <-
  variable_info %>%
  dplyr::filter(permutated_p_value < 0.05 & spearman_cor > 0) %>%
  dplyr::arrange(desc(spearman_cor)) %>%
  head(10) %>%
  dplyr::pull(variable_id)

top10_down_marker_name <-
  variable_info %>%
  dplyr::filter(permutated_p_value < 0.05 & spearman_cor < 0) %>%
  dplyr::arrange(desc(abs(spearman_cor))) %>%
  head(10) %>%
  dplyr::pull(variable_id)

volcano_plot2 <-
  variable_info %>%
  mutate(
    marker = case_when(
      permutated_p_value < 0.05 & spearman_cor > 0 ~ "Up",
      permutated_p_value < 0.05 &
        spearman_cor < 0 ~ "Down",
      TRUE ~ "No"
    )
  ) %>%
  ggplot(aes(spearman_cor, -log(permutated_p_value, 10))) +
  geom_point(aes(size = -log(permutated_p_value, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(permuated p-values)") +
  ggrepel::geom_text_repel(aes(label = ifelse(
    variable_id %in% top10_up_marker_name,
    variable_id, NA
  )), size = 3) +
  ggrepel::geom_text_repel(aes(label = ifelse(
    variable_id %in% top10_down_marker_name,
    variable_id, NA
  )), size = 3)

volcano_plot2
# ggsave(volcano_plot2,
#        filename = "volcano_plot2.pdf",
#        width = 9,
#        height = 7)


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
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))
# scale_x_continuous(limits = c(-0.2, 0.2)) +
# scale_y_continuous(limits = c(-0.2, 0.2))

plot
ggsave(plot,
       filename = "pca_age_plot_with_markers.pdf",
       width = 9,
       height = 7)

save(pca_object, file = "pca_object")

variable_info <-
  temp_object@variable_info

save(variable_info, file = "variable_info")
