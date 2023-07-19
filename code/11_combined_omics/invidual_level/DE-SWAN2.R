no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
setwd(masstools::get_project_wd())

source("code/tools.R")

dir.create("data_analysis/combined_omics/individual/DE_SWAN", recursive = TRUE)
setwd("data_analysis/combined_omics/individual/DE_SWAN")

load("object_loess")
object <- object_loess

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(GENETYPE)) %>%
  dplyr::filter(GENETYPE %in% c("GENETYPE", "protein-coding"))


object@sample_info$real_age <-
  as.numeric(object@sample_info$sample_id)

expression_data <-
  object@expression_data

sample_info <-
  object@sample_info

variable_info <-
  object@variable_info

range(sample_info$real_age)
plot(sort(sample_info$real_age))

temp_data <-
  do_se_swan(
    object = object,
    qt = "real_age",
    window_center = seq(56.4, 61.2, 0.2),
    buckets_size = 0.8
  )

save(temp_data, file = "temp_data")
load("temp_data")

# ##remove the background
# background <-
#   temp_data %>%
#   group_by(variable_id) %>%
#   dplyr::summarise(n = sum(p_value_adjust < 0.05))
# 
# sum(background$n > length(unique(temp_data$center)) * 0.8)
# 
# background <-
#   background %>%
#   dplyr::filter(n > length(unique(temp_data$center)) * 0.8)
# 
# temp_data <-
#   temp_data %>%
#   dplyr::filter(!variable_id %in% background$variable_id)

####all the molecules
plot <-
  temp_data %>%
  dplyr::group_by(center) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  ggplot(aes(center, number)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)")

plot

# ggsave(plot,
#        filename = "changed_molecules_total2.pdf",
#        width = 7,
#        height = 7)


####all the molecules
library(plyr)

temp_data <-
  temp_data %>%
  dplyr::left_join(variable_info[, c("variable_id", "class")],
                   by = c("variable_id"))

# temp_data$class[is.na(temp_data$class)]

plot <-
  temp_data %>%
  dplyr::group_by(class, center) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05)) %>%
  dplyr::mutate(class = factor(class, levels = names(omics_color))) %>%
  ggplot(aes(center, number)) +
  geom_point(aes(color = class),
             show.legend = FALSE) +
  geom_line(aes(group = 1,
                color = class),
            show.legend = FALSE) +
  scale_color_manual(values = omics_color) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)") +
  facet_wrap(facets = vars(class),
             scales = "free_y",
             nrow = 2)

plot

# ggsave(plot,
#        filename = "changed_molecules2.pdf",
#        width = 14,
#        height = 6)

####heatmap
library(ComplexHeatmap)
temp <-
  temp_data %>%
  dplyr::group_by(class, center) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05)) %>%
  dplyr::mutate(class = factor(class, levels = names(omics_color))) %>%
  tidyr::pivot_wider(names_from = "center", values_from = "number") %>%
  tibble::column_to_rownames(var = "class") %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t()

library(circlize)
col_fun = colorRamp2(seq(-2, 2, length.out = 11),
                     rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))

plot <-
  Heatmap(
    temp,
    col = col_fun,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = "Number (Z-score)",
                                direction = "horizontal"),
    border = TRUE,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    column_names_rot = 0
  )

plot <- ggplotify::as.ggplot(plot)
plot

# ggsave(plot, filename = "heatmap_for_changed_molecules.pdf",
#        width = 10, height = 4)


temp_data <-
  data.frame(
    class = names(omics_color),
    crest1 = c(44, 44, 47, 51, 45, 45, 43, 43, 41, 46),
    crest2 = c(61, 57, 57, 61, NA, NA, 61, 58, NA, 58)
  )

plot <-
  temp_data %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "crest",
                      values_to = "age") %>%
  # dplyr::filter(crest == "crest1") %>%
  # dplyr::mutate(age = as.character(age)) %>%
  dplyr::filter(!is.na(age)) %>%
  # dplyr::mutate(age = factor(age, levels = as.character(40:65))) %>%
  ggplot(aes(x = age)) +
  geom_bar(
    stat = "count",
    aes(fill = class),
    color = "black",
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = 41:65) +
  theme_base +
  scale_fill_manual(values = omics_color)
plot
# ggsave(plot, filename = "changed_molecules_distributation.pdf",
#        width = 10, height = 2)

######use different cutoff
load("temp_data")
plot <-
  temp_data %>%
  dplyr::group_by(center) %>%
  dplyr::summarise(
    "< 0.05" = sum(p_value_adjust < 0.05),
    "< 0.01" = sum(p_value_adjust < 0.01),
    "< 0.001" = sum(p_value_adjust < 0.001),
    "< 0.0001" = sum(p_value_adjust < 0.0001)
  ) %>%
  tidyr::pivot_longer(cols = -center,
                      names_to = "cutoff",
                      values_to = "number") %>%
  ggplot(aes(center, number)) +
  geom_point(aes(color = cutoff)) +
  geom_line(aes(group = cutoff, color = cutoff)) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)") +
  theme(legend.position = c(0, 0),
        legend.justification = c(0, 0)) +
  ggsci::scale_color_lancet()

plot

# ggsave(plot,
#        filename = "changed_molecules_total2_different_cutoff.pdf",
#        width = 7,
#        height = 7)

######use different cutoff
plot <-
  temp_data %>%
  dplyr::group_by(class, center) %>%
  dplyr::summarise(
    "< 0.05" = sum(p_value_adjust < 0.05),
    "< 0.01" = sum(p_value_adjust < 0.01),
    "< 0.001" = sum(p_value_adjust < 0.001),
    "< 0.0001" = sum(p_value_adjust < 0.0001)
  ) %>%
  tidyr::pivot_longer(
    cols = -c(class, center),
    names_to = "cutoff",
    values_to = "number"
  ) %>%
  dplyr::mutate(class = factor(class, levels = names(omics_color))) %>%
  ggplot(aes(center, number)) +
  # geom_point(aes(color = cutoff)) +
  geom_line(aes(group = cutoff, color = cutoff)) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)") +
  theme(legend.position = "bottom") +
  ggsci::scale_color_lancet() +
  facet_wrap(facets = vars(class),
             scales = "free_y",
             nrow = 2)

plot

# ggsave(plot,
#        filename = "changed_molecules2_different_cutoff.pdf",
#        width = 14,
#        height = 6)


#######use different windows
# temp_data_window15 <-
#   do_se_swan(
#     object = object,
#     qt = "age",
#     window_center = seq(40, 65, 1),
#     buckets_size = 15
#   )
#
# save(temp_data_window15, file = "temp_data_window15")
#
# temp_data_window25 <-
#   do_se_swan(
#     object = object,
#     qt = "age",
#     window_center = seq(40, 65, 1),
#     buckets_size = 25
#   )
#
# save(temp_data_window25, file = "temp_data_window25")
#
# temp_data_window30 <-
#   do_se_swan(
#     object = object,
#     qt = "age",
#     window_center = seq(40, 65, 1),
#     buckets_size = 30
#   )
#
# save(temp_data_window30, file = "temp_data_window30")


load("temp_data_window15")
load("temp_data")
load("temp_data_window25")
load("temp_data_window30")

temp15 <-
  temp_data_window15 %>%
  dplyr::group_by(center) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 15)

temp20 <-
  temp_data %>%
  dplyr::group_by(center) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 20)

temp25 <-
  temp_data_window25 %>%
  dplyr::group_by(center) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 25)

temp30 <-
  temp_data_window30 %>%
  dplyr::group_by(center) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 30)

temp <-
  rbind(temp15,
        temp20,
        temp25,
        temp30)

######use different cutoff
plot <-
  temp %>%
  dplyr::mutate(window = as.character(window)) %>%
  ggplot(aes(center, number)) +
  geom_point(aes(color = window)) +
  geom_line(aes(group = window, color = window)) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)") +
  theme(legend.position = c(0, 0),
        legend.justification = c(0, 0)) +
  ggsci::scale_color_nejm()

plot

# ggsave(plot,
#        filename = "changed_molecules_total2_different_window.pdf",
#        width = 7,
#        height = 7)


temp15 <-
  temp_data_window15 %>%
  dplyr::left_join(variable_info[, c("variable_id", "class")],
                   by = "variable_id") %>%
  dplyr::group_by(center, class) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 15)

temp20 <-
  temp_data %>%
  dplyr::left_join(variable_info[, c("variable_id", "class")],
                   by = "variable_id") %>%
  dplyr::group_by(center, class) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 20)

temp25 <-
  temp_data_window25 %>%
  dplyr::left_join(variable_info[, c("variable_id", "class")],
                   by = "variable_id") %>%
  dplyr::group_by(center, class) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 25)

temp30 <-
  temp_data_window30 %>%
  dplyr::left_join(variable_info[, c("variable_id", "class")],
                   by = "variable_id") %>%
  dplyr::group_by(center, class) %>%
  dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE)) %>%
  dplyr::mutate(window = 30)

temp <-
  rbind(temp15,
        temp20,
        temp25,
        temp30)

######use different cutoff
plot <-
  temp %>%
  dplyr::mutate(window = as.character(window)) %>%
  dplyr::mutate(class = factor(class, levels = names(omics_color))) %>%
  ggplot(aes(center, number)) +
  # geom_point(aes(color = window)) +
  geom_line(aes(group = window,
                color = window)) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)") +
  theme(legend.position = "bottom") +
  ggsci::scale_color_nejm() +
  facet_wrap(facets = vars(class),
             scales = "free_y",
             nrow = 2)

plot

# ggsave(plot,
#        filename = "changed_molecules2_different_window.pdf",
#        width = 14,
#        height = 6)


###Permutation test

# temp_object <-
#   object
#
# temp_object@sample_info$age <-
#   sample(temp_object@sample_info$age)
#
# result <-
#   do_se_swan(
#     object = temp_object,
#     qt = "age",
#     window_center = seq(40, 65, 1),
#     buckets_size = 20
#   )
#
# temp <-
#   result %>%
#   dplyr::group_by(center) %>%
#   dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE))
#
# plot <-
#   temp %>%
#   ggplot(aes(center, number)) +
#   geom_point() +
#   geom_line(aes(group = 1)) +
#   base_theme +
#   labs(x = "Age (years)",
#        y = "# significant molecules (q < 0.05)")
# plot
#
# ggsave(plot,
#        filename = "changed_molecules_total2_permutation.pdf",
#        width = 7,
#        height = 7)
#
# temp <-
#   result %>%
#   dplyr::left_join(variable_info[, c("variable_id", "class")],
#                    by = "variable_id") %>%
#   dplyr::group_by(center, class) %>%
#   dplyr::summarise(number = sum(p_value_adjust < 0.05, na.rm = TRUE))
#
# plot <-
#   temp %>%
#   dplyr::mutate(class = factor(class, levels = names(omics_color))) %>%
#   ggplot(aes(center, number)) +
#   geom_point(aes(color = class)) +
#   geom_line(aes(group = class, color = class)) +
#   base_theme +
#   labs(x = "Age (years)",
#        y = "# significant molecules (q < 0.05)") +
#   theme(legend.position = "bottom") +
#   scale_color_manual(values = omics_color) +
#   facet_wrap(facets = vars(class),
#              scales = "free_y",
#              nrow = 2)
# plot
# ggsave(plot,
#        filename = "changed_molecules_permutation.pdf",
#        width = 14,
#        height = 6)


###overlap between two different waves
crest_info <-
  data.frame(
    class = c(
      "all",
      "transcriptome",
      "proteomics",
      "metabolomics",
      "cytokine",
      "gut_microbiome",
      "skin_microbiome",
      "nasal_microbiome"
    ),
    crest1 = c(44, 44, 44, 47, 51, 43, 43, 46),
    crest2 = c(60, 60, 57, 57, 61, 61, 58, 58)
  )

library(ggVennDiagram)

# seq_len(nrow(crest_info)) %>%
#   purrr::map(function(i) {
#     if (crest_info$class[i] == "all") {
#       id1 <-
#         temp_data %>%
#         dplyr::filter(center == 44) %>%
#         dplyr::filter(p_value_adjust < 0.05) %>%
#         pull(variable_id)
#
#       id2 <-
#         temp_data %>%
#         dplyr::filter(center == 60) %>%
#         dplyr::filter(p_value_adjust < 0.05) %>%
#         pull(variable_id)
#     } else{
#       id1 <-
#         temp_data %>%
#         dplyr::left_join(variable_info[, c("variable_id", "class")],
#                          by = "variable_id") %>%
#         dplyr::filter(class == crest_info$class[i]) %>%
#         dplyr::filter(center == crest_info$crest1[i]) %>%
#         dplyr::filter(p_value_adjust < 0.05) %>%
#         pull(variable_id)
#
#       id2 <-
#         temp_data %>%
#         dplyr::left_join(variable_info[, c("variable_id", "class")],
#                          by = "variable_id") %>%
#         dplyr::filter(class == crest_info$class[i]) %>%
#         dplyr::filter(center == crest_info$crest2[i]) %>%
#         dplyr::filter(p_value_adjust < 0.05) %>%
#         pull(variable_id)
#     }
#
#     x <- list(crest1 = id1,
#               crest2 = id2)
#     plot <-
#       ggVennDiagram(x = x) +
#       scale_fill_gradient(low = "white",
#                           high = "red")
#
#     plot
#     ggsave(
#       plot,
#       filename = paste0("crest_overlap_", crest_info$class[i], ".pdf"),
#       width = 7,
#       height = 7
#     )
#   })
