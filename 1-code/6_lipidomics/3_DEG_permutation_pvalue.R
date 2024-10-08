no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3-data_analysis/plasma_lipidomics/data_preparation/object_cross_section_loess")

dir.create("3-data_analysis/plasma_lipidomics/DEG/cross_section_loess/",
           recursive = TRUE)

setwd("3-data_analysis/plasma_lipidomics/DEG/cross_section_loess")

object_cross_section_loess

dim(object_cross_section_loess)

###find marker which are change according to aging

###add age_range
age_index <-
  data.frame(from = c(25, 40, 45, 50, 55, 60, 65),
             to =   c(40, 45, 50, 55, 60, 65, 75.5))

apply(age_index, 1, function(x) {
  sum(
    object_cross_section_loess@sample_info$adjusted_age > x[1] &
      object_cross_section_loess@sample_info$adjusted_age <= x[2]
  )
})

age_range <-
  object_cross_section_loess@sample_info$sample_id %>%
  as.numeric() %>% 
  purrr::map(function(x) {
    idx <-
      apply(age_index, 1, function(y) {
        x > y[1] &
          x <= y[2]
      }) %>%
      which()
    paste(age_index$from, age_index$to, sep = "_")[idx]
  }) %>%
  unlist()

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age_range = age_range)

subject_data <-
  extract_expression_data(object_cross_section_loess) %>%
  apply(1, function(x) {
    (x - mean(x))
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

#####TP range
subject_data2 <-
  subject_data %>%
  t() %>%
  data.frame(
    .,
    age_range = object_cross_section_loess@sample_info$age_range,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
  mutate(
    age_range = factor(
      age_range,
      levels = object_cross_section_loess@sample_info$age_range %>%
        unique() %>%
        stringr::str_sort(numeric = TRUE)
    )
  ) %>%
  plyr::dlply(.variables = .(age_range))

subject_data2 <-
  lapply(subject_data2, function(x) {
    x <-
      x %>%
      dplyr::select(-age_range)
  })

#find all the peaks in different time points
original_difference <-
  2:length(subject_data2) %>%
  purrr::map(function(idx) {
    x <-
      subject_data2[[idx]]
    difference <-
      1:ncol(x) %>%
      purrr::map(function(i) {
        abs(mean(subject_data2[[1]][, i]) -
              mean(x[, i]))
      }) %>%
      unlist()
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

rownames(original_difference) <-
  colnames(subject_data2[[1]])

colnames(original_difference) <-
  names(subject_data2)[-1]

save(original_difference, file = "original_difference")

####permutation
dir.create("permutation")
# for (idx in 2:length(subject_data2)) {
#   cat(idx)
#   x <-
#     subject_data2[[idx]]
#   permutated_difference <-
#     purrr::map(1:100, function(j) {
#       cat(j, " ")
#       difference <-
#         1:ncol(x) %>%
#         purrr::map(function(i) {
#           value1 <-
#             subject_data[i, sample(1:ncol(subject_data),
#                                    length(subject_data2[[1]][, i]),
#                                    replace = FALSE,), drop = TRUE] %>%
#             unlist()
# 
#           value2 <-
#             subject_data[i, sample(1:ncol(subject_data),
#                                    length(subject_data2[[1]][, i]),
#                                    replace = FALSE,), drop = TRUE] %>%
#             unlist()
# 
#           abs(mean(value1) -
#                 mean(value2))
#         }) %>%
#         unlist()
#       difference
#     }) %>%
#     do.call(cbind, .) %>%
#     as.data.frame()
# 
#   save(permutated_difference,
#        file = file.path("permutation/",
#                         paste0("permutated_difference_", idx)))
# 
# }

###get the p value
load("original_difference")

# fc_p_value_permutation <-
# 2:length(subject_data2) %>%
#   purrr::map(function(idx) {
#     load(file.path("permutation/",
#                    paste0("permutated_difference_", idx)))
#     p_value <-
#       purrr::map(1:nrow(original_difference),
#                  function(i) {
#                    original_value <-
#                      original_difference[i, idx - 1]
#                    permutation_value <-
#                      permutated_difference[i,] %>%
#                      as.numeric()
#                    permutation_value <-
#                      sample(permutation_value, 10000, replace = TRUE)
#                    sum(permutation_value > original_value) / 10000
#                  }) %>%
#       unlist()
#     data.frame(p_value = p_value,
#                variable_id = rownames(original_difference))
#   })
# 
# names(fc_p_value_permutation) <-
#   names(original_difference)
# 
# save(fc_p_value_permutation, file = "fc_p_value_permutation")

load("fc_p_value_permutation")

names(fc_p_value_permutation)

dir.create("marker_in_different_points")

##find markers for each time points
marker_each_point_permutation <-
  lapply(fc_p_value_permutation, function(x) {
    idx1 <- which(p.adjust(x$p_value, method = "fdr") < 0.05)
    
    gene1 <-
      try(data.frame(x[idx1, ],
                     stringsAsFactors = FALSE),
          silent = TRUE)
    
    if (class(gene1) == "try-error") {
      gene1 <- NULL
    }
    gene1
  })

marker_each_point_permutation[[1]]

names(marker_each_point_permutation) <-
  names(fc_p_value_permutation)

save(marker_each_point_permutation, file = "marker_each_point_permutation")
load("marker_each_point_permutation")

#####a sankey
marker_each_point_permutation %>%
  lapply(nrow) %>%
  unlist()

all_marker_name_permutation <-
  lapply(marker_each_point_permutation, function(x) {
    x %>% 
      dplyr::filter(p_value * 6 < 0.05) %>% 
      pull(variable_id)
  }) %>%
  unlist() %>%
  unique()

length(all_marker_name_permutation)
getwd()

save(all_marker_name_permutation, file = "all_marker_name_permutation")

load("all_marker_name_permutation")

library(ggalluvial)

temp_data <-
  lapply(marker_each_point_permutation, function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x$class <- "changed"
    x <-
      data.frame(variable_id = all_marker_name_permutation,
                 stringsAsFactors = FALSE) %>%
      left_join(x, by = "variable_id") %>%
      dplyr::select(variable_id, class)
    
    x$class[is.na(x$class)] <- "no"
    x$freq <- 1
    x
    
  })


temp_data <-
  purrr::map2(
    .x = temp_data,
    .y = names(temp_data),
    .f = function(x, y) {
      if (is.null(x)) {
        return(NULL)
      }
      data.frame(x, point = y, stringsAsFactors = FALSE)
    }
  )

temp_data <-
  do.call(rbind, temp_data)

temp_data$point <-
  factor(temp_data$point, levels = unique(temp_data$point))

RColorBrewer::display.brewer.all()

plot1 <-
  ggplot(
    temp_data,
    aes(
      x = point,
      y = freq,
      stratum = class,
      alluvium = variable_id,
      fill = class,
      label = class
    )
  ) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_flow() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("changed" = unname(omics_color["lipidomics"]),
                               "no" = "grey")) +
  ggalluvial::geom_stratum(alpha = 1, color = "black") +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 2),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot1

ggsave(
  plot1,
  file = file.path(
    "marker_in_different_points",
    "gene_sankey_light_permutation.pdf"
  ),
  width = 14,
  height = 7,
  bg = "transparent"
)
