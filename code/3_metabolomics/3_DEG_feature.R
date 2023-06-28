no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("code/tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "data_analysis/plasma_metabolomics/data_preparation/peak/object_cross_section"
)

dir.create("data_analysis/plasma_metabolomics/DEG_feature/cross_section/",
           recursive = TRUE)

setwd("data_analysis/plasma_metabolomics/DEG_feature/cross_section")

object_cross_section

dim(object_cross_section)

###find marker which are change according to aging

###add age_range
age_index <-
  data.frame(from = c(25, 40, 45, 50, 55, 60, 65),
             to =   c(40, 45, 50, 55, 60, 65, 75.5))

apply(age_index, 1, function(x) {
  sum(
    object_cross_section@sample_info$adjusted_age > x[1] &
      object_cross_section@sample_info$adjusted_age <= x[2]
  )
})

age_range <-
  object_cross_section@sample_info$adjusted_age %>%
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

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age_range = age_range)

# ###ANOVA analysis
library(tidyverse)
library(ggpubr)
library(rstatix)

subject_data <-
  extract_expression_data(object_cross_section) %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

# anova_p <-
#   purrr::map(
#     as.data.frame(t(subject_data)),
#     .f = function(x) {
#       temp_data <-
#         data.frame(
#           subject_id = object_cross_section@sample_info$subject_id,
#           age_range = factor(
#             object_cross_section@sample_info$age_range,
#             levels = stringr::str_sort(
#               unique(object_cross_section@sample_info$age_range),
#               numeric = TRUE
#             )
#           ),
#           x = x,
#           stringsAsFactors = FALSE
#         )
#
#       res.aov <-
#         temp_data %>%
#         anova_test(x ~ age_range)
#
#       p <-
#         as.data.frame(get_anova_table(res.aov))$p
#
#       # pairwise comparisons
#       pwc <- temp_data %>%
#         pairwise_t_test(x ~ age_range,
#                         paired = FALSE,
#                         p.adjust.method = "fdr")
#       # pwc
#
#       p2 <-
#         pwc[, c(2, 3, 8, 9)]
#
#       name <- paste(pwc$group2, pwc$group1, sep = "_")
#
#       p2 <- matrix(pwc$p, nrow = 1)
#
#       colnames(p2) <- name
#
#       p <-
#         data.frame(p, p2, stringsAsFactors = FALSE, check.names = FALSE)
#
#       p
#
#     }
#   )
#
# anova_p <- anova_p %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(anova_p, file = "anova_p")
load("anova_p")

anova_p$p_adjust <-
  p.adjust(anova_p$p, method = "BH")

anova_marker_name <-
  rownames(anova_p)[which(anova_p$p_adjust < 0.05)]

# save(anova_marker_name, file = "anova_marker_name")

length(anova_marker_name)

#####TP range
subject_data2 <-
  subject_data %>%
  t() %>%
  data.frame(
    .,
    age_range = object_cross_section@sample_info$age_range,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
  mutate(
    age_range = factor(
      age_range,
      levels = object_cross_section@sample_info$age_range %>%
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
fc_p_value <-
  pbapply::pblapply(subject_data2[-1], function(x) {
    p_value <-
      lapply(1:ncol(x), function(idx) {
        t.test(x[, idx], subject_data2[[1]][, idx], paired = FALSE)$p.value
      }) %>%
      unlist() %>%
      p.adjust(method = "fdr")
    
    data.frame(
      p_value,
      variable_id = object_cross_section@variable_info$variable_id,
      stringsAsFactors = FALSE
    )
  })

# save(fc_p_value, file = "fc_p_value")

load("fc_p_value")

names(fc_p_value)

dir.create("marker_in_different_points")

##find markers for each time points
marker_each_point <-
  lapply(fc_p_value, function(x) {
    idx1 <- which(x$p_value < 0.05)
    
    gene1 <-
      try(data.frame(x[idx1,],
                     stringsAsFactors = FALSE),
          silent = TRUE)
    
    if (class(gene1) == "try-error") {
      gene1 <- NULL
    }
    gene1
  })

marker_each_point[[1]]

names(marker_each_point)

# save(marker_each_point, file = "marker_each_point")
load("marker_each_point")

#####a sankey
marker_each_point %>%
  lapply(nrow) %>%
  unlist()

all_marker_name <-
  lapply(marker_each_point, function(x) {
    x$variable_id
  }) %>%
  unlist() %>%
  unique()

length(all_marker_name)
getwd()

save(all_marker_name, file = "all_marker_name")

load("all_marker_name")

library(ggalluvial)

temp_data <-
  lapply(marker_each_point, function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x$class <- "changed"
    x <-
      data.frame(variable_id = all_marker_name,
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
  scale_fill_manual(values = c("changed" = unname(omics_color["metabolomics"]),
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
  file = file.path("marker_in_different_points", "gene_sankey_light.pdf"),
  width = 14,
  height = 7,
  bg = "transparent"
)
