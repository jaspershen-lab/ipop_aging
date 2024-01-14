no_source()
rm(list = ls())

setwd(r4projects::get_project_wd())

source("1-code/100-tools.R")

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("3-data_analysis/combined_omics/DE_SWAN")
setwd("3-data_analysis/combined_omics/DE_SWAN")

library("DEswan")

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id))

expression_data <-
  object_cross_section_loess@expression_data

sample_info <-
  object_cross_section_loess@sample_info

variable_info <-
  object_cross_section_loess@variable_info

temp_data <-
  data.frame(sample_info[, c("age"), drop = FALSE],
             t(expression_data))


head(temp_data[, 1:5])

range(temp_data$age)
range(agingplasmaproteome$Age)

# res.DEswan = DEswan(
#   data.df = temp_data[, -1],
#   qt = temp_data[, 1],
#   window.center = seq(40, 65, 1),
#   buckets.size = 20
# )
#
# save(res.DEswan, file = "res.DEswan")
load("res.DEswan")

head(res.DEswan$p)
head(res.DEswan$coeff)

res.DEswan.wide.p <-
  reshape.DEswan(res.DEswan, parameter = 1, factor = "qt")
head(res.DEswan.wide.p[, 1:5])

res.DEswan.wide.q = q.DEswan(res.DEswan.wide.p, method = "BH")
head(res.DEswan.wide.q[, 1:5])


####all the molecules
res.DEswan.wide.q %>%
  dplyr::select(-variable) %>%
  purrr::map(function(x) {
    sum(x < 0.05)
  }) %>%
  unlist() %>%
  data.frame(number = .) %>%
  tibble::rownames_to_column(var = "age") %>%
  dplyr::mutate(age = as.numeric(stringr::str_replace(age, "X", ""))) %>%
  ggplot(aes(age, number)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)")

####all the molecules
library(plyr)
temp_data <-
  res.DEswan.wide.q %>%
  dplyr::left_join(variable_info[, c("variable_id", "class")],
                   by = c("variable" = "variable_id"))

temp_data$class[is.na(temp_data$class)] <-
  "transcriptome"

plot <-
  temp_data %>%
  dplyr::select(-variable) %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(data) {
    class <- data$class[1]
    data <-
      data %>%
      dplyr::select(-class) %>%
      purrr::map(function(x) {
        sum(x < 0.05)
      }) %>%
      unlist() %>%
      data.frame(number = .) %>%
      tibble::rownames_to_column(var = "age") %>%
      dplyr::mutate(age = as.numeric(stringr::str_replace(age, "X", "")))
    data$class <- class
    data
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(class = factor(class, levels = names(omics_color))) %>%
  ggplot(aes(age, number)) +
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
             nrow = 3)

plot

# ggsave(plot,
#        filename = "changed_molecules.pdf",
#        width = 14,
#        height = 10)



plot <-
  temp_data %>%
  dplyr::select(-c(variable, class)) %>%
  purrr::map(function(x) {
    sum(x < 0.05)
  }) %>%
  unlist() %>%
  data.frame(number = .) %>%
  tibble::rownames_to_column(var = "age") %>%
  dplyr::mutate(age = as.numeric(stringr::str_replace(age, "X", ""))) %>%
  ggplot(aes(age, number)) +
  # geom_point() +
  geom_line(aes(group = 1)) +
  base_theme +
  labs(x = "Age (years)",
       y = "# significant molecules (q < 0.05)")

plot

# ggsave(plot,
#        filename = "changed_molecules_total.pdf",
#        width = 7,
#        height = 7)
