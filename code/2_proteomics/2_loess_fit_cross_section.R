no_function()

###
masstools::setwd_project()
rm(list = ls())

source("code/tools.R")
library(tidyverse)
library(tidymass)

setwd("data_analysis/plasma_proteomics/data_preparation")
load("object_cross_section")

dim(object_cross_section)

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(adjusted_age)

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
    data.frame(value = as.numeric(expression_data[temp_variable_id, ]),
               sample_info)
  
  optimize_span <-
    optimize_loess_span(
      x = temp_data$adjusted_age,
      y = temp_data$value,
      span_range = c(0.3, 0.4, 0.5, 0.6)
    )
  
  ###
  temp_data2 <-
    rbind(
      temp_data %>%
        dplyr::select(age = adjusted_age,
                      value) %>%
        dplyr::mutate(class = "real")
    )
  
  # plot <-
  #   temp_data2 %>%
  #   ggplot(aes(age, value)) +
  #   geom_point(size = 2) +
  #   geom_smooth(method = "loess",
  #               color = ggsci::pal_aaas()(n=9)[8],
  #               se = FALSE,
  #               span = 0.1) +
  #   geom_smooth(method = "loess",
  #               color = ggsci::pal_aaas()(n=9)[7],
  #               se = FALSE,
  #               span = 0.2) +
  #   geom_smooth(method = "loess",
  #               color = ggsci::pal_aaas()(n=9)[1],
  #               se = FALSE,
  #               span = 0.3) +
  #   geom_smooth(method = "loess",
  #               color = ggsci::pal_aaas()(n=9)[2],
  #               se = FALSE,
  #               span = 0.4) +
  #   geom_smooth(method = "loess",
  #               color = ggsci::pal_aaas()(n=9)[4],
  #               se = FALSE,
  #               span = 0.5) +
  #   geom_smooth(method = "loess",
  #               color = ggsci::pal_aaas()(n=9)[5],
  #               se = FALSE,
  #               span = 0.6) +
  #   geom_smooth(method = "loess",
  #               color = ggsci::pal_aaas()(n=9)[6],
  #               se = FALSE,
  #               span = 0.7) +
  #   theme_base +
  #   labs(x = "Age (years)", y = "Value", title = paste0("Span: ", span))
  #
  # ggsave(plot,
  #        filename = file.path("plot", paste0(temp_variable_id, "all_span.pdf")),
  #        width = 8, height = 7)
  
  span <-
    optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
  
  value <- temp_data$value
  adjusted_age <- temp_data$adjusted_age
  
  ls_reg <-
    loess(value ~ adjusted_age,
          span = span)
  
  prediction_value =
    predict(ls_reg,
            newdata = data.frame(adjusted_age = seq(26, 75, by = 0.5)))
  
  prediction_value2 =
    predict(ls_reg,
            newdata = data.frame(adjusted_age = temp_data$adjusted_age))
  
  # ###
  # temp_data2 <-
  #   rbind(
  #     temp_data %>%
  #       dplyr::select(age = adjusted_age,
  #                     value) %>%
  #       dplyr::mutate(class = "real"),
  #     data.frame(
  #       age = seq(26, 75, by = 0.5),
  #       value = prediction_value,
  #       class = "predicted"
  #     )
  #   )
  #
  # plot <-
  # temp_data2 %>%
  #   ggplot(aes(age, value)) +
  #   geom_point(aes(color = class), size = 2) +
  #   scale_color_manual(values = c("predicted" = "red",
  #                                 "real" = "darkblue")) +
  #   theme_base +
  #   labs(x = "Age (years)", y = "Value", title = paste0("Span: ", span))
  #
  # ggsave(plot,
  #        filename = file.path("plot", paste0(temp_variable_id, ".pdf")),
  #        width = 9, height = 7)
  
  
  new_expression_data[[i]] <- as.numeric(prediction_value)
}

new_expression_data <-
  new_expression_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(new_expression_data) <-
  as.character(seq(26, 75, by = 0.5))

rownames(new_expression_data) <- rownames(expression_data)
save(new_expression_data, file = "new_expression_data")

load("new_expression_data")

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

sum(new_expression_data < 0)
