no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/permutation",
  recursive = TRUE
)

setwd(
  "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/permutation"
)

object_cross_section_loess

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

####clustering
library(Mfuzz)

###permutation
# for (i in 1:1000) {
#   cat(i, " ")
#   index <-
#     sample(
#       1:ncol(object_cross_section_loess),
#       size = ncol(object_cross_section_loess),
#       replace = TRUE
#     ) %>%
#     unique() %>%
#     sort()
#
#   temp_data <-
#     object_cross_section_loess@expression_data[, index]
#
#   expression_data <-
#     object_cross_section_loess@expression_data[, index]
#
#   time <- colnames(temp_data)
#
#   temp_data <- rbind(time, temp_data)
#
#   row.names(temp_data)[1] <- "time"
#
#   write.table(
#     temp_data,
#     file = "temp_data.txt",
#     sep = '\t',
#     quote = FALSE,
#     col.names = NA
#   )
#
#   #read it back in as an expression set
#   data <- table2eset(filename = "temp_data.txt")
#   # data.s <- standardise(data)
#   data.s <- data
#   m1 <- mestimate(data.s)
#   m1
#
#   cluster_number <- 11
#
#   c <- mfuzz(data.s, c = cluster_number, m = m1)
#
#   membership_cutoff <- 0.5
#
#   cluster_info <-
#     data.frame(
#       variable_id = names(c$cluster),
#       c$membership,
#       cluster = c$cluster,
#       stringsAsFactors = FALSE
#     ) %>%
#     arrange(cluster)
#
#   save(cluster_info, file = paste("cluster_info", i, sep = "_"))
#
#   center <-
#     get_mfuzz_center(data = data.s,
#                      c = c,
#                      membership_cutoff = 0.5)
#
#   save(center, file = paste("center", i, sep = "_"))
# }

# ####calculate the p values for each molecule
# original_info <-
#   readxl::read_xlsx("../final_cluster_info.xlsx") %>%
#   as.data.frame()
# 
# load("../center")
# 
# original_center <-
#   center
# 
# info_permutation <- vector(mode = "list", length = 1000)
# center_permutation <- vector(mode = "list", length = 1000)
# for (i in 1:1000) {
#   cat(i, " ")
#   load(paste("cluster_info_", i, sep = ""))
#   info_permutation[[i]] <-
#     cluster_info %>%
#     as.data.frame()
#   
#   load(paste("center_", i, sep = ""))
#   center_permutation[[i]] <-
#     center %>%
#     as.data.frame()
# }
# 
# 
# ####cluster match
# matched_clusters <-
#   purrr::map(1:nrow(original_center), function(i) {
#     cat(i, " ")
#     temp_data <-
#       purrr::map(1:1000, function(j) {
#         idx <-
#           intersect(colnames(original_center[i, ]),
#                     colnames(center_permutation[[j]]))
#         cor_value <-
#           cor(as.numeric(original_center[i, idx]),
#               t(center_permutation[[j]][, idx]),
#               method = "spearman") %>%
#           as.numeric()
#         temp <-
#           data.frame(
#             original_cluster = i,
#             permutation_cluster = 1:nrow(center_permutation[[j]]),
#             cor = cor_value,
#             random = j
#           ) %>%
#           dplyr::filter(cor >= 0.6)
#         if (nrow(temp) == 0) {
#           temp <-
#             data.frame(
#               original_cluster = i,
#               permutation_cluster = NA,
#               cor = NA,
#               random = j
#             )
#         }
#         temp
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#     temp_data
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# #####cluster 2
# cluster2 <-
#   readxl::read_xlsx("../cluster_2/cluster2.xlsx")
# 
# library(tidyr)
# 
# matched_cluster <-
#   matched_clusters %>%
#   dplyr::filter(original_cluster == 2)
# 
# matched_cluster <-
#   split(matched_cluster, matched_cluster$random)
# 
# temp_data <-
#   purrr::map2(info_permutation,
#               matched_cluster, function(x, y) {
#                 temp <-
#                   cluster2[, c("variable_id")] %>%
#                   dplyr::left_join(x[, c("variable_id", "cluster")],
#                                    by = "variable_id") %>%
#                   as.data.frame()
#                 temp$same <-
#                   temp$cluster %in% c(y$permutation_cluster)
#                 temp[, -2]
#               })
# 
# temp_data2 <-
#   temp_data %>%
#   purrr::map(function(x) {
#     x[, 2, drop = FALSE]
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
# 
# rownames(temp_data2) <-
#   temp_data[[1]]$variable_id
# 
# p_value <-
#   cluster2$variable_id %>%
#   purrr::map(function(x) {
#     cat(x, " ")
#     
#     temp <-
#       temp_data2[which(rownames(temp_data2) == x),]
#     1 - sum(temp, na.rm = TRUE) / 1000
#   }) %>%
#   unlist()
# 
# p_value[p_value == 0] <-
#   min(p_value[p_value != 0])
# 
# plot(cluster2$membership, -log(p_value, 10))
# 
# 
# sum(p_value < 0.05)
# 
# cluster2$p_value <- p_value
# 
# openxlsx::write.xlsx(cluster2,
#                      "../cluster_2/cluster2_with_p_value.xlsx")

cluster2 <- readxl::read_xlsx("../cluster_2/cluster2_with_p_value.xlsx")

####cluster 4
cluster4 <-
  readxl::read_xlsx("../cluster_4/cluster4.xlsx")

library(tidyr)

matched_cluster <-
  matched_clusters %>%
  dplyr::filter(original_cluster == 4)

matched_cluster <-
  split(matched_cluster, matched_cluster$random)

temp_data <-
  purrr::map2(info_permutation,
              matched_cluster, function(x, y) {
                temp <-
                  cluster4[, c("variable_id")] %>%
                  dplyr::left_join(x[, c("variable_id", "cluster")],
                                   by = "variable_id") %>%
                  as.data.frame()
                temp$same <-
                  temp$cluster %in% c(y$permutation_cluster)
                temp[, -2]
              })

temp_data2 <-
  temp_data %>%
  purrr::map(function(x) {
    x[, 2, drop = FALSE]
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

rownames(temp_data2) <-
  temp_data[[1]]$variable_id

p_value <-
  cluster4$variable_id %>%
  purrr::map(function(x) {
    cat(x, " ")

    temp <-
      temp_data2[which(rownames(temp_data2) == x),]
    1 - sum(temp, na.rm = TRUE) / 1000
  }) %>%
  unlist()

p_value[p_value == 0] <-
  min(p_value[p_value != 0])

plot(cluster4$membership, -log(p_value, 10))

sum(p_value < 0.05)
sum(p_value < 0.2)
length(p_value)
cluster4$p_value <- p_value

openxlsx::write.xlsx(cluster4,
                     "../cluster_4/cluster4_with_p_value.xlsx")

cluster4 <- readxl::read_xlsx("../cluster_4/cluster4_with_p_value.xlsx")


#####cluster 5
# cluster5 <-
#   readxl::read_xlsx("../cluster_5/cluster5.xlsx")
# 
# library(tidyr)
# 
# matched_cluster <-
#   matched_clusters %>%
#   dplyr::filter(original_cluster == 5)
# 
# matched_cluster <-
#   split(matched_cluster, matched_cluster$random)
# 
# temp_data <-
#   purrr::map2(info_permutation,
#               matched_cluster, function(x, y) {
#                 temp <-
#                   cluster5[, c("variable_id")] %>%
#                   dplyr::left_join(x[, c("variable_id", "cluster")],
#                                    by = "variable_id") %>%
#                   as.data.frame()
#                 temp$same <-
#                   temp$cluster %in% c(y$permutation_cluster)
#                 temp[, -2]
#               })
# 
# temp_data2 <-
#   temp_data %>%
#   purrr::map(function(x) {
#     x[, 2, drop = FALSE]
#   }) %>%
#   do.call(cbind, .) %>%
#   as.data.frame()
# 
# rownames(temp_data2) <-
#   temp_data[[1]]$variable_id
# 
# p_value <-
#   cluster5$variable_id %>%
#   purrr::map(function(x) {
#     cat(x, " ")
#     
#     temp <-
#       temp_data2[which(rownames(temp_data2) == x),]
#     1 - sum(temp, na.rm = TRUE) / 1000
#   }) %>%
#   unlist()
# 
# p_value[p_value == 0] <-
#   min(p_value[p_value != 0])
# 
# plot(cluster5$membership, -log(p_value, 10))
# 
# sum(p_value < 0.05)
# sum(p_value < 0.2)
# length(p_value)
# cluster5$p_value <- p_value
# 
# openxlsx::write.xlsx(cluster5,
#                      "../cluster_5/cluster5_with_p_value.xlsx")


cluster5 <- readxl::read_xlsx("../cluster_5/cluster5_with_p_value.xlsx")

sum((1- cluster2$p_value > 0.8))/nrow(cluster2)
sum((1- cluster4$p_value > 0.8))/nrow(cluster4)
sum((1- cluster5$p_value > 0.8))/nrow(cluster5)
