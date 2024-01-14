no_function()

setwd(r4projects::get_project_wd())
rm(list = ls())

source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)
library(Mfuzz)

setwd(r4projects::get_project_wd())

setwd("3-data_analysis/combined_omics/")

####female
load("Female/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info")
final_cluster_info_female <- final_cluster_info

load("Female/fuzzy_c_means_clustering/cross_section_loess/c")
c_female <- c

data <-
  table2eset(filename = "Female/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt")
data.s <- data
female_center <-
  get_mfuzz_center(data = data.s,
                   c = c_female,
                   membership_cutoff = 0.5)
rownames(female_center) <-
  paste("Cluster", rownames(female_center), sep = ' ')

####male
load("Male/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info")
final_cluster_info_male <- final_cluster_info

load("Male/fuzzy_c_means_clustering/cross_section_loess/c")
c_male <- c

data <-
  table2eset(filename = "Male/fuzzy_c_means_clustering/cross_section_loess/temp_data.txt")
data.s <- data
male_center <-
  get_mfuzz_center(data = data.s,
                   c = c_male,
                   membership_cutoff = 0.5)
rownames(male_center) <-
  paste("Cluster", rownames(male_center), sep = ' ')

####all
load("fuzzy_c_means_clustering/cross_section_loess/final_cluster_info")
final_cluster_info_all <- final_cluster_info

load("fuzzy_c_means_clustering/cross_section_loess/c")
c_all <- c

data <-
  table2eset(filename = "fuzzy_c_means_clustering/cross_section_loess/temp_data.txt")
data.s <- data
all_center <-
  get_mfuzz_center(data = data.s,
                   c = c_all,
                   membership_cutoff = 0.5)
rownames(all_center) <-
  paste("Cluster", rownames(all_center), sep = ' ')

dir.create("summary_male_female_clustering",
           recursive = TRUE,
           showWarnings = FALSE)

setwd("summary_male_female_clustering/")

female_center
male_center
all_center

colnames(female_center)
colnames(male_center)

######calculate the correlation between clusters
#####male with all
# intersect_sample_id <-
#   intersect(colnames(male_center),
#             colnames(all_center))
#
# temp_data1 <-
#   male_center[, intersect_sample_id]
#
# temp_data2 <-
#   all_center[, intersect_sample_id]
#
#
# cor_data_male_all <-
#   seq_len(nrow(temp_data1)) %>%
#   purrr::map(function(i) {
#     seq_len(nrow(temp_data2)) %>%
#       purrr::map(function(j) {
#         temp <-
#           cor.test(as.numeric(temp_data1[i,]),
#                    as.numeric(temp_data2[j,]))
#         data.frame(
#           from_cluster = i,
#           to_cluster = j,
#           cor = temp$estimate,
#           p = temp$p.value
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame
#
# rownames(cor_data_male_all) <- NULL
#
# cor_data_male_all$from_data <- "male"
# cor_data_male_all$to_data <- "all"
#
#
#
# #####female with all
# intersect_sample_id <-
#   intersect(colnames(female_center),
#             colnames(all_center))
#
# temp_data1 <-
#   female_center[, intersect_sample_id]
#
# temp_data2 <-
#   all_center[, intersect_sample_id]
#
# cor_data_female_all <-
#   seq_len(nrow(temp_data1)) %>%
#   purrr::map(function(i) {
#     seq_len(nrow(temp_data2)) %>%
#       purrr::map(function(j) {
#         temp <-
#           cor.test(as.numeric(temp_data1[i,]),
#                    as.numeric(temp_data2[j,]))
#         data.frame(
#           from_cluster = i,
#           to_cluster = j,
#           cor = temp$estimate,
#           p = temp$p.value
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame
#
# rownames(cor_data_female_all) <- NULL
#
# cor_data_female_all$from_data <- "female"
# cor_data_female_all$to_data <- "all"
#
#
# #####female with male
# intersect_sample_id <-
#   intersect(colnames(female_center),
#             colnames(male_center))
#
# temp_data1 <-
#   female_center[, intersect_sample_id]
#
# temp_data2 <-
#   male_center[, intersect_sample_id]
#
# cor_data_female_male <-
#   seq_len(nrow(temp_data1)) %>%
#   purrr::map(function(i) {
#     seq_len(nrow(temp_data2)) %>%
#       purrr::map(function(j) {
#         temp <-
#           cor.test(as.numeric(temp_data1[i,]),
#                    as.numeric(temp_data2[j,]))
#         data.frame(
#           from_cluster = i,
#           to_cluster = j,
#           cor = temp$estimate,
#           p = temp$p.value
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame
#
# rownames(cor_data_female_male) <- NULL
#
# cor_data_female_male$from_data <- "female"
# cor_data_female_male$to_data <- "male"
#
# cor_data <-
#   rbind(
#     cor_data_female_all,
#     cor_data_male_all,
#     cor_data_female_male
#   )
# save(cor_data, file = "cor_data")

load("cor_data")

#####calculate the jaccard index
####male with all
# final_cluster_info_male
# final_cluster_info_all
#
# jaccard_index_male_all <-
# unique(final_cluster_info_male$cluster) %>%
#   purrr::map(function(i) {
#     unique(final_cluster_info_all$cluster) %>%
#       purrr::map(function(j) {
#         temp1 <-
#           final_cluster_info_male %>%
#           dplyr::filter(cluster == as.character(i))
#
#         temp2 <-
#           final_cluster_info_all %>%
#           dplyr::filter(cluster == as.character(j))
#
#         jaccard_index <-
#           length(intersect(temp1$variable_id, temp2$variable_id)) / length(unique(c(
#             temp1$variable_id, temp2$variable_id
#           )))
#         data.frame(from_cluster = i,
#                    to_cluster = j,
#                    jaccard_index = jaccard_index)
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# rownames(jaccard_index_male_all) <- NULL
#
# jaccard_index_male_all$from_data <- "male"
# jaccard_index_male_all$to_data <- "all"
#
#
#
#
#
# ####female with all
# jaccard_index_female_all <-
#   unique(final_cluster_info_female$cluster) %>%
#   purrr::map(function(i) {
#     unique(final_cluster_info_all$cluster) %>%
#       purrr::map(function(j) {
#         temp1 <-
#           final_cluster_info_female %>%
#           dplyr::filter(cluster == as.character(i))
#
#         temp2 <-
#           final_cluster_info_all %>%
#           dplyr::filter(cluster == as.character(j))
#
#         jaccard_index <-
#           length(intersect(temp1$variable_id, temp2$variable_id)) / length(unique(c(
#             temp1$variable_id, temp2$variable_id
#           )))
#         data.frame(from_cluster = i,
#                    to_cluster = j,
#                    jaccard_index = jaccard_index)
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# rownames(jaccard_index_female_all) <- NULL
#
# jaccard_index_female_all$from_data <- "female"
# jaccard_index_female_all$to_data <- "all"
#
#
# ####female with male
# jaccard_index_female_male <-
#   unique(final_cluster_info_female$cluster) %>%
#   purrr::map(function(i) {
#     unique(final_cluster_info_male$cluster) %>%
#       purrr::map(function(j) {
#         temp1 <-
#           final_cluster_info_female %>%
#           dplyr::filter(cluster == as.character(i))
#
#         temp2 <-
#           final_cluster_info_male %>%
#           dplyr::filter(cluster == as.character(j))
#
#         jaccard_index <-
#           length(intersect(temp1$variable_id, temp2$variable_id)) / length(unique(c(
#             temp1$variable_id, temp2$variable_id
#           )))
#         data.frame(from_cluster = i,
#                    to_cluster = j,
#                    jaccard_index = jaccard_index)
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# rownames(jaccard_index_female_male) <- NULL
#
# jaccard_index_female_male$from_data <- "female"
# jaccard_index_female_male$to_data <- "male"
#
#
# jaccard_index <-
#   rbind(
#     jaccard_index_female_all,
#     jaccard_index_male_all,
#     jaccard_index_female_male
#   )
# save(jaccard_index, file = "jaccard_index")

load("jaccard_index")

head(cor_data)
head(jaccard_index)

###find the cluster 2
cor_data %>%
  dplyr::filter(to_cluster == 2 & to_data == "all") %>%
  dplyr::group_by(from_data) %>%
  dplyr::arrange(dplyr::desc(abs(cor))) %>%
  dplyr::slice_head(n = 3)


###find the cluster 4
cor_data %>%
  dplyr::filter(to_cluster == 4 & to_data == "all") %>%
  dplyr::group_by(from_data) %>%
  dplyr::arrange(dplyr::desc(abs(cor))) %>%
  dplyr::slice_head(n = 3)

###find the cluster 5
cor_data %>%
  dplyr::filter(to_cluster == 5 & to_data == "all") %>%
  dplyr::group_by(from_data) %>%
  dplyr::arrange(dplyr::desc(abs(cor))) %>%
  dplyr::slice_head(n = 3)


####Network to show the relationships
temp_data <-
  rbind(
    cor_data %>%
      dplyr::filter(to_cluster == 2 & to_data == "all") %>%
      dplyr::group_by(from_data) %>%
      dplyr::arrange(dplyr::desc(abs(cor))),
    cor_data %>%
      dplyr::filter(to_cluster == 4 & to_data == "all") %>%
      dplyr::group_by(from_data) %>%
      dplyr::arrange(dplyr::desc(abs(cor))),
    cor_data %>%
      dplyr::filter(to_cluster == 5 & to_data == "all") %>%
      dplyr::group_by(from_data) %>%
      dplyr::arrange(dplyr::desc(abs(cor)))
  ) %>%
  dplyr::filter(cor > 0 & p < 0.05)

edge_data <-
  temp_data %>%
  dplyr::mutate(
    from = paste(from_data, from_cluster, sep = '_'),
    to = paste(to_data, to_cluster, sep = '_')
  ) %>%
  dplyr::ungroup(from_data) %>%
  dplyr::select(from, to, cor, p)

node = unique(c(edge_data$from, edge_data$to))

node_data <-
  data.frame(node = unique(c(edge_data$from, edge_data$to)),
             class = stringr::str_replace(node, "_[0-9]{1,2}", "")) %>%
  dplyr::mutate(cluster = stringr::str_extract(node, "[0-9]{1,2}"))

library(ggraph)
library(igraph)
library(tidygraph)

network <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data)
####
V(network)$type <- bipartite_mapping(network)$type

coords <-
  create_layout(network, layout = "bipartite") %>%
  dplyr::select(x, y, everything())

coords$y[coords$class == "all"] = 0.3
coords$y[coords$class != "all"] = 1

table(coords$class)

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

network_new <-
  create_layout(
    graph = network,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot <-
  ggraph(network_new,
         layout = 'bipartite') +
  geom_edge_link(
    aes(width = cor),
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "grey"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(x = x,
        y = y,
        label = cluster),
    color = "black",
    bg.color = "white",
    size = 7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "female" = unname(sex_color["Female"]),
    "male" = unname(sex_color["Male"]),
    "all" = "black"
  )) +
  scale_size_manual(values = c(
    "all" = 12,
    "female" = 7,
    "male" = 7
  )) +
  scale_edge_width(range = c(0.5, 3)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot
library(extrafont)
extrafont::loadfonts()
# ggsave(plot = plot,
#        filename = "cluster_correlation.pdf",
#        width = 8.2,
#        height = 7)


temp_data %>% 
  dplyr::filter(to_data == "all" & from_data == "female") %>% 
  dplyr::filter(to_cluster == "2") %>% 
  dplyr::arrange(desc(abs(cor)))

temp_data %>% 
  dplyr::filter(to_data == "all" & from_data == "male") %>% 
  dplyr::filter(to_cluster == "2") %>% 
  dplyr::arrange(desc(abs(cor)))


temp_data %>% 
  dplyr::filter(to_data == "all" & from_data == "female") %>% 
  dplyr::filter(to_cluster == "4") %>% 
  dplyr::arrange(desc(abs(cor)))

temp_data %>% 
  dplyr::filter(to_data == "all" & from_data == "male") %>% 
  dplyr::filter(to_cluster == "4") %>% 
  dplyr::arrange(desc(abs(cor)))



temp_data %>% 
  dplyr::filter(to_data == "all" & from_data == "female") %>% 
  dplyr::filter(to_cluster == "5") %>% 
  dplyr::arrange(desc(abs(cor)))

temp_data %>% 
  dplyr::filter(to_data == "all" & from_data == "male") %>% 
  dplyr::filter(to_cluster == "5") %>% 
  dplyr::arrange(desc(abs(cor)))
