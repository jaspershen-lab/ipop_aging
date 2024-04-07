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

setwd("3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess")

object_cross_section_loess

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id)) %>%
  dplyr::arrange(age)

load("cor_data")

load("final_cluster_info")

head(final_cluster_info)

dim(final_cluster_info)

dim(object_cross_section_loess)

object_cross_section_loess@variable_info

variable_info <-
  object_cross_section_loess@variable_info %>%
  dplyr::select(-cluster)

variable_info <-
  variable_info %>%
  dplyr::left_join(cor_data, by = "variable_id")

sample_info <-
  object_cross_section_loess@sample_info

expression_data <-
  object_cross_section_loess@expression_data

final_cluster_info <-
  final_cluster_info %>%
  dplyr::left_join(variable_info, by = "variable_id")

library(org.Hs.eg.db)
library(clusterProfiler)

###Cluster 1
result_all_cluster1 <-
  tryCatch(
    expr = {
      load('cluster_1/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )

load('cluster_1/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster1 <- transcriptome_go
load('cluster_1/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster1 <- transcriptome_kegg
load('cluster_1/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster1 <- transcriptome_reactome

result_all_cluster1$pathway_id

###Cluster 2
result_all_cluster2 <-
  tryCatch(
    expr = {
      load('cluster_2/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )

result_all_cluster2 <- result_all
load('cluster_2/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster2 <- transcriptome_go
load('cluster_2/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster2 <- transcriptome_kegg
load('cluster_2/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster2 <- transcriptome_reactome

result_all_cluster2$pathway_id

###Cluster 3
result_all_cluster3 <-
  tryCatch(
    expr = {
      load('cluster_3/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_3/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster3 <- transcriptome_go
load('cluster_3/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster3 <- transcriptome_kegg
load('cluster_3/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster3 <- transcriptome_reactome

result_all_cluster3$pathway_id

###Cluster 4
result_all_cluster4 <-
  tryCatch(
    expr = {
      load('cluster_4/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_4/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster4 <- transcriptome_go
load('cluster_4/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster4 <- transcriptome_kegg
load('cluster_4/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster4 <- transcriptome_reactome

result_all_cluster4$pathway_id

###Cluster 5
result_all_cluster5 <-
  tryCatch(
    expr = {
      load('cluster_5/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_5/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster5 <- transcriptome_go
load('cluster_5/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster5 <- transcriptome_kegg
load('cluster_5/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster5 <- transcriptome_reactome

result_all_cluster5$pathway_id

###Cluster 6
result_all_cluster6 <-
  tryCatch(
    expr = {
      load('cluster_6/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_6/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster6 <- transcriptome_go
load('cluster_6/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster6 <- transcriptome_kegg
load('cluster_6/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster6 <- transcriptome_reactome

result_all_cluster6$pathway_id

###Cluster 7
result_all_cluster7 <-
  tryCatch(
    expr = {
      load('cluster_7/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_7/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster7 <- transcriptome_go
load('cluster_7/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster7 <- transcriptome_kegg
load('cluster_7/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster7 <- transcriptome_reactome

result_all_cluster7$pathway_id

###Cluster 8
result_all_cluster8 <-
  tryCatch(
    expr = {
      load('cluster_8/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_8/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster8 <- transcriptome_go
load('cluster_8/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster8 <- transcriptome_kegg
load('cluster_8/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster8 <- transcriptome_reactome

result_all_cluster8$pathway_id

###Cluster 9
result_all_cluster9 <-
  tryCatch(
    expr = {
      load('cluster_9/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_9/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster9 <- transcriptome_go
load('cluster_9/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster9 <- transcriptome_kegg
load('cluster_9/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster9 <- transcriptome_reactome

result_all_cluster9$pathway_id

###Cluster 10
result_all_cluster10 <-
  tryCatch(
    expr = {
      load('cluster_10/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_10/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster10 <- transcriptome_go
load('cluster_10/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster10 <- transcriptome_kegg
load('cluster_10/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster10 <- transcriptome_reactome

result_all_cluster10$pathway_id

###Cluster 11
result_all_cluster11 <-
  tryCatch(
    expr = {
      load('cluster_11/transcriptome_pathway/result_all')
      result_all
    },
    error = function(e) {
      NULL
    }
  )
load('cluster_11/transcriptome_pathway/GO_result/transcriptome_go')
transcriptome_go_cluster11 <- transcriptome_go
load('cluster_11/transcriptome_pathway/KEGG_result/transcriptome_kegg')
transcriptome_kegg_cluster11 <- transcriptome_kegg
load('cluster_11/transcriptome_pathway/Reactome_result/transcriptome_reactome')
transcriptome_reactome_cluster11 <- transcriptome_reactome

result_all_cluster11$pathway_id

###Cluter 3, 7 and 8 have no enrichmed pathways
result_all_cluster1$Description
result_all_cluster7$Description
result_all_cluster4$Description
result_all_cluster11$Description
result_all_cluster8$Description
result_all_cluster3$Description
result_all_cluster6$Description
result_all_cluster10$Description
result_all_cluster2$Description
result_all_cluster5$Description
result_all_cluster9$Description

result_all_list <-
  list(
    cluster1 = result_all_cluster1 %>% dplyr::arrange(p.adjust),
    # cluster4 = result_all_cluster4 %>% dplyr::arrange(p.adjust),
    cluster11 = result_all_cluster11 %>% dplyr::arrange(p.adjust),
    cluster6 = result_all_cluster6 %>% dplyr::arrange(p.adjust),
    cluster10 = result_all_cluster10 %>% dplyr::arrange(p.adjust),
    # cluster2 = result_all_cluster2 %>% dplyr::arrange(p.adjust),
    # cluster5 = result_all_cluster5 %>% dplyr::arrange(p.adjust),
    cluster9 = result_all_cluster9 %>% dplyr::arrange(p.adjust)
  )

word_cloud_list <-
  list(
    cluster1 = result_all_cluster1,
    # cluster4 = result_all_cluster4,
    cluster11 = result_all_cluster11,
    cluster6 = result_all_cluster6,
    cluster10 = result_all_cluster10,
    # cluster2 = result_all_cluster2,
    # cluster5 = result_all_cluster5,
    cluster9 = result_all_cluster9
  ) %>%
  purrr::map(function(x) {
    y <-
      stringr::str_split(x$Description, ";") %>%
      unlist() %>%
      stringr::str_split(pattern = " ") %>%
      unlist() %>%
      # stringr::str_to_lower() %>%
      stringr::str_replace_all("\\(", "") %>%
      stringr::str_replace_all("\\)", "") %>%
      stringr::str_replace_all("\\,", "") %>%
      stringr::str_replace_all("\\.", "") %>%
      stringr::str_replace_all("\\:", "")
    y <-
      combine_workds(y = y)
    y <- y[!y %in% remove_words] %>%
      sort()
    combine_workds(y = y)
  })

plot <-
  pahtway_heatmap_all_cluster(
    result_all_list = result_all_list,
    word_cloud_list = word_cloud_list,
    gene_marker = variable_info,
    expression_data = expression_data,
    top_pathway = 10,
    font_size = c(5, 10),
    word_count_breaks = c(1, 30),
    add_text = TRUE
  )

plot

# ggsave(plot,
#        filename = "transcriptome_pathway_heatmap.pdf",
#        width = 9,
#        height = 5)


####only cluster 2, cluster 4 and cluster 5

###Cluter 3, 7 and 8 have no enrichmed pathways
result_all_list <-
  list(
    cluster2 = result_all_cluster2 %>% dplyr::arrange(p.adjust),
    cluster4 = result_all_cluster4 %>% dplyr::arrange(p.adjust),
    cluster5 = result_all_cluster5 %>% dplyr::arrange(p.adjust)
  )

word_cloud_list <-
  list(cluster2 = result_all_cluster2,
       cluster4 = result_all_cluster4,
       cluster5 = result_all_cluster5) %>%
  purrr::map(function(x) {
    y <-
      stringr::str_split(x$Description, ";") %>%
      unlist() %>%
      stringr::str_split(pattern = " ") %>%
      unlist() %>%
      # stringr::str_to_lower() %>%
      stringr::str_replace_all("\\(", "") %>%
      stringr::str_replace_all("\\)", "") %>%
      stringr::str_replace_all("\\,", "") %>%
      stringr::str_replace_all("\\.", "") %>%
      stringr::str_replace_all("\\:", "")
    y <-
      combine_workds(y = y)
    y <- y[!y %in% remove_words] %>%
      sort()
    combine_workds(y = y)
  })

plot <-
  pahtway_heatmap_all_cluster(
    result_all_list = result_all_list,
    word_cloud_list = word_cloud_list,
    gene_marker = variable_info,
    expression_data = expression_data,
    top_pathway = 10,
    font_size = c(5, 10),
    word_count_breaks = c(1, 30),
    add_text = TRUE
  )

plot

# ggsave(plot,
#        filename = "transcriptome_pathway_heatmap2.pdf",
#        width = 9,
#        height = 3)


#####cluster 2
temp_data <-
  result_all_cluster2 %>%
  dplyr::arrange(p.adjust) %>%
  head(10)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = temp_data,
    result_kegg_cluster = NULL,
    result_reactome_cluster = NULL,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.1) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  temp_data

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)

subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

new_cluster <-
  purrr::map(cluster, function(x) {
    paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

new_cluster[new_cluster == "Module NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster)


node_result <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame()

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "result")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = node_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "cluster_2/transcriptome_pathway/node_result.xlsx",
             overwrite = TRUE)

###cluster different Reactome terms
result_all <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

library(plyr)
result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, fill = module)) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = module_annotation)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
  ) +
  ggsci::scale_fill_nejm()

plot

extrafont::loadfonts(device = "pdf")

# ggsave(plot,
#        filename = "cluster_2/transcriptome_pathway/network.pdf",
#        width = 10,
#        height = 6)


#####cluster 4
temp_data <-
  result_all_cluster4 %>%
  dplyr::arrange(p.adjust) %>%
  head(10)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = temp_data,
    result_kegg_cluster = NULL,
    result_reactome_cluster = NULL,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.1) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  temp_data

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)

subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

new_cluster <-
  purrr::map(cluster, function(x) {
    paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

new_cluster[new_cluster == "Module NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster)

node_result <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame()

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "result")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = node_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "cluster_4/transcriptome_pathway/node_result.xlsx",
             overwrite = TRUE)

###cluster different Reactome terms
result_all <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)



###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

library(plyr)
result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, fill = module)) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = module_annotation)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
  ) +
  ggsci::scale_fill_nejm()

plot

extrafont::loadfonts(device = "pdf")

# ggsave(plot,
#        filename = "cluster_4/transcriptome_pathway/network.pdf",
#        width = 10,
#        height = 6)

# write.csv(result_all, "cluster_4/transcriptome_pathway/result_all.csv")


#####cluster 5
temp_data <-
  result_all_cluster5 %>%
  dplyr::arrange(p.adjust) %>%
  head(10)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = temp_data,
    result_kegg_cluster = NULL,
    result_reactome_cluster = NULL,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.1) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  temp_data

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)

subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

new_cluster <-
  purrr::map(cluster, function(x) {
    paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

new_cluster[new_cluster == "Module NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster)

node_result <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame()

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "result")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = node_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "cluster_5/transcriptome_pathway/node_result.xlsx",
             overwrite = TRUE)

###cluster different Reactome terms
result_all <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

library(plyr)
result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()



plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, fill = module)) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = module_annotation)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
  ) +
  ggsci::scale_fill_nejm()

plot

extrafont::loadfonts(device = "pdf")

# ggsave(plot,
#        filename = "cluster_5/transcriptome_pathway/network.pdf",
#        width = 10,
#        height = 6)










#####cluster 10
temp_data <-
  result_all_cluster10 %>%
  dplyr::arrange(p.adjust) %>%
  head(10)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = temp_data,
    result_kegg_cluster = NULL,
    result_reactome_cluster = NULL,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.1) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  temp_data

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)

subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

new_cluster <-
  purrr::map(cluster, function(x) {
    paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

new_cluster[new_cluster == "Module NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster)

node_result <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame()

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "result")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = node_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "cluster_10/transcriptome_pathway/node_result.xlsx",
             overwrite = TRUE)

###cluster different Reactome terms
result_all <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

library(plyr)
result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, fill = module)) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = module_annotation)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
  ) +
  ggsci::scale_fill_nejm()

plot

extrafont::loadfonts(device = "pdf")

ggsave(plot,
       filename = "cluster_10/transcriptome_pathway/network.pdf",
       width = 10,
       height = 3)








#####cluster 11
temp_data <-
  result_all_cluster11 %>%
  dplyr::arrange(p.adjust) %>%
  head(10)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = temp_data,
    result_kegg_cluster = NULL,
    result_reactome_cluster = NULL,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.1) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  temp_data

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)

subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

new_cluster <-
  purrr::map(cluster, function(x) {
    paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

new_cluster[new_cluster == "Module NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster)

node_result <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame()

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "result")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = node_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "cluster_11/transcriptome_pathway/node_result.xlsx",
             overwrite = TRUE)

###cluster different Reactome terms
result_all <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

library(plyr)
result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, fill = module)) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = module_annotation)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
  ) +
  ggsci::scale_fill_nejm()

plot

extrafont::loadfonts(device = "pdf")

ggsave(plot,
       filename = "cluster_11/transcriptome_pathway/network.pdf",
       width = 10,
       height = 3)








#####cluster 6
temp_data <-
  result_all_cluster6 %>%
  dplyr::arrange(p.adjust) %>%
  head(10)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = temp_data,
    result_kegg_cluster = NULL,
    result_reactome_cluster = NULL,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.1) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  temp_data

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)

subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

new_cluster <-
  purrr::map(cluster, function(x) {
    paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

new_cluster[new_cluster == "Module NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster)

node_result <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame()

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "result")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = node_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "cluster_6/transcriptome_pathway/node_result.xlsx",
             overwrite = TRUE)

###cluster different Reactome terms
result_all <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

library(plyr)
result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, fill = module)) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = module_annotation)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
  ) +
  ggsci::scale_fill_nejm()

plot

extrafont::loadfonts(device = "pdf")

ggsave(plot,
       filename = "cluster_6/transcriptome_pathway/network.pdf",
       width = 10,
       height = 3)









#####cluster 9
temp_data <-
  result_all_cluster9 %>%
  dplyr::arrange(p.adjust) %>%
  head(10)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = temp_data,
    result_kegg_cluster = NULL,
    result_reactome_cluster = NULL,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.1) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  temp_data

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)

subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

new_cluster <-
  purrr::map(cluster, function(x) {
    paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

new_cluster[new_cluster == "Module NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster)

node_result <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame()

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "result")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

writeDataTable(
  wb,
  sheet = 1,
  x = node_result,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "cluster_9/transcriptome_pathway/node_result.xlsx",
             overwrite = TRUE)

###cluster different Reactome terms
result_all <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(module, p.adjust)

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

library(plyr)
result_all =
  result_all %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x %>% dplyr::select(-degree))
    } else{
      x =
        x %>%
        dplyr::arrange(p.adjust)
      x$module_annotation = x$module_annotation[1]
      x$Description = paste(x$Description, collapse = ";")
      x$p.adjust = x$p.adjust[1]
      x$database = paste(x$database, collapse = ";")
      x$geneID =
        x$geneID %>%
        stringr::str_split("/") %>%
        unlist() %>%
        unique() %>%
        paste(collapse = "/")
      x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
      x$pathway_id = paste(x$pathway_id, collapse = ";")
      x$class = paste(x$class, collapse = ";")
      x$degree = paste(x$degree, collapse = ";")
      x$module = x$module[1]
      x %>%
        dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
        dplyr::select(-degree)
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, fill = module)) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = module_annotation)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA),
  ) +
  ggsci::scale_fill_nejm()

plot

extrafont::loadfonts(device = "pdf")

ggsave(plot,
       filename = "cluster_9/transcriptome_pathway/network.pdf",
       width = 10,
       height = 3)

