no_source()
rm(list = ls())
library(tidyverse)
library(tidymass)
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

load("temp_data")

##remove the background
background <-
  temp_data %>%
  group_by(variable_id) %>%
  dplyr::summarise(n = sum(p_value_adjust < 0.05))

sum(background$n > length(unique(temp_data$center)) * 0.8)

background <-
  background %>%
  dplyr::filter(n > length(unique(temp_data$center)) * 0.8)

temp_data <-
  temp_data %>%
  dplyr::filter(!variable_id %in% background$variable_id)

temp_data <-
  temp_data %>%
  dplyr::left_join(variable_info,
                   by = c("variable_id"))

###Pathway enrichment
temp_data_proteomics <-
  temp_data %>%
  dplyr::filter(class == "proteomics" & p_value_adjust < 0.05)

####proteomics
dir.create("proteomics_pathway")
dir.create("proteomics_pathway/GO_result/")
dir.create("proteomics_pathway/KEGG_result/")
dir.create("proteomics_pathway/Reactome_result/")

####crest 1
temp_data_proteomics_crest1 <-
  temp_data_proteomics %>%
  dplyr::filter(center == 44)

library(org.Hs.eg.db)
library(clusterProfiler)

# ###GO enrichment
# proteomics_crest1_go <-
#   enrichGO(
#     gene = temp_data_proteomics_crest1$ENTREZID[!is.na(temp_data_proteomics_crest1$ENTREZID)],
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "fdr",
#     qvalueCutoff = 0.05
#   )
#
# proteomics_crest1_go@result <-
#   proteomics_crest1_go@result %>%
#   dplyr::filter(p.adjust < 0.05)
# save(proteomics_crest1_go, file = "proteomics_pathway/GO_result/proteomics_crest1_go")
load("proteomics_pathway/GO_result/proteomics_crest1_go")


#enrichment KEGG analysis
# proteomics_crest1_kegg <-
#   enrichKEGG(
#     gene = temp_data_proteomics_crest1$UNIPROT[!is.na(temp_data_proteomics_crest1$UNIPROT)],
#     organism = "hsa",
#     keyType = "uniprot",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "fdr",
#     qvalueCutoff = 0.05
#   )
#
# if(!is.null(proteomics_crest1_kegg)){
#   proteomics_crest1_kegg@result =
#     proteomics_crest1_kegg@result %>%
#     dplyr::filter(p.adjust < 0.05)
#
#   result <- proteomics_crest1_kegg@result
#
#   kegg_class <-
#     purrr::map(
#       .x = result$ID,
#       .f = function(x) {
#         temp <- KEGGREST::keggGet(dbentries = x)
#         class <- temp[[1]]$CLASS
#         if (is.null(class)) {
#           class <- "no"
#         }
#         class
#       }
#     )  %>%
#     unlist()
#
#   proteomics_crest1_kegg@result$kegg_class <- kegg_class
#
# }
#
# save(proteomics_crest1_kegg, file = "proteomics_pathway/KEGG_result/proteomics_crest1_kegg")
load("proteomics_pathway/KEGG_result/proteomics_crest1_kegg")

##reactome pathway enrichment
# library(ReactomePA)
# proteomics_crest1_reactome <-
#   ReactomePA::enrichPathway(gene = temp_data_proteomics_crest1$ENTREZID[!is.na(temp_data_proteomics_crest1$ENTREZID)],
#                             organism = "human",
#                             pvalueCutoff = 0.05,
#                             pAdjustMethod = "fdr",
#                             qvalueCutoff = 0.05)
#
# proteomics_crest1_reactome@result <-
#   proteomics_crest1_reactome@result %>%
#   dplyr::filter(p.adjust < 0.05)
#
# ##add disease class to Reactome pathways
# disease <-
#   purrr::map(proteomics_crest1_reactome@result$ID, function(x) {
#     temp <- ReactomeContentService4R::getPathways(x)$isInDisease
#     if(is.null(temp)){
#       return(FALSE)
#     }
#     temp
#   }) %>%
#   unlist()
#
# proteomics_crest1_reactome@result$disease <-
#   disease
#
# save(proteomics_crest1_reactome, file = "proteomics_pathway/Reactome_result/proteomics_crest1_reactome")
load("proteomics_pathway/Reactome_result/proteomics_crest1_reactome")

##------------------------------------------------------------------------------
##If want to repeat, from here
##-------------------------------------------------------------------------------
##get the enriched results and only remain the pathways with P < 0.05 and Count >= 3
####GO
load("proteomics_pathway/GO_result/proteomics_crest1_go")
proteomics_crest1_go <-
  proteomics_crest1_go@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(Count >= 3) %>%
  dplyr::arrange(p.adjust)

###KEGG
load("proteomics_pathway/KEGG_result/proteomics_crest1_kegg")
if (!is.null(proteomics_crest1_kegg)) {
  proteomics_crest1_kegg <-
    proteomics_crest1_kegg@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::filter(Count >= 3)
}

####enrichment reactome
load("proteomics_pathway/Reactome_result/proteomics_crest1_reactome")

proteomics_crest1_reactome <-
  proteomics_crest1_reactome@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  # dplyr::filter(!disease) %>%
  dplyr::filter(Count >= 3) %>%
  dplyr::arrange(p.adjust)

# #######get the similarity between GO terms
# go_sim_matrix =
#   get_go_result_sim(result = proteomics_crest1_go, sim.cutoff = 0)
#
# save(go_sim_matrix, file = "proteomics_pathway/GO_result/go_sim_matrix")

load("proteomics_pathway/GO_result/go_sim_matrix")

go_sim_matrix =
  go_sim_matrix %>%
  dplyr::filter(sim > 0.7)

edge_data <-
  rbind(go_sim_matrix) %>%
  dplyr::rename(from = name1, to = name2)

node_data <-
  rbind(proteomics_crest1_go) %>%
  as.data.frame() %>%
  dplyr::select(ID, everything()) %>%
  dplyr::rename(node = ID)

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)
library(igraph)

# subnetwork <-
#   igraph::cluster_edge_betweenness(graph = temp_graph,
#                                    weights = abs(edge_attr(temp_graph,
#                                                            "sim")))
#
# save(subnetwork, file = "proteomics_pathway/GO_result/subnetwork")
load("proteomics_pathway/GO_result/subnetwork")

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

library(tidygraph)

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(module = new_cluster) %>%
  activate(what = "nodes")
# dplyr::filter(module != "Other")

###clustered different GO terms
result_go <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(ONTOLOGY, module, p.adjust)

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "cluster 1")
#
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
#
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = result_go,
#   colNames = TRUE,
#   rowNames = FALSE
# )
#
# saveWorkbook(wb, "proteomics_pathway/GO_result/result_go.xlsx", overwrite = TRUE)

###plot to show the clusters of GO terms
cluster_label1 =
  igraph::as_data_frame(temp_graph, what = "vertices") %>%
  # dplyr::filter(module != "Other") %>%
  dplyr::group_by(module) %>%
  dplyr::filter(p.adjust == min(p.adjust)) %>%
  pull(Description)

cluster_label2 =
  igraph::as_data_frame(temp_graph, what = "vertices") %>%
  dplyr::filter(module == "Other") %>%
  pull(Description)

cluster_label = c(cluster_label1, cluster_label2)

plot <-
  temp_graph %>%
  tidygraph::filter(module != "Other") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(fill = module,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_node_text(aes(
    x = x,
    y = y,
    label = ifelse(Description %in% cluster_label, Description, NA)
  ),
  size = 3,
  repel = TRUE) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 7)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

# ggplot2::ggsave(plot,
#        filename =
#          "proteomics_pathway/GO_result/GO_sim_plot.pdf",
#        width = 7,
#        height = 7)

dir.create("proteomics_pathway/GO_result/GO_sim_plot")

# for(temp_cluster in unique(new_cluster)) {
#   cat(temp_cluster, " ")
#   if(temp_cluster == "Other"){
#     next()
#   }
#   plot1 <-
#     temp_graph %>%
#     tidygraph::filter(module == temp_cluster) %>%
#     ggraph(layout = 'kk',
#            circular = FALSE) +
#     geom_edge_link(
#       aes(width = sim),
#       strength = 1,
#       color = "black",
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     geom_node_point(
#       aes(fill = -log(p.adjust, 10),
#           size = Count),
#       shape = 21,
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     shadowtext::geom_shadowtext(aes(x = x,
#                                     y = y,
#                                     label = ifelse(module == "Other", NA, Description)),
#                                 check_overlap = TRUE,
#                                 size = 3,
#                                 color = "black",
#                                 bg.color = "white") +
#     guides(fill = guide_legend(ncol = 1)) +
#     scale_edge_width_continuous(range = c(0.1, 2)) +
#     scale_size_continuous(range = c(3, 10)) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "left",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     )
#
#   # plot1
#
#   plot2 <-
#     result_go %>%
#     dplyr::filter(module == temp_cluster) %>%
#     dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#     dplyr::arrange(p.adjust) %>%
#     dplyr::mutate(Description = factor(Description, levels = Description)) %>%
#     ggplot(aes(p.adjust, Description)) +
#     geom_bar(stat = "identity", color = "grey") +
#     geom_text(aes(x = 0, Description, label = Description),
#               hjust = 0,
#               size = 5, color = "red") +
#     theme_bw() +
#     labs(y = "", x = "-log10(FDR adjusted P value)") +
#     scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
#     theme(
#       panel.grid = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank()
#     )
#
#   library(ggwordcloud)
#
#   temp_data =
#     result_go %>%
#     dplyr::filter(module == temp_cluster) %>%
#     dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#     dplyr::select(Description, p.adjust) %>%
#     dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
#     plyr::dlply(.variables = .(Description)) %>%
#     purrr::map(function(x) {
#       data.frame(word = stringr::str_split(x$Description, " ")[[1]],
#                  p.adjust = x$p.adjust)
#     }) %>%
#     do.call(rbind, .) %>%
#     as.data.frame() %>%
#     plyr::dlply(.variables = .(word)) %>%
#     purrr::map(function(x) {
#       x$p.adjust <- sum(x$p.adjust)
#       x %>%
#         dplyr::distinct(word, .keep_all = TRUE)
#     }) %>%
#     do.call(rbind, .) %>%
#     as.data.frame() %>%
#     dplyr::filter(!word %in% remove_words)
#
#   plot3 =
#     temp_data %>%
#     ggplot(aes(label = word, size = p.adjust)) +
#     geom_text_wordcloud() +
#     scale_radius(range = c(5, 15), limits = c(0, NA)) +
#     theme_minimal()
#
#   library(patchwork)
#   plot =
#     plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)
#
#   ggsave(
#     plot,
#     filename = file.path(
#       "proteomics_pathway/GO_result/GO_sim_plot",
#       paste(temp_cluster, "sim_plot.pdf", sep = "_")
#     ),
#     width = 21,
#     height = 7
#   )
# }

##matrix tow show the cluster GO terms
result_go %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(n = n * 10 / max(n) + 2)

library(igraph)
library(simplifyEnrichment)

# for(ont in c('MF', "BP", "CC")) {
#   cat(ont, " ")
#   show_matrix_cluster(
#     result = result_go %>% dplyr::mutate(Direction = "UP") %>% dplyr::rename(cluster = module),
#     ont = ont,
#     measure = "Wang",
#     remove_words = remove_words,
#     margin = 15,
#     width = 16,
#     height = 14,
#     path = "proteomics_pathway/GO_result",
#     top = 15
#   )
# }

#####output GO result as clusters
library(plyr)
result_go_cluster <-
  result_go %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }
    
    if (x$module[1] == 'Other') {
      return(x)
    }
    
    x =
      x %>%
      dplyr::arrange(p.adjust)
    
    x$node <-
      paste(x$node, collapse = ";")
    
    x$Description <-
      paste(x$Description, collapse = ";")
    
    x$BgRatio <-
      paste(x$BgRatio, collapse = ";")
    
    x$pvalue <- min(as.numeric(x$pvalue))
    x$p.adjust <- min(as.numeric(x$p.adjust))
    x$qvalue <- min(as.numeric(x$qvalue))
    x$geneID =
      x$geneID %>%
      stringr::str_split(pattern = "/") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = '/')
    
    x$Count <-
      length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
    
    x =
      x %>%
      dplyr::select(module, everything()) %>%
      dplyr::distinct(module, .keep_all = TRUE)
    
    x
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(
    module_annotation = case_when(
      module == "Other" ~ Description,
      module != "Other" ~ stringr::str_split(Description, ";")[[1]][1]
    )
  ) %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(module_annotation, everything())

result_go_cluster$module_annotation =
  stringr::str_split(result_go_cluster$Description, ";") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "cluster 1")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
#
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = result_go_cluster,
#   colNames = TRUE,
#   rowNames = FALSE
# )
#
# saveWorkbook(wb, "proteomics_pathway/GO_result/result_go_cluster.xlsx", overwrite = TRUE)

###output the cluster annotation for each cluster
dir.create("proteomics_pathway/GO_result/GO_graph")
# unique(result_go_cluster$module) %>%
#   purrr::map(
#     .f = function(x) {
#       cat(x, " ")
#       if (x == "Other") {
#         return(NULL)
#       }
#
#       temp_id =
#         result_go_cluster %>%
#         dplyr::filter(module == x) %>%
#         dplyr::pull(node) %>%
#         stringr::str_split(";") %>%
#         `[[`(1) %>%
#         pRoloc::goIdToTerm(keepNA = FALSE) %>%
#         data.frame(id = ., class = "YES") %>%
#         tibble::rownames_to_column(var = "name")
#
#       temp_plot =
#         GOSim::getGOGraph(term = temp_id$name, prune = Inf) %>%
#         igraph::igraph.from.graphNEL() %>%
#         tidygraph::as_tbl_graph() %>%
#         tidygraph::left_join(temp_id, by = "name") %>%
#         dplyr::mutate(class = case_when(is.na(class) ~ "NO",
#                                         TRUE ~ class))
#
#       plot =
#         temp_plot %>%
#         ggraph(layout = 'kk',
#                circular = FALSE) +
#         geom_edge_link(
#           color = ggsci::pal_aaas()(n = 10)[1],
#           alpha = 1,
#           arrow = grid::arrow(
#             angle = 10,
#             length = unit(0.2, "inches"),
#             type = "closed"
#           ),
#           show.legend = FALSE
#         ) +
#         geom_node_point(
#           aes(fill = class),
#           shape = 21,
#           alpha = 1,
#           size = 6,
#           show.legend = FALSE
#         ) +
#         geom_node_text(aes(
#           x = x,
#           y = y,
#           label = ifelse(class == "YES", id, NA)
#         ),
#         size = 3,
#         repel = TRUE) +
#         scale_fill_manual(values = c('YES' = "red", 'NO' = "white")) +
#         ggraph::theme_graph() +
#         theme(
#           plot.background = element_rect(fill = "transparent", color = NA),
#           panel.background = element_rect(fill = "transparent", color = NA),
#           legend.position = "left",
#           legend.background = element_rect(fill = "transparent", color = NA)
#         )
#       plot
#       ggsave(
#         plot,
#         filename = file.path(
#           "proteomics_pathway/GO_result/GO_graph",
#           paste(x, "_GO graph.pdf", sep = "")
#         ),
#         width = 7,
#         height = 7
#       )
#     }
#   )

# ######-----------------------------------------------------------------------
### KEGG
###calculate the similarity of KEGG pathway
# if(!is.null(proteomics_crest1_kegg)){
#   kegg_sim_matrix <-
#     simplifyEnrichment::term_similarity_from_KEGG(term_id = c(proteomics_crest1_kegg$ID),
#                                                   method = "jaccard") %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column(var = "name1") %>%
#     tidyr::pivot_longer(cols = -name1,
#                         names_to = "name2",
#                         values_to = "sim") %>%
#     dplyr::filter(name1 != name2)
# }else{
#   kegg_sim_matrix <- NULL
# }
#
# save(kegg_sim_matrix, file = "proteomics_pathway/KEGG_result/kegg_sim_matrix")

load("proteomics_pathway/KEGG_result/kegg_sim_matrix")

if (is.null(kegg_sim_matrix)) {
  kegg_sim_matrix <-
    data.frame(name1 = "a",
               name2 = "b",
               sim = 0)
}

kegg_sim_matrix =
  kegg_sim_matrix %>%
  dplyr::filter(sim > 0.5)

edge_data <-
  kegg_sim_matrix %>%
  dplyr::rename(from = name1, to = name2)

node_data <-
  rbind(proteomics_crest1_kegg) %>%
  as.data.frame() %>%
  dplyr::select(ID, everything()) %>%
  dplyr::rename(node = ID) %>%
  dplyr::filter(node %in% c(edge_data$from, edge_data$to))

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

library(igraph)

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

####clustered different KEGG terms
result_kegg <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(p.adjust)

# save(result_kegg, file = "proteomics_pathway/KEGG_result/result_kegg")
load("proteomics_pathway/KEGG_result/result_kegg")

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "cluster 1")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
#
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = result_kegg,
#   colNames = TRUE,
#   rowNames = FALSE
# )
#
# saveWorkbook(wb, "KEGG_result/result_kegg.xlsx", overwrite = TRUE)

plot <-
  temp_graph %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = cluster,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_text(aes(
    x = x,
    y = y,
    label = ifelse(cluster == "Other", Description, Description)
  ),
  size = 3,
  repel = TRUE) +
  guides(fill = guide_legend(ncol = 1)) +
  # scale_fill_manual(values = c(
  #   "Other" = "grey"
  # )) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 7)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
# facet_nodes(Direction~ONTOLOGY)

plot

# ggsave(plot,
#        filename = "proteomics_pathway/KEGG_result/kegg_sim_plot.pdf",
#        width = 9,
#        height = 7)


dir.create("proteomics_pathway/KEGG_result/KEGG_sim_plot")

# for(temp_cluster in unique(new_cluster)){
#   cat(temp_cluster, " ")
#   plot1 <-
#     temp_graph %>%
#     tidygraph::filter(module == temp_cluster) %>%
#     ggraph(layout = 'fr',
#            circular = FALSE) +
#     geom_edge_link(
#       aes(width = sim),
#       strength = 1,
#       color = "black",
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     geom_node_point(
#       aes(fill = cluster,
#           size = -log(p.adjust, 10)),
#       shape = 21,
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     geom_node_text(aes(
#       x = x,
#       y = y,
#       label = ifelse(module == "Other", NA, Description)
#     ),
#     size = 3,
#     repel = TRUE) +
#     guides(fill = guide_legend(ncol = 1)) +
#     scale_edge_width_continuous(range = c(0.1,2)) +
#     scale_size_continuous(range = c(1, 7)) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "left",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     )
#
#   plot1
#
#   plot2 <-
#   result_kegg %>%
#     dplyr::filter(module == temp_cluster) %>%
#     dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#     dplyr::arrange(p.adjust) %>%
#     dplyr::mutate(Description = factor(Description, levels = Description)) %>%
#     ggplot(aes(p.adjust, Description)) +
#     geom_bar(stat = "identity") +
#     geom_text(aes(x = 0, Description, label = Description),
#               hjust = 0,
#               size = 5) +
#     theme_bw() +
#     labs(y = "", x = "-log10(FDR adjusted P value)") +
#     scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
#     theme(
#       panel.grid = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank()
#     )
#
#   library(ggwordcloud)
#
#   temp_data =
#     result_kegg %>%
#     dplyr::filter(module == temp_cluster) %>%
#     dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#     dplyr::select(Description, p.adjust) %>%
#     dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
#     plyr::dlply(.variables = .(Description)) %>%
#     purrr::map(function(x){
#       data.frame(word = stringr::str_split(x$Description, " ")[[1]],
#                  p.adjust = x$p.adjust)
#     }) %>%
#     do.call(rbind, .) %>%
#     as.data.frame() %>%
#     plyr::dlply(.variables = .(word)) %>%
#     purrr::map(function(x){
#       x$p.adjust <- sum(x$p.adjust)
#       x %>%
#         dplyr::distinct(word, .keep_all = TRUE)
#     }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#     dplyr::filter(!word %in% remove_words)
#
#   plot3 =
#     temp_data %>%
#     ggplot(aes(label = word, size = p.adjust)) +
#     geom_text_wordcloud() +
#     scale_radius(range = c(5, 15), limits = c(0, NA)) +
#     theme_minimal()
#
#   library(patchwork)
#   plot =
#     plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)
#
#   ggsave(
#     plot,
#     filename = file.path(
#       "proteomics_pathway/KEGG_result/KEGG_sim_plot",
#       paste(temp_cluster, "sim_plot.pdf", sep = "_")
#     ),
#     width = 21,
#     height = 7
#   )
# }


##matrix tow show the cluster KEGG terms

# show_matrix_cluster(
#   result = result_kegg %>% dplyr::mutate(Direction = "UP") %>% dplyr::rename(cluster = module),
#   ont = NULL,
#   measure = 'jaccard',
#   remove_words = remove_words,
#   margin = 15,
#   height = 8,
#   path = "KEGG_result"
# )

#####output KEGG result
result_kegg_cluster <-
  result_kegg %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }
    
    if (x$cluster[1] == 'Other') {
      return(x)
    }
    
    x =
      x %>%
      dplyr::arrange(p.adjust)
    
    x$node <-
      paste(x$node, collapse = ";")
    
    x$Description <-
      paste(x$Description, collapse = ";")
    
    x$BgRatio <-
      paste(x$BgRatio, collapse = ";")
    
    x$pvalue <- min(as.numeric(x$pvalue))
    x$p.adjust <- min(as.numeric(x$p.adjust))
    x$qvalue <- min(as.numeric(x$qvalue))
    x$geneID =
      x$geneID %>%
      stringr::str_split(pattern = "/") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = '/')
    
    x$Count <-
      length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
    
    x =
      x %>%
      dplyr::select(module, everything()) %>%
      dplyr::distinct(module, .keep_all = TRUE)
    
    x
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(module_annotation = case_when(module == "Other" ~ Description,
                                              module != "Other" ~ cluster)) %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(module_annotation, everything())

result_kegg_cluster$module_annotation =
  stringr::str_split(result_kegg_cluster$Description, ";") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "cluster 1")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
#
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = result_kegg_cluster,
#   colNames = TRUE,
#   rowNames = FALSE
# )
#
# saveWorkbook(wb, "proteomics_pathway/KEGG_result/result_kegg_cluster.xlsx", overwrite = TRUE)



#########-----------------------------------------------------------------------
cat("reactome")
####Reactome
##calculate the similarity of Reactome pathway
# reactome_sim_matrix <-
#   tryCatch(
#     expr =  simplifyEnrichment::term_similarity_from_Reactome(term_id = c(proteomics_crest1_reactome$ID),
#                                                                method = "jaccard") %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column(var = "name1") %>%
#       tidyr::pivot_longer(
#         cols = -name1,
#         names_to = "name2",
#         values_to = "sim"
#       ) %>%
#       dplyr::filter(name1 != name2),
#     error = function(e) {
#       data.frame(name1 = character(), name2 = character(),
#                  sim = numeric())
#     }
#   )
#
# save(reactome_sim_matrix, file = "proteomics_pathway/Reactome_result/reactome_sim_matrix")

load("proteomics_pathway/Reactome_result/reactome_sim_matrix")

reactome_sim_matrix =
  reactome_sim_matrix %>%
  dplyr::filter(sim > 0.5)

edge_data <-
  reactome_sim_matrix %>%
  dplyr::rename(from = name1, to = name2)

node_data <-
  rbind(proteomics_crest1_reactome) %>%
  as.data.frame() %>%
  dplyr::select(ID, everything()) %>%
  dplyr::rename(node = ID)

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
  tidygraph::mutate(module = new_cluster) %>%
  activate(what = "nodes")
# dplyr::filter(module != "Other")

# ###cluster different Reactome terms
# result_reactome <-
#   igraph::vertex_attr(temp_graph) %>%
#   do.call(cbind, .) %>%
#   as.data.frame() %>%
#   dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#   dplyr::arrange(p.adjust)
#
# save(result_reactome, file = "proteomics_pathway/Reactome_result/result_reactome")
load("proteomics_pathway/Reactome_result/result_reactome")

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "cluster 1")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
#
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = result_reactome,
#   colNames = TRUE,
#   rowNames = FALSE
# )
#
# saveWorkbook(wb, "proteomics_pathway/Reactome_result/result_reactome.xlsx", overwrite = TRUE)

cluster_label1 =
  igraph::as_data_frame(temp_graph, what = "vertices") %>%
  dplyr::filter(module != "Other") %>%
  dplyr::group_by(module) %>%
  dplyr::filter(p.adjust == min(p.adjust)) %>%
  pull(Description)

cluster_label2 =
  igraph::as_data_frame(temp_graph, what = "vertices") %>%
  dplyr::filter(module == "Other") %>%
  pull(Description)

cluster_label = c(cluster_label1, cluster_label2)

plot <-
  temp_graph %>%
  tidygraph::filter(module != "Other") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
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
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(
        Description %in% cluster_label,
        stringr::str_replace(Description, "Homo sapiens\\\r: ", ""),
        NA
      )
    ),
    check_overlap = TRUE,
    size = 3,
    color = "black",
    bg.color = "white"
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 7)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
# facet_nodes(Direction~ONTOLOGY)

plot

# ggsave(plot,
#        filename =
#          "proteomics_pathway/Reactome_result/Reactome_sim_plot.pdf",
#        width = 8.5,
#        height = 7)

dir.create("proteomics_pathway/Reactome_result/Reactome_sim_plot")

# for(temp_cluster in unique(new_cluster)) {
#   cat(temp_cluster, " ")
#   if(temp_cluster == "Other"){
#     next()
#   }
#   plot1 <-
#     temp_graph %>%
#     tidygraph::filter(module == temp_cluster) %>%
#     ggraph(layout = 'fr',
#            circular = FALSE) +
#     geom_edge_link(
#       aes(width = sim),
#       strength = 1,
#       color = "black",
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     geom_node_point(
#       aes(fill = -log(p.adjust, 10),
#           size = Count),
#       shape = 21,
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     shadowtext::geom_shadowtext(aes(x = x,
#                                     y = y,
#                                     label = ifelse(module == "Other", NA, Description)),
#                                 check_overlap = TRUE,
#                                 size = 3,
#                                 color = "black",
#                                 bg.color = "white") +
#     guides(fill = guide_legend(ncol = 1)) +
#     scale_edge_width_continuous(range = c(0.1, 2)) +
#     scale_size_continuous(range = c(3, 10)) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "left",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     )
#
#   plot1
#
# plot2 <-
#     result_reactome %>%
#   dplyr::filter(module == temp_cluster) %>%
#   dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#   dplyr::arrange(p.adjust) %>%
#   dplyr::mutate(Description = factor(Description, levels = Description)) %>%
#   ggplot(aes(p.adjust, Description)) +
#   geom_bar(stat = "identity", color = "grey") +
#   geom_text(aes(x = 0, Description, label = Description),
#             hjust = 0,
#             size = 5, color = "red") +
#   theme_bw() +
#   labs(y = "", x = "-log10(FDR adjusted P value)") +
#   scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )
#
#   library(ggwordcloud)
#
#   temp_data =
#     result_reactome %>%
#     dplyr::filter(module == temp_cluster) %>%
#     dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#     dplyr::select(Description, p.adjust) %>%
#     dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
#     dplyr::mutate(Description = stringr::str_replace_all(Description, "\\(", "")) %>%
#     dplyr::mutate(Description = stringr::str_replace_all(Description, "\\)", "")) %>%
#     dplyr::mutate(Description = stringr::str_replace_all(Description, "\\:", "")) %>%
#     plyr::dlply(.variables = .(Description)) %>%
#     purrr::map(function(x) {
#       data.frame(word = stringr::str_split(x$Description, " ")[[1]],
#                  p.adjust = x$p.adjust)
#     }) %>%
#     do.call(rbind, .) %>%
#     as.data.frame() %>%
#     plyr::dlply(.variables = .(word)) %>%
#     purrr::map(function(x) {
#       x$p.adjust <- sum(x$p.adjust)
#       x %>%
#         dplyr::distinct(word, .keep_all = TRUE)
#     }) %>%
#     do.call(rbind, .) %>%
#     as.data.frame() %>%
#     dplyr::filter(!word %in% remove_words)
#
#   plot3 =
#     temp_data %>%
#     ggplot(aes(label = word, size = p.adjust)) +
#     geom_text_wordcloud() +
#     scale_radius(range = c(5, 15), limits = c(0, NA)) +
#     theme_minimal()
#
#   library(patchwork)
#   plot =
#     plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)
#
#   ggsave(
#     plot,
#     filename = file.path(
#       "proteomics_pathway/Reactome_result/Reactome_sim_plot",
#       paste(temp_cluster, "sim_plot.pdf", sep = "_")
#     ),
#     width = 21,
#     height = 7
#   )
# }

#matrix tow show the cluster Reactome terms
library(simplifyEnrichment)

# show_matrix_cluster(
#   result = result_reactome %>% dplyr::mutate(Direction = "UP") %>% dplyr::rename(cluster = module),
#   ont = NULL,
#   measure = "jaccard",
#   remove_words = remove_words,
#   margin = 15,
#   width = 16,
#   height = 12,
#   path = "proteomics_pathway/Reactome_result"
# )

#####output Reactome result
result_reactome_cluster <-
  result_reactome %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    }
    
    if (x$module[1] == 'Other') {
      return(x)
    }
    
    x =
      x %>%
      dplyr::arrange(p.adjust)
    
    x$node <-
      paste(x$node, collapse = ";")
    
    x$Description <-
      paste(x$Description, collapse = ";")
    
    x$BgRatio <-
      paste(x$BgRatio, collapse = ";")
    
    x$pvalue <- min(as.numeric(x$pvalue))
    x$p.adjust <- min(as.numeric(x$p.adjust))
    x$qvalue <- min(as.numeric(x$qvalue))
    x$geneID =
      x$geneID %>%
      stringr::str_split(pattern = "/") %>%
      unlist() %>%
      unique() %>%
      paste(collapse = '/')
    
    x$Count <-
      length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
    
    x =
      x %>%
      dplyr::select(module, everything()) %>%
      dplyr::distinct(module, .keep_all = TRUE)
    
    x
    
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(module_annotation = case_when(module == "Other" ~ Description,
                                              module != "Other" ~ module)) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(module_annotation, everything())

result_reactome_cluster$module_annotation =
  stringr::str_split(result_reactome_cluster$Description, ";") %>%
  purrr::map(function(x) {
    x[1]
  }) %>%
  unlist()

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "cluster 1")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
#
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = result_reactome_cluster,
#   colNames = TRUE,
#   rowNames = FALSE
# )
#
# saveWorkbook(
#   wb,
#   "proteomics_pathway/Reactome_result/result_reactome_cluster.xlsx",
#   overwrite = TRUE
# )

####---------------------------------------------------------------------------
########combine three databases together
####barplot
cat("stop here")

result_go_cluster <-
  rbind(
    readxl::read_xlsx(
      "proteomics_pathway/GO_result/result_go_cluster.xlsx",
      sheet = 1
    ) %>%
      dplyr::filter(ONTOLOGY != "CC") %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(database = "GO") %>%
      dplyr::select(
        module_annotation,
        Description,
        p.adjust,
        Count,
        database,
        geneID,
        pathway_id = node
      )
  ) %>%
  dplyr::mutate(Count = as.numeric(Count)) %>%
  dplyr::filter(!is.na(module_annotation))

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
addWorksheet(wb, sheetName = "cluster 1")
freezePane(
  wb = wb,
  sheet = 1,
  firstRow = TRUE,
  firstCol = TRUE
)

proteomics_crest1_kegg

result_kegg_cluster <-
  proteomics_crest1_kegg

result_kegg_cluster <-
  result_kegg_cluster %>%
  dplyr::mutate(module_annotation = Description,
                database = "KEGG") %>%
  dplyr::rename(pathway_id = ID) %>%
  dplyr::select(module_annotation,
                Description,
                p.adjust,
                Count,
                database,
                geneID,
                pathway_id)

writeDataTable(
  wb,
  sheet = 1,
  x = result_kegg_cluster,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(wb,
             "proteomics_pathway/KEGG_result/result_kegg_cluster.xlsx",
             overwrite = TRUE)

# result_kegg_cluster <-
#   rbind(
#     readxl::read_xlsx(
#       "proteomics_pathway/KEGG_result/result_kegg_cluster.xlsx",
#       sheet = 1
#     ) %>%
#       dplyr::arrange(p.adjust) %>%
#       dplyr::mutate(database = "KEGG") %>%
#       dplyr::select(
#         module_annotation,
#         Description,
#         p.adjust,
#         Count,
#         database,
#         geneID,
#         pathway_id = node
#       )
#   ) %>%
#   dplyr::mutate(Count = as.numeric(Count)) %>%
#   dplyr::filter(!is.na(module_annotation))

result_reactome_cluster <-
  rbind(
    readxl::read_xlsx(
      "proteomics_pathway/Reactome_result/result_reactome_cluster.xlsx",
      sheet = 1
    ) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(database = "Reactome") %>%
      dplyr::select(
        module_annotation,
        Description,
        p.adjust,
        Count,
        database,
        geneID,
        pathway_id = node
      )
  ) %>%
  dplyr::mutate(Count = as.numeric(Count)) %>%
  dplyr::filter(!is.na(module_annotation)) %>%
  dplyr::mutate(module_annotation = stringr::str_replace(module_annotation, "Homo sapiens\\\r: ", ""))

dim(result_go_cluster)
dim(result_kegg_cluster)
dim(result_reactome_cluster)

######calculate the similarity (jaccard index) between all the pathways
jaccard_index <-
  get_jaccard_index_for_three_databases(
    result_go_cluster = result_go_cluster,
    result_kegg_cluster = result_kegg_cluster,
    result_reactome_cluster = result_reactome_cluster,
    variable_info = variable_info
  )

head(jaccard_index)

edge_data <-
  jaccard_index %>%
  dplyr::filter(value > 0.5) %>%
  dplyr::rename(from = name1, to = name2, sim = value)

node_data <-
  rbind(
    data.frame(result_go_cluster, class = "GO"),
    # data.frame(result_kegg_cluster, class = "KEGG"),
    data.frame(result_reactome_cluster, class = "Reactome")
  )
# dplyr::filter(module_annotation %in% c(edge_data$from, edge_data$to))

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

###cluster different Reactome terms
# result_all <-
#   igraph::vertex_attr(temp_graph) %>%
#   do.call(cbind, .) %>%
#   as.data.frame() %>%
#   dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#   dplyr::arrange(module, p.adjust)
#
# ##result_all is the pahtway enrichment for all three databases
# save(result_all, file = "proteomics_pathway/result_all")
load("proteomics_pathway/result_all")

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

# library(plyr)
# result_all =
#   result_all %>%
#   plyr::dlply(.variables = .(module)) %>%
#   purrr::map(function(x) {
#     if (nrow(x) == 1) {
#       return(x %>% dplyr::select(-degree))
#     } else{
#       x =
#         x %>%
#         dplyr::arrange(p.adjust)
#       x$module_annotation = x$module_annotation[1]
#       x$Description = paste(x$Description, collapse = ";")
#       x$p.adjust = x$p.adjust[1]
#       x$database = paste(x$database, collapse = ";")
#       x$geneID =
#         x$geneID %>%
#         stringr::str_split("/") %>%
#         unlist() %>%
#         unique() %>%
#         paste(collapse = "/")
#       x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
#       x$pathway_id = paste(x$pathway_id, collapse = ";")
#       x$class = paste(x$class, collapse = ";")
#       x$degree = paste(x$degree, collapse = ";")
#       x$module = x$module[1]
#       x %>%
#         dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
#         dplyr::select(-degree)
#     }
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
# addWorksheet(wb, sheetName = "cluster 1")
# freezePane(
#   wb = wb,
#   sheet = 1,
#   firstRow = TRUE,
#   firstCol = TRUE
# )
#
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = result_all,
#   colNames = TRUE,
#   rowNames = FALSE
# )
#
# saveWorkbook(wb, "proteomics_pathway/result_all.xlsx", overwrite = TRUE)

cluster_label1 <-
  igraph::as_data_frame(temp_graph, what = "vertices") %>%
  dplyr::filter(module != "Other") %>%
  dplyr::group_by(module) %>%
  dplyr::filter(p.adjust == min(p.adjust)) %>%
  pull(module_annotation)

cluster_label2 <-
  igraph::as_data_frame(temp_graph, what = "vertices") %>%
  dplyr::filter(module == "Other") %>%
  pull(module_annotation)

cluster_label = c(cluster_label1, cluster_label2)

plot <-
  temp_graph %>%
  activate(what = "nodes") %>%
  tidygraph::filter(module != "Other") %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = class,
        size = -log(p.adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(module_annotation %in% cluster_label1, module_annotation, NA)
    ),
    check_overlap = TRUE,
    size = 3,
    color = "black",
    bg.color = "white"
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 10)) +
  scale_fill_manual(values = database_color) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
  )

plot

extrafont::loadfonts(device = "pdf")

ggsave(plot,
       filename = "proteomics_pathway/GO_KEGG_Reactome_sim_plot.pdf",
       width = 6,
       height = 7)

result_all
