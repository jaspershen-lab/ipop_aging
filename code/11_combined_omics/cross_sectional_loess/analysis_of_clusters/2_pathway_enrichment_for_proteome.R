no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("code/tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

# 1:11 %>%
#   purrr::walk(function(i) {
#     message(i)
#     name <-
#       paste0(
#         "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_",
#         i
#       )
#     dir.create(name,
#                recursive = TRUE)
#   })

setwd("data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess")

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

final_cluster_info <-
  final_cluster_info %>%
  dplyr::left_join(variable_info, by = "variable_id")

table(final_cluster_info$cluster, final_cluster_info$class)

library(org.Hs.eg.db)
library(clusterProfiler)

###proteomics
# for (i in 1:11) {
#   path <- paste0("cluster_", i)
#   dir.create(file.path(path, "proteomics_pathway"))
#   dir.create(file.path(path, "proteomics_pathway/GO_result/"))
#   dir.create(file.path(path, "proteomics_pathway/KEGG_result/"))
#   dir.create(file.path(path, "proteomics_pathway/Reactome_result/"))
#   message(i)
#   cluster_info <-
#     final_cluster_info %>%
#     dplyr::filter(cluster == i)
# 
#   ####proteomics
#   temp_data_proteomics <-
#     cluster_info %>%
#     dplyr::filter(class == "proteomics")
# 
#   if (nrow(temp_data_proteomics) == 0) {
#     next()
#   }
# 
#   ###GO enrichment
#   if (all(dir(file.path(path, "proteomics_pathway/GO_result")) != "proteomics_go")) {
#     proteomics_go <-
#       enrichGO(
#         gene = temp_data_proteomics$ENTREZID[!is.na(temp_data_proteomics$ENTREZID)],
#         OrgDb = org.Hs.eg.db,
#         keyType = "ENTREZID",
#         ont = "ALL",
#         pvalueCutoff = 0.05,
#         pAdjustMethod = "fdr",
#         qvalueCutoff = 0.05
#       )
# 
#     if (!is.null(proteomics_go)) {
#       proteomics_go@result <-
#         proteomics_go@result %>%
#         dplyr::filter(p.adjust < 0.05)
#     }
# 
#     save(proteomics_go,
#          file = file.path(path,
#                           "proteomics_pathway/GO_result/proteomics_go"))
#   } else{
#     load(file.path(path,
#                    "proteomics_pathway/GO_result/proteomics_go"))
#   }
# 
#   #enrichment KEGG analysis
#   if (all(dir(file.path(
#     path, "proteomics_pathway/KEGG_result"
#   )) != "proteomics_kegg")) {
#     proteomics_kegg <-
#       enrichKEGG(
#         gene = temp_data_proteomics$UNIPROT[!is.na(temp_data_proteomics$UNIPROT)],
#         organism = "hsa",
#         keyType = "uniprot",
#         pvalueCutoff = 0.05,
#         pAdjustMethod = "fdr",
#         qvalueCutoff = 0.05
#       )
# 
#     if (!is.null(proteomics_kegg)) {
#       proteomics_kegg@result =
#         proteomics_kegg@result %>%
#         dplyr::filter(p.adjust < 0.05)
# 
#       result <- proteomics_kegg@result
# 
#       kegg_class <-
#         purrr::map(
#           .x = result$ID,
#           .f = function(x) {
#             temp <- KEGGREST::keggGet(dbentries = x)
#             class <- temp[[1]]$CLASS
#             if (is.null(class)) {
#               class <- "no"
#             }
#             class
#           }
#         )  %>%
#         unlist()
# 
#       proteomics_kegg@result$kegg_class <- kegg_class
#     }
# 
#     save(
#       proteomics_kegg,
#       file = file.path(path,
#                        "proteomics_pathway/KEGG_result/proteomics_kegg")
#     )
#   } else{
#     load(file.path(path,
#                    "proteomics_pathway/KEGG_result/proteomics_kegg"))
#   }
# 
#   ##reactome pathway enrichment
#   library(ReactomePA)
#   if (all(dir(file.path(
#     path, "proteomics_pathway/Reactome_result"
#   )) != "proteomics_reactome")) {
#     proteomics_reactome <-
#       ReactomePA::enrichPathway(
#         gene = temp_data_proteomics$ENTREZID[!is.na(temp_data_proteomics$ENTREZID)],
#         organism = "human",
#         pvalueCutoff = 0.05,
#         pAdjustMethod = "fdr",
#         qvalueCutoff = 0.05
#       )
# 
#     if (!is.null(proteomics_reactome)) {
#       proteomics_reactome@result <-
#         proteomics_reactome@result %>%
#         dplyr::filter(p.adjust < 0.05)
# 
#       ##add disease class to Reactome pathways
#       disease <-
#         purrr::map(proteomics_reactome@result$ID, function(x) {
#           temp <- ReactomeContentService4R::getPathways(x)$isInDisease
#           if (is.null(temp)) {
#             return(FALSE)
#           }
#           temp
#         }) %>%
#         unlist()
# 
#       proteomics_reactome@result$disease <-
#         disease
# 
#     }
# 
#     save(
#       proteomics_reactome,
#       file = file.path(
#         path,
#         "proteomics_pathway/Reactome_result/proteomics_reactome"
#       )
#     )
#   } else{
#     load(file.path(
#       path,
#       "proteomics_pathway/Reactome_result/proteomics_reactome"
#     ))
#   }
# 
#   ##get the enriched results and only remain the pathways with P < 0.05 and Count >= 3
#   if (!is.null(proteomics_go)) {
#     proteomics_go <-
#       proteomics_go@result %>%
#       dplyr::filter(p.adjust < 0.05) %>%
#       dplyr::filter(Count >= 3) %>%
#       dplyr::arrange(p.adjust)
#   } else{
#     proteomics_go <-
#       data.frame(
#         ONTOLOGY = character(),
#         ID = character(),
#         Description = character(),
#         GeneRatio = numeric(),
#         BgRatio = numeric(),
#         pvalue = numeric(),
#         p.adjust = numeric(),
#         qvalue = numeric(),
#         geneID = character(),
#         Count = numeric()
#       )
#   }
# 
#   dim(proteomics_go)
# 
#   ###KEGG
#   if (!is.null(proteomics_kegg)) {
#     proteomics_kegg <-
#       proteomics_kegg@result %>%
#       dplyr::filter(p.adjust < 0.05) %>%
#       dplyr::arrange(p.adjust) %>%
#       dplyr::filter(Count >= 3)
#   } else{
#     proteomics_kegg <-
#       data.frame(
#         ID = character(),
#         Description = character(),
#         GeneRatio = numeric(),
#         BgRatio = numeric(),
#         pvalue = numeric(),
#         p.adjust = numeric(),
#         qvalue = numeric(),
#         geneID = character(),
#         Count = numeric()
#       )
#   }
# 
# 
#   ####enrichment reactome
# 
#   if (!is.null(proteomics_reactome)) {
#     proteomics_reactome <-
#       proteomics_reactome@result %>%
#       dplyr::filter(p.adjust < 0.05) %>%
#       # dplyr::filter(!disease) %>%
#       dplyr::filter(Count >= 3) %>%
#       dplyr::arrange(p.adjust)
#   } else{
#     proteomics_reactome <-
#       data.frame(
#         ID = character(),
#         Description = character(),
#         GeneRatio = numeric(),
#         BgRatio = numeric(),
#         pvalue = numeric(),
#         p.adjust = numeric(),
#         qvalue = numeric(),
#         geneID = character(),
#         Count = numeric()
#       )
#   }
# 
#   # #######get the similarity between GO terms
#   if (all(dir(file.path(path, "proteomics_pathway/GO_result")) != "go_sim_matrix")) {
#     go_sim_matrix =
#       get_go_result_sim(result = proteomics_go, sim.cutoff = 0)
# 
#     save(go_sim_matrix,
#          file = file.path(path, "proteomics_pathway/GO_result/go_sim_matrix"))
#   } else{
#     load(file.path(path, "proteomics_pathway/GO_result/go_sim_matrix"))
#   }
# 
#   go_sim_matrix =
#     go_sim_matrix %>%
#     dplyr::filter(sim > 0.7)
# 
#   edge_data <-
#     rbind(go_sim_matrix) %>%
#     dplyr::rename(from = name1, to = name2)
# 
#   node_data <-
#     rbind(proteomics_go) %>%
#     as.data.frame() %>%
#     dplyr::select(ID, everything()) %>%
#     dplyr::rename(node = ID)
# 
#   temp_graph <-
#     tidygraph::tbl_graph(nodes = node_data,
#                          edges = edge_data,
#                          directed = FALSE) %>%
#     dplyr::mutate(degree = tidygraph::centrality_degree())
# 
#   library(ggraph)
#   library(igraph)
# 
#   if (all(dir(file.path(path, "proteomics_pathway/GO_result")) != "subnetwork")) {
#     subnetwork <-
#       igraph::cluster_edge_betweenness(graph = temp_graph,
#                                        weights = abs(edge_attr(temp_graph,
#                                                                "sim")))
#     save(subnetwork,
#          file = file.path(path, "proteomics_pathway/GO_result/subnetwork"))
# 
#   } else{
#     load(file.path(path, "proteomics_pathway/GO_result/subnetwork"))
#   }
# 
# 
#   cluster <-
#     as.character(membership(subnetwork)) %>%
#     purrr::map(function(x) {
#       if (sum(x == as.character(membership(subnetwork))) == 1) {
#         return("Other")
#       } else{
#         return(x)
#       }
#     }) %>%
#     unlist()
# 
#   new_cluster <-
#     purrr::map(cluster, function(x) {
#       paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
#     }) %>%
#     unlist()
# 
#   new_cluster[new_cluster == "Module NA"] <- "Other"
# 
#   temp_graph <-
#     temp_graph %>%
#     tidygraph::mutate(module = new_cluster) %>%
#     activate(what = "nodes")
#   # dplyr::filter(module != "Other")
# 
#   ###clustered different GO terms
#   result_go <-
#     igraph::vertex_attr(temp_graph) %>%
#     do.call(cbind, .) %>%
#     as.data.frame() %>%
#     dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#     dplyr::arrange(ONTOLOGY, module, p.adjust)
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "cluster 1")
# 
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = result_go,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(
#     wb,
#     file = file.path(path, "proteomics_pathway/GO_result/result_go.xlsx"),
#     overwrite = TRUE
#   )
# 
#   ###plot to show the clusters of GO terms
#   cluster_label1 =
#     igraph::as_data_frame(temp_graph, what = "vertices") %>%
#     # dplyr::filter(module != "Other") %>%
#     dplyr::group_by(module) %>%
#     dplyr::filter(p.adjust == min(p.adjust)) %>%
#     pull(Description)
# 
#   cluster_label2 =
#     igraph::as_data_frame(temp_graph, what = "vertices") %>%
#     dplyr::filter(module == "Other") %>%
#     pull(Description)
# 
#   cluster_label = c(cluster_label1, cluster_label2)
# 
#   plot <-
#     temp_graph %>%
#     tidygraph::filter(module != "Other") %>%
#     ggraph(layout = 'kk',
#            circular = FALSE) +
#     geom_edge_link(
#       aes(width = sim),
#       strength = 1,
#       color = "black",
#       alpha = 1,
#       show.legend = FALSE
#     ) +
#     geom_node_point(
#       aes(fill = module,
#           size = -log(p.adjust, 10)),
#       shape = 21,
#       alpha = 0.7,
#       show.legend = FALSE
#     ) +
#     geom_node_text(aes(
#       x = x,
#       y = y,
#       label = ifelse(Description %in% cluster_label, Description, NA)
#     ),
#     size = 3,
#     repel = TRUE) +
#     guides(fill = guide_legend(ncol = 1)) +
#     scale_edge_width_continuous(range = c(0.1, 2)) +
#     scale_size_continuous(range = c(1, 7)) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "right",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     )
# 
#   plot
# 
#   extrafont::loadfonts()
# 
#   ggplot2::ggsave(
#     plot,
#     filename =
#       file.path(path, "proteomics_pathway/GO_result/GO_sim_plot.pdf"),
#     width = 7,
#     height = 7
#   )
# 
#   dir.create(file.path(path, "proteomics_pathway/GO_result/GO_sim_plot"))
# 
#   if (length(dir(
#     file.path(path,
#               "proteomics_pathway/GO_result/GO_sim_plot")
#   )) == 0) {
#     for (temp_cluster in unique(new_cluster)) {
#       cat(temp_cluster, " ")
#       if (temp_cluster == "Other") {
#         next()
#       }
#       plot1 <-
#         temp_graph %>%
#         tidygraph::filter(module == temp_cluster) %>%
#         ggraph(layout = 'kk',
#                circular = FALSE) +
#         geom_edge_link(
#           aes(width = sim),
#           strength = 1,
#           color = "black",
#           alpha = 1,
#           show.legend = TRUE
#         ) +
#         geom_node_point(
#           aes(fill = -log(p.adjust, 10),
#               size = Count),
#           shape = 21,
#           alpha = 1,
#           show.legend = TRUE
#         ) +
#         shadowtext::geom_shadowtext(
#           aes(
#             x = x,
#             y = y,
#             label = ifelse(module == "Other", NA, Description)
#           ),
#           check_overlap = TRUE,
#           size = 3,
#           color = "black",
#           bg.color = "white"
#         ) +
#         guides(fill = guide_legend(ncol = 1)) +
#         scale_edge_width_continuous(range = c(0.1, 2)) +
#         scale_size_continuous(range = c(3, 10)) +
#         ggraph::theme_graph() +
#         theme(
#           plot.background = element_rect(fill = "transparent", color = NA),
#           panel.background = element_rect(fill = "transparent", color = NA),
#           legend.position = "left",
#           legend.background = element_rect(fill = "transparent", color = NA)
#         )
# 
#       # plot1
# 
#       plot2 <-
#         result_go %>%
#         dplyr::filter(module == temp_cluster) %>%
#         dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#         dplyr::arrange(p.adjust) %>%
#         dplyr::mutate(Description = factor(Description, levels = Description)) %>%
#         ggplot(aes(p.adjust, Description)) +
#         geom_bar(stat = "identity", color = "grey") +
#         geom_text(
#           aes(x = 0, Description, label = Description),
#           hjust = 0,
#           size = 5,
#           color = "red"
#         ) +
#         theme_bw() +
#         labs(y = "", x = "-log10(FDR adjusted P value)") +
#         scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
#         theme(
#           panel.grid = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank()
#         )
# 
#       library(ggwordcloud)
#       library(plyr)
#       temp_data =
#         result_go %>%
#         dplyr::filter(module == temp_cluster) %>%
#         dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#         dplyr::select(Description, p.adjust) %>%
#         dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
#         plyr::dlply(.variables = .(Description)) %>%
#         purrr::map(function(x) {
#           data.frame(
#             word = stringr::str_split(x$Description, " ")[[1]],
#             p.adjust = x$p.adjust
#           )
#         }) %>%
#         do.call(rbind, .) %>%
#         as.data.frame() %>%
#         plyr::dlply(.variables = .(word)) %>%
#         purrr::map(function(x) {
#           x$p.adjust <- sum(x$p.adjust)
#           x %>%
#             dplyr::distinct(word, .keep_all = TRUE)
#         }) %>%
#         do.call(rbind, .) %>%
#         as.data.frame() %>%
#         dplyr::filter(!word %in% remove_words)
# 
#       plot3 =
#         temp_data %>%
#         ggplot(aes(label = word, size = p.adjust)) +
#         geom_text_wordcloud() +
#         scale_radius(range = c(5, 15), limits = c(0, NA)) +
#         theme_minimal()
# 
#       library(patchwork)
#       plot =
#         plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)
# 
#       ggsave(
#         plot,
#         filename = file.path(
#           path,
#           "proteomics_pathway/GO_result/GO_sim_plot",
#           paste(temp_cluster, "sim_plot.pdf", sep = "_")
#         ),
#         width = 21,
#         height = 7
#       )
#     }
#   }
# 
#   ##matrix tow show the cluster GO terms
#   result_go %>%
#     dplyr::group_by(ONTOLOGY) %>%
#     dplyr::summarise(n = n()) %>%
#     dplyr::mutate(n = n * 10 / max(n) + 2)
# 
#   library(igraph)
#   library(simplifyEnrichment)
# 
#   if (length(dir(file.path(path, "proteomics_pathway/GO_result"),
#                  pattern = "sim_matrix.pdf")) == 0) {
#     for (ont in c('MF', "BP", "CC")) {
#       cat(ont, " ")
#       show_matrix_cluster(
#         result = result_go %>% dplyr::mutate(Direction = "UP") %>% dplyr::rename(cluster = module),
#         ont = ont,
#         measure = "Wang",
#         remove_words = remove_words,
#         margin = 15,
#         width = 16,
#         height = 14,
#         path = file.path(path, "proteomics_pathway/GO_result"),
#         top = 15
#       )
#     }
#   }
# 
#   #####output GO result as clusters
#   library(plyr)
#   if (nrow(result_go) > 0) {
#     result_go_cluster <-
#       result_go %>%
#       plyr::dlply(.variables = .(module)) %>%
#       purrr::map(function(x) {
#         if (nrow(x) == 1) {
#           return(x)
#         }
# 
#         if (x$module[1] == 'Other') {
#           return(x)
#         }
# 
#         x =
#           x %>%
#           dplyr::arrange(p.adjust)
# 
#         x$node <-
#           paste(x$node, collapse = ";")
# 
#         x$Description <-
#           paste(x$Description, collapse = ";")
# 
#         x$BgRatio <-
#           paste(x$BgRatio, collapse = ";")
# 
#         x$pvalue <- min(as.numeric(x$pvalue))
#         x$p.adjust <- min(as.numeric(x$p.adjust))
#         x$qvalue <- min(as.numeric(x$qvalue))
#         x$geneID =
#           x$geneID %>%
#           stringr::str_split(pattern = "/") %>%
#           unlist() %>%
#           unique() %>%
#           paste(collapse = '/')
# 
#         x$Count <-
#           length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
# 
#         x =
#           x %>%
#           dplyr::select(module, everything()) %>%
#           dplyr::distinct(module, .keep_all = TRUE)
# 
#         x
# 
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(
#         module_annotation = case_when(
#           module == "Other" ~ Description,
#           module != "Other" ~ stringr::str_split(Description, ";")[[1]][1]
#         )
#       ) %>%
#       dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#       dplyr::arrange(p.adjust) %>%
#       dplyr::select(module_annotation, everything())
# 
#     result_go_cluster$module_annotation =
#       stringr::str_split(result_go_cluster$Description, ";") %>%
#       purrr::map(function(x) {
#         x[1]
#       }) %>%
#       unlist()
#   } else{
#     result_go_cluster <-
#       proteomics_go
# 
#     result_go_cluster <-
#       result_go_cluster %>%
#       dplyr::mutate(module_annotation = Description,
#                     database = "GO") %>%
#       dplyr::rename(pathway_id = ID) %>%
#       dplyr::select(module_annotation,
#                     Description,
#                     p.adjust,
#                     Count,
#                     database,
#                     geneID,
#                     pathway_id)
#   }
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "cluster 1")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = result_go_cluster,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(
#     wb,
#     file = file.path(
#       path,
#       "proteomics_pathway/GO_result/result_go_cluster.xlsx"
#     ),
#     overwrite = TRUE
#   )
# 
#   ###output the cluster annotation for each cluster
#   dir.create(file.path(path, "proteomics_pathway/GO_result/GO_graph"))
#   unique(result_go_cluster$module) %>%
#     purrr::map(
#       .f = function(x) {
#         cat(x, " ")
#         if (x == "Other") {
#           return(NULL)
#         }
# 
#         temp_id =
#           result_go_cluster %>%
#           dplyr::filter(module == x) %>%
#           dplyr::pull(node) %>%
#           stringr::str_split(";") %>%
#           `[[`(1) %>%
#           pRoloc::goIdToTerm(keepNA = FALSE) %>%
#           data.frame(id = ., class = "YES") %>%
#           tibble::rownames_to_column(var = "name")
# 
#         temp_plot =
#           GOSim::getGOGraph(term = temp_id$name, prune = Inf) %>%
#           igraph::igraph.from.graphNEL() %>%
#           tidygraph::as_tbl_graph() %>%
#           tidygraph::left_join(temp_id, by = "name") %>%
#           dplyr::mutate(class = case_when(is.na(class) ~ "NO",
#                                           TRUE ~ class))
# 
#         plot =
#           temp_plot %>%
#           ggraph(layout = 'kk',
#                  circular = FALSE) +
#           geom_edge_link(
#             color = ggsci::pal_aaas()(n = 10)[1],
#             alpha = 1,
#             arrow = grid::arrow(
#               angle = 10,
#               length = unit(0.2, "inches"),
#               type = "closed"
#             ),
#             show.legend = FALSE
#           ) +
#           geom_node_point(
#             aes(fill = class),
#             shape = 21,
#             alpha = 1,
#             size = 6,
#             show.legend = FALSE
#           ) +
#           geom_node_text(aes(
#             x = x,
#             y = y,
#             label = ifelse(class == "YES", id, NA)
#           ),
#           size = 3,
#           repel = TRUE) +
#           scale_fill_manual(values = c('YES' = "red", 'NO' = "white")) +
#           ggraph::theme_graph() +
#           theme(
#             plot.background = element_rect(fill = "transparent", color = NA),
#             panel.background = element_rect(fill = "transparent", color = NA),
#             legend.position = "left",
#             legend.background = element_rect(fill = "transparent", color = NA)
#           )
#         plot
#         ggsave(
#           plot,
#           filename = file.path(
#             path,
#             "proteomics_pathway/GO_result/GO_graph",
#             paste(x, "_GO graph.pdf", sep = "")
#           ),
#           width = 7,
#           height = 7
#         )
#       }
#     )
# 
#   #######-----------------------------------------------------------------------
#   #### KEGG
#   ###calculate the similarity of KEGG pathway
#   message('KEGG.....')
#   kegg_sim_matrix <-
#     tryCatch(
#       simplifyEnrichment::term_similarity_from_KEGG(term_id = c(proteomics_kegg$ID),
#                                                     method = "jaccard") %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column(var = "name1") %>%
#         tidyr::pivot_longer(
#           cols = -name1,
#           names_to = "name2",
#           values_to = "sim"
#         ) %>%
#         dplyr::filter(name1 != name2),
#       error = function(e) {
#         data.frame(name1 = "x",
#                    name2 = "y",
#                    sim = 0) %>%
#           dplyr::filter(sim > 9)
#       }
#     )
# 
#   save(kegg_sim_matrix,
#        file = file.path(path, "proteomics_pathway/KEGG_result/kegg_sim_matrix"))
# 
#   kegg_sim_matrix =
#     kegg_sim_matrix %>%
#     dplyr::filter(sim > 0.5)
# 
#   edge_data <-
#     kegg_sim_matrix %>%
#     dplyr::rename(from = name1, to = name2)
# 
#   node_data <-
#     rbind(proteomics_kegg) %>%
#     as.data.frame() %>%
#     dplyr::select(ID, everything()) %>%
#     dplyr::rename(node = ID) %>%
#     dplyr::filter(node %in% c(edge_data$from, edge_data$to))
# 
#   temp_graph <-
#     tidygraph::tbl_graph(nodes = node_data,
#                          edges = edge_data,
#                          directed = FALSE) %>%
#     dplyr::mutate(degree = tidygraph::centrality_degree())
# 
#   library(ggraph)
# 
#   subnetwork <-
#     igraph::cluster_edge_betweenness(graph = temp_graph,
#                                      weights = abs(edge_attr(temp_graph,
#                                                              "sim")))
# 
#   library(igraph)
# 
#   cluster <-
#     as.character(membership(subnetwork)) %>%
#     purrr::map(function(x) {
#       if (sum(x == as.character(membership(subnetwork))) == 1) {
#         return("Other")
#       } else{
#         return(x)
#       }
#     }) %>%
#     unlist()
# 
#   new_cluster <-
#     purrr::map(cluster, function(x) {
#       paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
#     }) %>%
#     unlist()
# 
#   new_cluster[new_cluster == "Module NA"] <- "Other"
# 
#   temp_graph <-
#     temp_graph %>%
#     tidygraph::mutate(module = new_cluster)
# 
#   ####clustered different KEGG terms
#   result_kegg <-
#     igraph::vertex_attr(temp_graph) %>%
#     do.call(cbind, .) %>%
#     as.data.frame() %>%
#     dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#     dplyr::arrange(p.adjust)
# 
#   save(result_kegg,
#        file = file.path(path, "proteomics_pathway/KEGG_result/result_kegg"))
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "cluster 1")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = result_kegg,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(wb,
#                file = file.path(path, "KEGG_result/result_kegg.xlsx"),
#                overwrite = TRUE)
# 
#   plot <-
#     temp_graph %>%
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
#       label = ifelse(cluster == "Other", Description, Description)
#     ),
#     size = 3,
#     repel = TRUE) +
#     guides(fill = guide_legend(ncol = 1)) +
#     # scale_fill_manual(values = c(
#     #   "Other" = "grey"
#     # )) +
#     scale_edge_width_continuous(range = c(0.1, 2)) +
#     scale_size_continuous(range = c(1, 7)) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "right",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     )
#   # facet_nodes(Direction~ONTOLOGY)
# 
#   plot
# 
#   ggsave(
#     plot,
#     filename = file.path(path,
#                          "proteomics_pathway/KEGG_result/kegg_sim_plot.pdf"),
#     width = 9,
#     height = 7
#   )
# 
#   dir.create(file.path(path, "proteomics_pathway/KEGG_result/KEGG_sim_plot"))
# 
#   for (temp_cluster in unique(new_cluster)) {
#     cat(temp_cluster, " ")
#     plot1 <-
#       temp_graph %>%
#       tidygraph::filter(module == temp_cluster) %>%
#       ggraph(layout = 'fr',
#              circular = FALSE) +
#       geom_edge_link(
#         aes(width = sim),
#         strength = 1,
#         color = "black",
#         alpha = 1,
#         show.legend = TRUE
#       ) +
#       geom_node_point(
#         aes(fill = cluster,
#             size = -log(p.adjust, 10)),
#         shape = 21,
#         alpha = 1,
#         show.legend = TRUE
#       ) +
#       geom_node_text(aes(
#         x = x,
#         y = y,
#         label = ifelse(module == "Other", NA, Description)
#       ),
#       size = 3,
#       repel = TRUE) +
#       guides(fill = guide_legend(ncol = 1)) +
#       scale_edge_width_continuous(range = c(0.1, 2)) +
#       scale_size_continuous(range = c(1, 7)) +
#       ggraph::theme_graph() +
#       theme(
#         plot.background = element_rect(fill = "transparent", color = NA),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         legend.position = "left",
#         legend.background = element_rect(fill = "transparent", color = NA)
#       )
# 
#     plot1
# 
#     plot2 <-
#       result_kegg %>%
#       dplyr::filter(module == temp_cluster) %>%
#       dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#       dplyr::arrange(p.adjust) %>%
#       dplyr::mutate(Description = factor(Description, levels = Description)) %>%
#       ggplot(aes(p.adjust, Description)) +
#       geom_bar(stat = "identity") +
#       geom_text(aes(x = 0, Description, label = Description),
#                 hjust = 0,
#                 size = 5) +
#       theme_bw() +
#       labs(y = "", x = "-log10(FDR adjusted P value)") +
#       scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
#       theme(
#         panel.grid = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#       )
# 
#     library(ggwordcloud)
# 
#     temp_data =
#       result_kegg %>%
#       dplyr::filter(module == temp_cluster) %>%
#       dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#       dplyr::select(Description, p.adjust) %>%
#       dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
#       plyr::dlply(.variables = .(Description)) %>%
#       purrr::map(function(x) {
#         data.frame(
#           word = stringr::str_split(x$Description, " ")[[1]],
#           p.adjust = x$p.adjust
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       plyr::dlply(.variables = .(word)) %>%
#       purrr::map(function(x) {
#         x$p.adjust <- sum(x$p.adjust)
#         x %>%
#           dplyr::distinct(word, .keep_all = TRUE)
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::filter(!word %in% remove_words)
# 
#     plot3 =
#       temp_data %>%
#       ggplot(aes(label = word, size = p.adjust)) +
#       geom_text_wordcloud() +
#       scale_radius(range = c(5, 15), limits = c(0, NA)) +
#       theme_minimal()
# 
#     library(patchwork)
#     plot =
#       plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)
# 
#     ggsave(
#       plot,
#       filename = file.path(
#         path,
#         "proteomics_pathway/KEGG_result/KEGG_sim_plot",
#         paste(temp_cluster, "sim_plot.pdf", sep = "_")
#       ),
#       width = 21,
#       height = 7
#     )
#   }
# 
# 
#   #matrix tow show the cluster KEGG terms
#   show_matrix_cluster(
#     result = result_kegg %>% dplyr::mutate(Direction = "UP") %>% dplyr::rename(cluster = module),
#     ont = NULL,
#     measure = 'jaccard',
#     remove_words = remove_words,
#     margin = 15,
#     height = 8,
#     path = file.path(path, "KEGG_result")
#   )
# 
#   #####output KEGG result
#   if (nrow(result_kegg) > 0) {
#     result_kegg_cluster <-
#       result_kegg %>%
#       plyr::dlply(.variables = .(module)) %>%
#       purrr::map(function(x) {
#         if (nrow(x) == 1) {
#           return(x)
#         }
# 
#         if (x$cluster[1] == 'Other') {
#           return(x)
#         }
# 
#         x =
#           x %>%
#           dplyr::arrange(p.adjust)
# 
#         x$node <-
#           paste(x$node, collapse = ";")
# 
#         x$Description <-
#           paste(x$Description, collapse = ";")
# 
#         x$BgRatio <-
#           paste(x$BgRatio, collapse = ";")
# 
#         x$pvalue <- min(as.numeric(x$pvalue))
#         x$p.adjust <- min(as.numeric(x$p.adjust))
#         x$qvalue <- min(as.numeric(x$qvalue))
#         x$geneID =
#           x$geneID %>%
#           stringr::str_split(pattern = "/") %>%
#           unlist() %>%
#           unique() %>%
#           paste(collapse = '/')
# 
#         x$Count <-
#           length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
# 
#         x =
#           x %>%
#           dplyr::select(module, everything()) %>%
#           dplyr::distinct(module, .keep_all = TRUE)
# 
#         x
# 
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(module_annotation = case_when(module == "Other" ~ Description,
#                                                   module != "Other" ~ cluster)) %>%
#       dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#       dplyr::arrange(p.adjust) %>%
#       dplyr::select(module_annotation, everything())
# 
#     result_kegg_cluster$module_annotation =
#       stringr::str_split(result_kegg_cluster$Description, ";") %>%
#       purrr::map(function(x) {
#         x[1]
#       }) %>%
#       unlist()
#   } else{
#     result_kegg_cluster <-
#       proteomics_kegg
# 
#     result_kegg_cluster <-
#       result_kegg_cluster %>%
#       dplyr::mutate(module_annotation = Description,
#                     database = "KEGG") %>%
#       dplyr::rename(pathway_id = ID) %>%
#       dplyr::select(module_annotation,
#                     Description,
#                     p.adjust,
#                     Count,
#                     database,
#                     geneID,
#                     pathway_id)
# 
#   }
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "cluster 1")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = result_kegg_cluster,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(
#     wb,
#     file = file.path(
#       path,
#       "proteomics_pathway/KEGG_result/result_kegg_cluster.xlsx"
#     ),
#     overwrite = TRUE
#   )
# 
# 
#   #########-----------------------------------------------------------------------
#   message("Reactome")
#   ####Reactome
#   ##calculate the similarity of Reactome pathway
#   reactome_sim_matrix <-
#     tryCatch(
#       expr =  simplifyEnrichment::term_similarity_from_Reactome(
#         term_id = c(proteomics_reactome$ID),
#         method = "jaccard"
#       ) %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column(var = "name1") %>%
#         tidyr::pivot_longer(
#           cols = -name1,
#           names_to = "name2",
#           values_to = "sim"
#         ) %>%
#         dplyr::filter(name1 != name2),
#       error = function(e) {
#         data.frame(name1 = character(),
#                    name2 = character(),
#                    sim = numeric())
#       }
#     )
# 
#   save(
#     reactome_sim_matrix,
#     file = file.path(
#       path,
#       "proteomics_pathway/Reactome_result/reactome_sim_matrix"
#     )
#   )
# 
#   reactome_sim_matrix =
#     reactome_sim_matrix %>%
#     dplyr::filter(sim > 0.5)
# 
#   edge_data <-
#     reactome_sim_matrix %>%
#     dplyr::rename(from = name1, to = name2)
# 
#   node_data <-
#     rbind(proteomics_reactome) %>%
#     as.data.frame() %>%
#     dplyr::select(ID, everything()) %>%
#     dplyr::rename(node = ID)
# 
#   temp_graph <-
#     tidygraph::tbl_graph(nodes = node_data,
#                          edges = edge_data,
#                          directed = FALSE) %>%
#     dplyr::mutate(degree = tidygraph::centrality_degree())
# 
#   library(ggraph)
# 
#   subnetwork <-
#     igraph::cluster_edge_betweenness(graph = temp_graph,
#                                      weights = abs(edge_attr(temp_graph,
#                                                              "sim")))
#   cluster <-
#     as.character(membership(subnetwork)) %>%
#     purrr::map(function(x) {
#       if (sum(x == as.character(membership(subnetwork))) == 1) {
#         return("Other")
#       } else{
#         return(x)
#       }
#     }) %>%
#     unlist()
# 
#   new_cluster <-
#     purrr::map(cluster, function(x) {
#       paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
#     }) %>%
#     unlist()
# 
#   new_cluster[new_cluster == "Module NA"] <- "Other"
# 
#   temp_graph <-
#     temp_graph %>%
#     tidygraph::mutate(module = new_cluster) %>%
#     activate(what = "nodes")
#   # dplyr::filter(module != "Other")
# 
#   ###cluster different Reactome terms
#   result_reactome <-
#     igraph::vertex_attr(temp_graph) %>%
#     do.call(cbind, .) %>%
#     as.data.frame() %>%
#     dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#     dplyr::arrange(p.adjust)
# 
#   save(
#     result_reactome,
#     file = file.path(
#       path,
#       "proteomics_pathway/Reactome_result/result_reactome"
#     )
#   )
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "cluster 1")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = result_reactome,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(
#     wb,
#     file.path(
#       path,
#       "proteomics_pathway/Reactome_result/result_reactome.xlsx"
#     ),
#     overwrite = TRUE
#   )
# 
#   cluster_label1 =
#     igraph::as_data_frame(temp_graph, what = "vertices") %>%
#     dplyr::filter(module != "Other") %>%
#     dplyr::group_by(module) %>%
#     dplyr::filter(p.adjust == min(p.adjust)) %>%
#     pull(Description)
# 
#   cluster_label2 =
#     igraph::as_data_frame(temp_graph, what = "vertices") %>%
#     dplyr::filter(module == "Other") %>%
#     pull(Description)
# 
#   cluster_label = c(cluster_label1, cluster_label2)
# 
#   plot <-
#     temp_graph %>%
#     tidygraph::filter(module != "Other") %>%
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
#       aes(fill = module,
#           size = -log(p.adjust, 10)),
#       shape = 21,
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     shadowtext::geom_shadowtext(
#       aes(
#         x = x,
#         y = y,
#         label = ifelse(
#           Description %in% cluster_label,
#           stringr::str_replace(Description, "Homo sapiens\\\r: ", ""),
#           NA
#         )
#       ),
#       check_overlap = TRUE,
#       size = 3,
#       color = "black",
#       bg.color = "white"
#     ) +
#     guides(fill = guide_legend(ncol = 1)) +
#     scale_edge_width_continuous(range = c(0.1, 2)) +
#     scale_size_continuous(range = c(1, 7)) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "right",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     )
#   # facet_nodes(Direction~ONTOLOGY)
# 
#   plot
# 
#   ggsave(
#     plot,
#     filename =
#       file.path(
#         path,
#         "proteomics_pathway/Reactome_result/Reactome_sim_plot.pdf"
#       ),
#     width = 8.5,
#     height = 7
#   )
# 
#   dir.create(file.path(
#     path,
#     "proteomics_pathway/Reactome_result/Reactome_sim_plot"
#   ))
# 
#   for (temp_cluster in unique(new_cluster)) {
#     cat(temp_cluster, " ")
#     if (temp_cluster == "Other") {
#       next()
#     }
#     plot1 <-
#       temp_graph %>%
#       tidygraph::filter(module == temp_cluster) %>%
#       ggraph(layout = 'fr',
#              circular = FALSE) +
#       geom_edge_link(
#         aes(width = sim),
#         strength = 1,
#         color = "black",
#         alpha = 1,
#         show.legend = TRUE
#       ) +
#       geom_node_point(
#         aes(fill = -log(p.adjust, 10),
#             size = Count),
#         shape = 21,
#         alpha = 1,
#         show.legend = TRUE
#       ) +
#       shadowtext::geom_shadowtext(
#         aes(
#           x = x,
#           y = y,
#           label = ifelse(module == "Other", NA, Description)
#         ),
#         check_overlap = TRUE,
#         size = 3,
#         color = "black",
#         bg.color = "white"
#       ) +
#       guides(fill = guide_legend(ncol = 1)) +
#       scale_edge_width_continuous(range = c(0.1, 2)) +
#       scale_size_continuous(range = c(3, 10)) +
#       ggraph::theme_graph() +
#       theme(
#         plot.background = element_rect(fill = "transparent", color = NA),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         legend.position = "left",
#         legend.background = element_rect(fill = "transparent", color = NA)
#       )
# 
#     plot1
# 
#     plot2 <-
#       result_reactome %>%
#       dplyr::filter(module == temp_cluster) %>%
#       dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#       dplyr::arrange(p.adjust) %>%
#       dplyr::mutate(Description = factor(Description, levels = Description)) %>%
#       ggplot(aes(p.adjust, Description)) +
#       geom_bar(stat = "identity", color = "grey") +
#       geom_text(
#         aes(x = 0, Description, label = Description),
#         hjust = 0,
#         size = 5,
#         color = "red"
#       ) +
#       theme_bw() +
#       labs(y = "", x = "-log10(FDR adjusted P value)") +
#       scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
#       theme(
#         panel.grid = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#       )
# 
#     library(ggwordcloud)
# 
#     temp_data =
#       result_reactome %>%
#       dplyr::filter(module == temp_cluster) %>%
#       dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
#       dplyr::select(Description, p.adjust) %>%
#       dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
#       dplyr::mutate(Description = stringr::str_replace_all(Description, "\\(", "")) %>%
#       dplyr::mutate(Description = stringr::str_replace_all(Description, "\\)", "")) %>%
#       dplyr::mutate(Description = stringr::str_replace_all(Description, "\\:", "")) %>%
#       plyr::dlply(.variables = .(Description)) %>%
#       purrr::map(function(x) {
#         data.frame(
#           word = stringr::str_split(x$Description, " ")[[1]],
#           p.adjust = x$p.adjust
#         )
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       plyr::dlply(.variables = .(word)) %>%
#       purrr::map(function(x) {
#         x$p.adjust <- sum(x$p.adjust)
#         x %>%
#           dplyr::distinct(word, .keep_all = TRUE)
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::filter(!word %in% remove_words)
# 
#     plot3 =
#       temp_data %>%
#       ggplot(aes(label = word, size = p.adjust)) +
#       geom_text_wordcloud() +
#       scale_radius(range = c(5, 15), limits = c(0, NA)) +
#       theme_minimal()
# 
#     library(patchwork)
#     plot =
#       plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)
# 
#     ggsave(
#       plot,
#       filename = file.path(
#         path,
#         "proteomics_pathway/Reactome_result/Reactome_sim_plot",
#         paste(temp_cluster, "sim_plot.pdf", sep = "_")
#       ),
#       width = 21,
#       height = 7
#     )
#   }
# 
#   #matrix tow show the cluster Reactome terms
#   library(simplifyEnrichment)
# 
#   show_matrix_cluster(
#     result = result_reactome %>%
#       dplyr::mutate(Direction = "UP") %>%
#       dplyr::rename(cluster = module),
#     ont = NULL,
#     measure = "jaccard",
#     remove_words = remove_words,
#     margin = 15,
#     width = 16,
#     height = 12,
#     path = file.path(path, "proteomics_pathway/Reactome_result")
#   )
# 
#   #####output Reactome result
#   if (nrow(result_reactome) > 0) {
#     result_reactome_cluster <-
#       result_reactome %>%
#       plyr::dlply(.variables = .(module)) %>%
#       purrr::map(function(x) {
#         if (nrow(x) == 1) {
#           return(x)
#         }
# 
#         if (x$module[1] == 'Other') {
#           return(x)
#         }
# 
#         x =
#           x %>%
#           dplyr::arrange(p.adjust)
# 
#         x$node <-
#           paste(x$node, collapse = ";")
# 
#         x$Description <-
#           paste(x$Description, collapse = ";")
# 
#         x$BgRatio <-
#           paste(x$BgRatio, collapse = ";")
# 
#         x$pvalue <- min(as.numeric(x$pvalue))
#         x$p.adjust <- min(as.numeric(x$p.adjust))
#         x$qvalue <- min(as.numeric(x$qvalue))
#         x$geneID =
#           x$geneID %>%
#           stringr::str_split(pattern = "/") %>%
#           unlist() %>%
#           unique() %>%
#           paste(collapse = '/')
# 
#         x$Count <-
#           length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
# 
#         x =
#           x %>%
#           dplyr::select(module, everything()) %>%
#           dplyr::distinct(module, .keep_all = TRUE)
# 
#         x
# 
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(module_annotation = case_when(module == "Other" ~ Description,
#                                                   module != "Other" ~ module)) %>%
#       dplyr::arrange(p.adjust) %>%
#       dplyr::select(module_annotation, everything())
# 
#     result_reactome_cluster$module_annotation =
#       stringr::str_split(result_reactome_cluster$Description, ";") %>%
#       purrr::map(function(x) {
#         x[1]
#       }) %>%
#       unlist()
#   } else{
#     result_reactome_cluster <-
#       proteomics_reactome
# 
#     result_reactome_cluster <-
#       result_reactome_cluster %>%
#       dplyr::mutate(module_annotation = Description,
#                     database = "KEGG") %>%
#       dplyr::rename(pathway_id = ID) %>%
#       dplyr::select(module_annotation,
#                     Description,
#                     p.adjust,
#                     Count,
#                     database,
#                     geneID,
#                     pathway_id)
#   }
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "cluster 1")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = result_reactome_cluster,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(
#     wb,
#     file =
#       file.path(
#         path,
#         "proteomics_pathway/Reactome_result/result_reactome_cluster.xlsx"
#       ),
#     overwrite = TRUE
#   )
# 
#   ####---------------------------------------------------------------------------
#   ########combine three databases together
#   ####barplot
#   cat("stop here")
# 
#   if (nrow(result_go) > 0) {
#     result_go_cluster <-
#       rbind(
#         readxl::read_xlsx(
#           file.path(
#             path,
#             "proteomics_pathway/GO_result/result_go_cluster.xlsx"
#           ),
#           sheet = 1
#         ) %>%
#           dplyr::filter(ONTOLOGY != "CC") %>%
#           dplyr::arrange(p.adjust) %>%
#           dplyr::mutate(database = "GO") %>%
#           dplyr::select(
#             module_annotation,
#             Description,
#             p.adjust,
#             Count,
#             database,
#             geneID,
#             pathway_id = node
#           )
#       ) %>%
#       dplyr::mutate(Count = as.numeric(Count)) %>%
#       dplyr::filter(!is.na(module_annotation))
#   }
#   #
#   # library(openxlsx)
#   # wb <- createWorkbook()
#   # modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   # addWorksheet(wb, sheetName = "cluster 1")
#   # freezePane(
#   #   wb = wb,
#   #   sheet = 1,
#   #   firstRow = TRUE,
#   #   firstCol = TRUE
#   # )
#   #
#   # proteomics_kegg
# 
#   if (nrow(result_kegg) > 0) {
#     result_kegg_cluster <-
#       rbind(
#         readxl::read_xlsx(
#           file.path(
#             path,
#             "proteomics_pathway/KEGG_result/result_kegg_cluster.xlsx"
#           ),
#           sheet = 1
#         ) %>%
#           dplyr::arrange(p.adjust) %>%
#           dplyr::mutate(database = "KEGG") %>%
#           dplyr::select(
#             module_annotation,
#             Description,
#             p.adjust,
#             Count,
#             database,
#             geneID,
#             pathway_id = node
#           )
#       ) %>%
#       dplyr::mutate(Count = as.numeric(Count)) %>%
#       dplyr::filter(!is.na(module_annotation))
#   }
# 
#   if (nrow(result_reactome) > 0) {
#     result_reactome_cluster <-
#       rbind(
#         readxl::read_xlsx(
#           file.path(
#             path,
#             "proteomics_pathway/Reactome_result/result_reactome_cluster.xlsx"
#           ),
#           sheet = 1
#         ) %>%
#           dplyr::arrange(p.adjust) %>%
#           dplyr::mutate(database = "Reactome") %>%
#           dplyr::select(
#             module_annotation,
#             Description,
#             p.adjust,
#             Count,
#             database,
#             geneID,
#             pathway_id = node
#           )
#       ) %>%
#       dplyr::mutate(Count = as.numeric(Count)) %>%
#       dplyr::filter(!is.na(module_annotation)) %>%
#       dplyr::mutate(module_annotation = stringr::str_replace(module_annotation, "Homo sapiens\\\r: ", ""))
# 
#   }
# 
# 
#   dim(result_go_cluster)
#   dim(result_kegg_cluster)
#   dim(result_reactome_cluster)
# 
#   if (nrow(result_go_cluster) == 0 &
#       nrow(result_kegg_cluster) == 0 &
#       nrow(result_reactome_cluster) == 0) {
#     next
#   }
# 
#   ######calculate the similarity (jaccard index) between all the pathways
#   jaccard_index <-
#     tryCatch(
#       expr = get_jaccard_index_for_three_databases(
#         result_go_cluster = result_go_cluster,
#         result_kegg_cluster = result_kegg_cluster,
#         result_reactome_cluster = result_reactome_cluster,
#         variable_info = variable_info
#       ),
#       error = function(e) {
#         data.frame(name1 = character(),
#                    name2 = character(),
#                    value = numeric())
#       }
#     )
# 
#   head(jaccard_index)
# 
#   edge_data <-
#     jaccard_index %>%
#     dplyr::filter(value > 0.5) %>%
#     dplyr::rename(from = name1, to = name2, sim = value)
# 
# 
# 
#   if (nrow(result_go_cluster) == 0) {
#     temp_data1 <- NULL
#   } else{
#     temp_data1 <-
#       data.frame(result_go_cluster, class = "GO")
#   }
# 
#   if (nrow(result_kegg_cluster) == 0) {
#     temp_data2 <- NULL
#   } else{
#     temp_data2 <-
#       data.frame(result_kegg_cluster, class = "KEGG")
#   }
# 
#   if (nrow(result_reactome_cluster) == 0) {
#     temp_data3 <- NULL
#   } else{
#     temp_data3 <-
#       data.frame(result_reactome_cluster, class = "Reactome")
#   }
# 
#   node_data <-
#     rbind(temp_data1,
#           temp_data2,
#           temp_data3)
#   # dplyr::filter(module_annotation %in% c(edge_data$from, edge_data$to))
# 
#   edge_data <-
#     edge_data %>%
#     dplyr::filter(from %in% node_data$module_annotation &
#                     to %in% node_data$module_annotation)
# 
#   temp_graph <-
#     tidygraph::tbl_graph(nodes = node_data,
#                          edges = edge_data,
#                          directed = FALSE) %>%
#     dplyr::mutate(degree = tidygraph::centrality_degree())
# 
#   library(ggraph)
# 
#   subnetwork <-
#     igraph::cluster_edge_betweenness(graph = temp_graph,
#                                      weights = abs(edge_attr(temp_graph,
#                                                              "sim")))
#   cluster <-
#     as.character(membership(subnetwork)) %>%
#     purrr::map(function(x) {
#       if (sum(x == as.character(membership(subnetwork))) == 1) {
#         return("Other")
#       } else{
#         return(x)
#       }
#     }) %>%
#     unlist()
# 
#   new_cluster <-
#     purrr::map(cluster, function(x) {
#       paste("Module", match(x, unique(cluster)[unique(cluster) != "Other"]))
#     }) %>%
#     unlist()
# 
#   new_cluster[new_cluster == "Module NA"] <- "Other"
# 
#   temp_graph <-
#     temp_graph %>%
#     tidygraph::mutate(module = new_cluster)
# 
#   ###cluster different Reactome terms
#   result_all <-
#     igraph::vertex_attr(temp_graph) %>%
#     do.call(cbind, .) %>%
#     as.data.frame() %>%
#     dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
#     dplyr::arrange(module, p.adjust)
# 
#   ##result_all is the pahtway enrichment for all three databases
#   save(result_all,
#        file = file.path(path, "proteomics_pathway/result_all"))
# 
#   ###add new information to
#   result_all$module[result_all$module == "Other"] =
#     paste("Other", 1:sum(result_all$module == "Other"))
# 
#   library(plyr)
#   result_all =
#     result_all %>%
#     plyr::dlply(.variables = .(module)) %>%
#     purrr::map(function(x) {
#       if (nrow(x) == 1) {
#         return(x %>% dplyr::select(-degree))
#       } else{
#         x =
#           x %>%
#           dplyr::arrange(p.adjust)
#         x$module_annotation = x$module_annotation[1]
#         x$Description = paste(x$Description, collapse = ";")
#         x$p.adjust = x$p.adjust[1]
#         x$database = paste(x$database, collapse = ";")
#         x$geneID =
#           x$geneID %>%
#           stringr::str_split("/") %>%
#           unlist() %>%
#           unique() %>%
#           paste(collapse = "/")
#         x$Count = length(stringr::str_split(x$geneID, pattern = "/")[[1]])
#         x$pathway_id = paste(x$pathway_id, collapse = ";")
#         x$class = paste(x$class, collapse = ";")
#         x$degree = paste(x$degree, collapse = ";")
#         x$module = x$module[1]
#         x %>%
#           dplyr::distinct(module_annotation, .keep_all = TRUE) %>%
#           dplyr::select(-degree)
#       }
#     }) %>%
#     do.call(rbind, .) %>%
#     as.data.frame()
# 
#   library(openxlsx)
#   wb <- createWorkbook()
#   modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")
#   addWorksheet(wb, sheetName = "cluster 1")
#   freezePane(
#     wb = wb,
#     sheet = 1,
#     firstRow = TRUE,
#     firstCol = TRUE
#   )
# 
#   writeDataTable(
#     wb,
#     sheet = 1,
#     x = result_all,
#     colNames = TRUE,
#     rowNames = FALSE
#   )
# 
#   saveWorkbook(wb,
#                file.path(path, "proteomics_pathway/result_all.xlsx"),
#                overwrite = TRUE)
# 
#   cluster_label1 <-
#     igraph::as_data_frame(temp_graph, what = "vertices") %>%
#     dplyr::filter(module != "Other") %>%
#     dplyr::group_by(module) %>%
#     dplyr::filter(p.adjust == min(p.adjust)) %>%
#     pull(module_annotation)
# 
#   cluster_label2 <-
#     igraph::as_data_frame(temp_graph, what = "vertices") %>%
#     dplyr::filter(module == "Other") %>%
#     pull(module_annotation)
# 
#   cluster_label = c(cluster_label1, cluster_label2)
# 
#   plot <-
#     temp_graph %>%
#     activate(what = "nodes") %>%
#     tidygraph::filter(module != "Other") %>%
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
#       aes(fill = class,
#           size = -log(p.adjust, 10)),
#       shape = 21,
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     shadowtext::geom_shadowtext(
#       aes(
#         x = x,
#         y = y,
#         label = ifelse(module_annotation %in% cluster_label1, module_annotation, NA)
#       ),
#       check_overlap = TRUE,
#       size = 3,
#       color = "black",
#       bg.color = "white"
#     ) +
#     guides(fill = guide_legend(ncol = 1)) +
#     scale_edge_width_continuous(range = c(0.1, 2)) +
#     scale_size_continuous(range = c(1, 10)) +
#     scale_fill_manual(values = database_color) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "bottom",
#       legend.background = element_rect(fill = "transparent", color = NA),
#     )
# 
#   plot
# 
#   extrafont::loadfonts(device = "pdf")
# 
#   ggsave(
#     plot,
#     filename = file.path(path,
#                          "proteomics_pathway/GO_KEGG_Reactome_sim_plot.pdf"),
#     width = 6,
#     height = 7
#   )
#   # result_all
# }


###cluster 2
