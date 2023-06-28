# no_source()
rm(list = ls())
gc()
library(tidyverse)
library(tidymass)
library(plyr)
setwd(masstools::get_project_wd())

source("code/tools.R")

dir.create("data_analysis/combined_omics/DE_SWAN")
dir.create(
  "data_analysis/combined_omics/DE_SWAN/transcriptome_proteome_metabolome_summary"
)
setwd("data_analysis/combined_omics/DE_SWAN")

#####Proteomics
load("proteomics_pathway/result_all")

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

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

result_all_crest1 <- result_all

####crest2
load("proteomics_pathway_crest2/result_all")

###add new information to
result_all$module[result_all$module == "Other"] =
  paste("Other", 1:sum(result_all$module == "Other"))

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

result_all_crest2 <- result_all

result_all_crest1 <-
  result_all_crest1 %>%
  dplyr::arrange(p.adjust)

result_all_crest2 <-
  result_all_crest2 %>%
  dplyr::arrange(p.adjust)

temp_crest1 <-
  result_all_crest1 %>%
  head(20)

temp_crest2 <-
  result_all_crest2 %>%
  head(20)

temp_crest1_proteomics <-
  temp_crest1

temp_crest2_proteomics <-
  temp_crest2

library(clusterProfiler)

load("proteomics_pathway/GO_result/proteomics_crest1_go")
load("proteomics_pathway/KEGG_result/proteomics_crest1_kegg")
load("proteomics_pathway/Reactome_result/proteomics_crest1_reactome")

load("proteomics_pathway_crest2/GO_result/proteomics_crest2_go")
load("proteomics_pathway_crest2/KEGG_result/proteomics_crest2_kegg")
load("proteomics_pathway_crest2/Reactome_result/proteomics_crest2_reactome")

#####create a network
####tanscriptomics
temp_crest1_proteomics$module_annotation

###rename
###crest1
proteomics_node3_go <-
  proteomics_crest1_go@result %>%
  dplyr::select(ID,
                ONTOLOGY,
                Description,
                p.adjust,
                geneID,
                Count) %>%
  dplyr::rename(node_id = ID) %>%
  dplyr::mutate(level = 3) %>%
  dplyr::mutate(
    pathway_id = node_id,
    node_id = paste("proteomics_pathway_crest1", node_id, sep = "_"),
    contents = node_id,
    module_annotation = Description
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node3_kegg <-
  proteomics_crest1_kegg@result %>%
  dplyr::select(ID,
                Description,
                p.adjust,
                geneID,
                Count) %>%
  dplyr::rename(node_id = ID) %>%
  dplyr::mutate(level = 3) %>%
  dplyr::mutate(
    pathway_id = node_id,
    node_id = paste("proteomics_pathway_crest1", node_id, sep = "_"),
    contents = node_id,
    module_annotation = Description,
    ONTOLOGY = NA
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node3_reactome <-
  proteomics_crest1_reactome@result %>%
  dplyr::select(ID,
                Description,
                p.adjust,
                geneID,
                Count) %>%
  dplyr::rename(node_id = ID) %>%
  dplyr::mutate(level = 3) %>%
  dplyr::mutate(
    pathway_id = node_id,
    node_id = paste("proteomics_pathway_crest1", node_id, sep = "_"),
    contents = node_id,
    module_annotation = Description,
    ONTOLOGY = NA
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node3 <-
  rbind(proteomics_node3_go,
        proteomics_node3_kegg,
        proteomics_node3_reactome)

proteomics_node2 <-
  temp_crest1_proteomics %>%
  dplyr::select(
    module_annotation,
    p.adjust,
    Count,
    database,
    class,
    module,
    Description,
    pathway_id,
    geneID
  ) %>%
  dplyr::mutate(contents = pathway_id) %>%
  dplyr::rename(node_id = module) %>%
  dplyr::mutate(level = 2) %>%
  dplyr::mutate(
    node_id = paste("proteomics_module_crest1", node_id, sep = "_"),
    ONTOLOGY = NA
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node1 <-
  data.frame(
    node_id = c(
      "CVD_crest1",
      "Skin-Muscle_crest1",
      "lipid_metabolism_crest1",
      "Alcohol_metabolism_crest1"
    ),
    contents = c(
      "Complement and coagulation cascades;Platelet degranulation;blood coagulation;fibrinolysis;HDL remodeling;Complement cascade;protease binding;protein activation cascade",
      "extracellular matrix structural constituent;Scavenging of heme from plasma;glycosaminoglycan binding;cytolysis",
      "Plasma lipoprotein remodeling;Chylomicron assembly",
      "alcohol binding"
    )
  )

match("alcohol binding",
      proteomics_node2$module_annotation)

#####remove some node
proteomics_node2 <-
  proteomics_node2 %>%
  dplyr::filter(module_annotation %in% unlist(stringr::str_split(proteomics_node1$contents, ";")))

proteomics_node3 <-
  proteomics_node3 %>%
  dplyr::filter(pathway_id %in% unlist(stringr::str_split(proteomics_node2$pathway_id, ";")))

###get edge data
###pathway to protein
proteomics_edge_crest1_3 <-
  purrr::map2(proteomics_node3$node_id,
              proteomics_node3$geneID,
              function(x, y) {
                y <- stringr::str_split(y, "\\/")[[1]]
                temp <- data.frame(from = x,
                                   to = y)
                return(temp)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::distinct(.keep_all = TRUE)

proteomics_edge_crest1_3$to <-
  paste("proteomics", proteomics_edge_crest1_3$to, sep = "_")

####module to pathway
proteomics_edge_crest1_2 <-
  purrr::map2(proteomics_node2$node_id,
              proteomics_node2$contents,
              function(x, y) {
                y <- stringr::str_split(y, ";")[[1]]
                temp <- data.frame(from = x,
                                   to = y)
                return(temp)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::distinct(.keep_all = TRUE)

proteomics_edge_crest1_2 <-
  proteomics_edge_crest1_2 %>%
  dplyr::left_join(proteomics_node3[, c("pathway_id", "node_id")],
                   by = c("to" = "pathway_id")) %>%
  dplyr::select(from, node_id) %>%
  dplyr::rename(to = node_id)

proteomics_edge_crest1_1 <-
  purrr::map2(proteomics_node1$node_id,
              proteomics_node1$contents,
              function(x, y) {
                y <- stringr::str_split(y, ";")[[1]]
                temp <- data.frame(from = x,
                                   to = y)
                return(temp)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::distinct(.keep_all = TRUE)

proteomics_edge_crest1_1 <-
  proteomics_edge_crest1_1 %>%
  dplyr::left_join(proteomics_node2[, c("node_id", "module_annotation")],
                   by = c("to" = "module_annotation")) %>%
  dplyr::select(from, node_id) %>%
  dplyr::rename(to = node_id)


proteomics_edge_crest1 <-
  rbind(proteomics_edge_crest1_1,
        proteomics_edge_crest1_2,
        proteomics_edge_crest1_3)

proteomics_node_crest1_1 <-
  data.frame(
    node_id = proteomics_node1$node_id,
    class = "function",
    Count = stringr::str_split(proteomics_node1$contents, ";") %>%
      lapply(length) %>%
      unlist()
  ) %>%
  dplyr::mutate(
    ONTOLOGY = NA,
    p.adjust = NA,
    module_annotation = NA,
    level = 1
  ) %>%
  dplyr::select(node_id,
                ONTOLOGY,
                p.adjust,
                module_annotation,
                Count,
                level,
                class) %>%
  dplyr::mutate(crest = 1)

proteomics_node_crest1_2 <-
  data.frame(
    node_id = proteomics_node2$node_id,
    ONTOLOGY = proteomics_node2$ONTOLOGY,
    p.adjust = proteomics_node2$p.adjust,
    module_annotation = proteomics_node2$module_annotation,
    Count = as.numeric(proteomics_node2$Count),
    level = proteomics_node2$level,
    class = "module"
  ) %>%
  dplyr::select(node_id,
                ONTOLOGY,
                p.adjust,
                module_annotation,
                Count,
                level,
                class) %>%
  dplyr::mutate(crest = 1)

proteomics_node_crest1_3 <-
  rbind(
    data.frame(
      node_id = proteomics_node3$node_id,
      ONTOLOGY = proteomics_node3$ONTOLOGY,
      p.adjust = proteomics_node3$p.adjust,
      module_annotation = proteomics_node3$module_annotation,
      Count = as.numeric(proteomics_node3$Count),
      level = proteomics_node3$level,
      class = "pathway"
    ),
    data.frame(
      node_id = unique(unlist(
        stringr::str_split(proteomics_node3$geneID, "\\/")
      )),
      ONTOLOGY = NA,
      p.adjust = NA,
      module_annotation = NA,
      Count = NA,
      level = 4,
      class = "protein"
    ) %>%
      dplyr::mutate(node_id = paste("proteomics", node_id, sep = "_"))
  ) %>%
  dplyr::select(node_id,
                ONTOLOGY,
                p.adjust,
                module_annotation,
                Count,
                level,
                class) %>%
  dplyr::mutate(crest = 1)

proteomics_node_crest1 <-
  rbind(proteomics_node_crest1_1,
        proteomics_node_crest1_2,
        proteomics_node_crest1_3)

library(ggraph)
library(tidygraph)
library(igraph)

proteomics_node_crest1$node_id %in% c(proteomics_edge_crest1$from, proteomics_edge_crest1$to)
c(proteomics_edge_crest1$from, proteomics_edge_crest1$to) %in% proteomics_node_crest1$node_id

save(proteomics_node_crest1, file = "transcriptome_proteome_metabolome_summary/proteomics_node_crest1")
save(proteomics_edge_crest1, file = "transcriptome_proteome_metabolome_summary/proteomics_edge_crest1")

###crest2
proteomics_node3_go <-
  proteomics_crest2_go@result %>%
  dplyr::select(ID,
                ONTOLOGY,
                Description,
                p.adjust,
                geneID,
                Count) %>%
  dplyr::rename(node_id = ID) %>%
  dplyr::mutate(level = 3) %>%
  dplyr::mutate(
    pathway_id = node_id,
    node_id = paste("proteomics_pathway_crest2", node_id, sep = "_"),
    contents = node_id,
    module_annotation = Description
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node3_kegg <-
  proteomics_crest2_kegg@result %>%
  dplyr::select(ID,
                Description,
                p.adjust,
                geneID,
                Count) %>%
  dplyr::rename(node_id = ID) %>%
  dplyr::mutate(level = 3) %>%
  dplyr::mutate(
    pathway_id = node_id,
    node_id = paste("proteomics_pathway_crest2", node_id, sep = "_"),
    contents = node_id,
    module_annotation = Description,
    ONTOLOGY = NA
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node3_reactome <-
  proteomics_crest2_reactome@result %>%
  dplyr::select(ID,
                Description,
                p.adjust,
                geneID,
                Count) %>%
  dplyr::rename(node_id = ID) %>%
  dplyr::mutate(level = 3) %>%
  dplyr::mutate(
    pathway_id = node_id,
    node_id = paste("proteomics_pathway_crest2", node_id, sep = "_"),
    contents = node_id,
    module_annotation = Description,
    ONTOLOGY = NA
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node3 <-
  rbind(proteomics_node3_go,
        proteomics_node3_kegg,
        proteomics_node3_reactome)

proteomics_node2 <-
  temp_crest2_proteomics %>%
  dplyr::select(
    module_annotation,
    p.adjust,
    Count,
    database,
    class,
    module,
    Description,
    pathway_id,
    geneID
  ) %>%
  dplyr::mutate(contents = pathway_id) %>%
  dplyr::rename(node_id = module) %>%
  dplyr::mutate(level = 2) %>%
  dplyr::mutate(
    node_id = paste("proteomics_module_crest2", node_id, sep = "_"),
    ONTOLOGY = NA
  ) %>%
  dplyr::select(
    node_id,
    ONTOLOGY,
    p.adjust,
    module_annotation,
    Description,
    Count,
    level,
    contents,
    geneID,
    pathway_id
  )

proteomics_node1 <-
  data.frame(
    node_id = c("CVD_crest2",
                "Skin-Muscle_crest2",
                "Kindey_function_crest2",
                "Carbohydrate_metabolism_crest2",
                "Immune_crest2"),
    contents = c(
      "Complement and coagulation cascades;Platelet degranulation;blood coagulation;fibrinolysis;HDL remodeling",
      "extracellular matrix structural constituent;Scavenging of heme from plasma;structural constituent of muscle;Regulation of actin cytoskeleton",
      "glomerular filtration",
      "carbohydrate binding",
      "acute-phase response;zymogen activation;antimicrobial humoral response;complement binding"
    )
  )

match("carbohydrate binding",
      proteomics_node2$module_annotation)

#####remove some node
proteomics_node2 <-
  proteomics_node2 %>%
  dplyr::filter(module_annotation %in% unlist(stringr::str_split(proteomics_node1$contents, ";")))

proteomics_node3 <-
  proteomics_node3 %>%
  dplyr::filter(pathway_id %in% unlist(stringr::str_split(proteomics_node2$pathway_id, ";")))

###get edge data
###pathway to protein
proteomics_edge_crest2_3 <-
  purrr::map2(proteomics_node3$node_id,
              proteomics_node3$geneID,
              function(x, y) {
                y <- stringr::str_split(y, "\\/")[[1]]
                temp <- data.frame(from = x,
                                   to = y)
                return(temp)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::distinct(.keep_all = TRUE)

proteomics_edge_crest2_3$to <-
  paste("proteomics", proteomics_edge_crest2_3$to, sep = "_")

####module to pathway
proteomics_edge_crest2_2 <-
  purrr::map2(proteomics_node2$node_id,
              proteomics_node2$contents,
              function(x, y) {
                y <- stringr::str_split(y, ";")[[1]]
                temp <- data.frame(from = x,
                                   to = y)
                return(temp)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::distinct(.keep_all = TRUE)

proteomics_edge_crest2_2 <-
  proteomics_edge_crest2_2 %>%
  dplyr::left_join(proteomics_node3[, c("pathway_id", "node_id")],
                   by = c("to" = "pathway_id")) %>%
  dplyr::select(from, node_id) %>%
  dplyr::rename(to = node_id)

proteomics_edge_crest2_1 <-
  purrr::map2(proteomics_node1$node_id,
              proteomics_node1$contents,
              function(x, y) {
                y <- stringr::str_split(y, ";")[[1]]
                temp <- data.frame(from = x,
                                   to = y)
                return(temp)
              }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::distinct(.keep_all = TRUE)

proteomics_edge_crest2_1 <-
  proteomics_edge_crest2_1 %>%
  dplyr::left_join(proteomics_node2[, c("node_id", "module_annotation")],
                   by = c("to" = "module_annotation")) %>%
  dplyr::select(from, node_id) %>%
  dplyr::rename(to = node_id)


proteomics_edge_crest2 <-
  rbind(proteomics_edge_crest2_1,
        proteomics_edge_crest2_2,
        proteomics_edge_crest2_3)

proteomics_node_crest2_1 <-
  data.frame(
    node_id = proteomics_node1$node_id,
    class = "function",
    Count = stringr::str_split(proteomics_node1$contents, ";") %>%
      lapply(length) %>%
      unlist()
  ) %>%
  dplyr::mutate(
    ONTOLOGY = NA,
    p.adjust = NA,
    module_annotation = NA,
    level = 1
  ) %>%
  dplyr::select(node_id,
                ONTOLOGY,
                p.adjust,
                module_annotation,
                Count,
                level,
                class) %>%
  dplyr::mutate(crest = 2)

proteomics_node_crest2_2 <-
  data.frame(
    node_id = proteomics_node2$node_id,
    ONTOLOGY = proteomics_node2$ONTOLOGY,
    p.adjust = proteomics_node2$p.adjust,
    module_annotation = proteomics_node2$module_annotation,
    Count = as.numeric(proteomics_node2$Count),
    level = proteomics_node2$level,
    class = "module"
  ) %>%
  dplyr::select(node_id,
                ONTOLOGY,
                p.adjust,
                module_annotation,
                Count,
                level,
                class) %>%
  dplyr::mutate(crest = 2)

proteomics_node_crest2_3 <-
  rbind(
    data.frame(
      node_id = proteomics_node3$node_id,
      ONTOLOGY = proteomics_node3$ONTOLOGY,
      p.adjust = proteomics_node3$p.adjust,
      module_annotation = proteomics_node3$module_annotation,
      Count = as.numeric(proteomics_node3$Count),
      level = proteomics_node3$level,
      class = "pathway"
    ),
    data.frame(
      node_id = unique(unlist(
        stringr::str_split(proteomics_node3$geneID, "\\/")
      )),
      ONTOLOGY = NA,
      p.adjust = NA,
      module_annotation = NA,
      Count = NA,
      level = 4,
      class = "protein"
    ) %>%
      dplyr::mutate(node_id = paste("proteomics", node_id, sep = "_"))
  ) %>%
  dplyr::select(node_id,
                ONTOLOGY,
                p.adjust,
                module_annotation,
                Count,
                level,
                class) %>%
  dplyr::mutate(crest = 2)

proteomics_node_crest2 <-
  rbind(proteomics_node_crest2_1,
        proteomics_node_crest2_2,
        proteomics_node_crest2_3)

library(ggraph)
library(tidygraph)
library(igraph)

proteomics_node_crest2$node_id %in% c(proteomics_edge_crest2$from, proteomics_edge_crest2$to)
c(proteomics_edge_crest2$from, proteomics_edge_crest1$to) %in% proteomics_node_crest1$node_id

save(proteomics_node_crest2, file = "transcriptome_proteome_metabolome_summary/proteomics_node_crest2")
save(proteomics_edge_crest2, file = "transcriptome_proteome_metabolome_summary/proteomics_edge_crest2")

#####combine crest1 and crest2
proteomics_node <-
  rbind(proteomics_node_crest1,
        proteomics_node_crest2) %>%
  dplyr::distinct(node_id, .keep_all = TRUE) %>%
  dplyr::mutate(label = module_annotation)

proteomics_node$label[proteomics_node$class == "function"] <-
  proteomics_node$node_id[proteomics_node$class == "function"] %>%
  stringr::str_replace("_crest[0-1]{1,2}", "")

proteomics_node$label[proteomics_node$class == "pathway"] <-
  NA

proteomics_node$Count[proteomics_node$class == "protein"] <- 1

proteomics_node$Count[proteomics_node$class == "function"] <- 300

proteomics_edge <-
  rbind(proteomics_edge_crest1,
        proteomics_edge_crest2)

save(proteomics_node, file = "transcriptome_proteome_metabolome_summary/proteomics_node")
save(proteomics_edge, file = "transcriptome_proteome_metabolome_summary/proteomics_edge")

graph_data <-
  tidygraph::tbl_graph(nodes = proteomics_node,
                       edges = proteomics_edge,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

g <- graph_data

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite")
# dplyr::select(x, y)

coords$y[coords$crest == 1 & coords$class == "function"] <- 0
coords$y[coords$crest == 1 & coords$class == "module"] <- 1
coords$y[coords$crest == 1 & coords$class == "pathway"] <- 2
coords$y[coords$class == "protein"] <- 3
coords$y[coords$crest == 2 & coords$class == "pathway"] <- 4
coords$y[coords$crest == 2 & coords$class == "module"] <- 5
coords$y[coords$crest == 2 & coords$class == "function"] <- 6

range(coords$x[coords$y == 3])

##function
coords$x[coords$y == 0] <-
  new_coords(range_left = 250,
             range_right = 750,
             old_x = coords$x[coords$y == 0])
coords$x[coords$y == 6] <-
  new_coords(range_left = 250,
             range_right = 750,
             old_x = coords$x[coords$y == 6])

##module
coords$x[coords$y == 1] <-
  new_coords(range_left = 100,
             range_right = 900,
             old_x = coords$x[coords$y == 1])
coords$x[coords$y == 5] <-
  new_coords(range_left = 100,
             range_right = 900,
             old_x = coords$x[coords$y == 5])

##pathway
coords$x[coords$y == 2] <-
  new_coords(range_left = 50,
             range_right = 950,
             old_x = coords$x[coords$y == 2])
coords$x[coords$y == 4] <-
  new_coords(range_left = 50,
             range_right = 950,
             old_x = coords$x[coords$y == 4])

###protein
coords$x[coords$y == 3] <-
  new_coords(range_left = 0,
             range_right = 1000,
             old_x = coords$x[coords$y == 3])

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
  )

# plot <-
ggraph(my_graph,
       layout = 'bipartite') +
  geom_edge_diagonal(color = "grey",
                     show.legend = TRUE) +
  geom_node_point(shape = 21,
                  aes(fill = class,
                      size = Count),
                  show.legend = TRUE) +
  scale_fill_manual(values = network_class_color) +
  geom_node_text(
    aes(x = x,
        y = y,
        label = label),
    # size = 3,
    alpha = 1,
    show.legend = FALSE,
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ) +
  scale_size_continuous(range = c(1, 5)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom"
  )

