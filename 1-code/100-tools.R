metabolite_super_class <-
  c(
    "Acetylides",
    "Alkaloids and derivatives",
    "Allenes",
    "Benzenoids",
    "Homogeneous metal compounds",
    "Homogeneous non-metal compounds",
    "Hydrocarbon derivatives",
    "Hydrocarbons",
    "Lignans, neolignans and related compounds",
    "Lipids and lipid-like molecules",
    "Miscellaneous inorganic compounds",
    "Mixed metal/non-metal compounds",
    "Nucleosides, nucleotides, and analogues",
    "Organic 1,3-dipolar compounds",
    "Organic acids and derivatives",
    "Organic compounds",
    "Organic nitrogen compounds",
    "Organic oxygen compounds",
    "Organic Polymers",
    "Organic salts",
    "Organohalogen compounds",
    "Organoheterocyclic compounds",
    "Organometallic compounds",
    "Organophosphorus compounds",
    "Phenylpropanoids and polyketides"
  )

metabolite_super_class_color <-
  colorRampPalette(colors = ggsci::pal_lancet()(n = 9)[1:9])(n = length(metabolite_super_class))

names(metabolite_super_class_color) <- metabolite_super_class

lipid_class <- c(
  "Fatty Acyls",
  "Glycerolipids",
  "Glycerophospholipids",
  "Prenol lipids",
  "Sphingolipids",
  "Steroids and steroid derivatives"
)

lipid_class_color <- RColorBrewer::brewer.pal(n = 6, name = "Dark2")
names(lipid_class_color) <- lipid_class

sex_color <-
  c("F" = ggsci::pal_aaas()(n = 5)[2],
    "M" = ggsci::pal_aaas()(n = 5)[1])

library(ggplot2)

theme_base <-
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank()
  )


sample_class_color <-
  c(
    "Blank" = "#6F99ADFF",
    "QC" = "#BC3C29FF",
    "Subject" = "#E18727FF"
  )

marker_color <-
  c(
    "Up" = ggsci::pal_futurama()(n = 9)[6],
    "Down" = ggsci::pal_futurama()(n = 9)[3],
    "No" = ggsci::pal_futurama()(n = 9)[9]
  )


young_color <-
  c(young = ggsci::pal_tron()(n = 9)[2],
    old = ggsci::pal_material()(n = 9)[9])

# from <- "108F"
# to <- "26B"
#
# get_family_score(graph = graph,
#                  from = "72X",
#                  to = "DBL7")

get_family_score <-
  function(graph, from, to) {
    idx <- match(c(from, to), igraph::vertex_attr(graph = graph)$node)
    paths <-
      igraph::all_simple_paths(
        graph = graph,
        from = idx[1],
        to = idx[2],
        mode = "all"
      )
    
    if (length(paths) == 0) {
      return(Inf)
    }
    
    paths %>%
      purrr::map(function(x) {
        x <- as.numeric(x)
        if (length(x) == 2) {
          return(1)
        } else{
          break_point <-
            lapply(
              2:(length(x) - 1),
              FUN = function(z) {
                get_node_is_children_from_two_nodes(
                  graph = graph,
                  child_node_idx = x[z],
                  parient_node_idx1 = x[z -
                                          1],
                  parient_node_idx2 = x[z +
                                          1]
                )
              }
            ) %>%
            unlist()
          if (any(break_point)) {
            return(Inf)
          } else{
            return(length(x) - 1)
          }
          
        }
      }) %>%
      unlist() %>%
      min()
  }

get_node_is_children_from_two_nodes <-
  function(graph,
           child_node_idx,
           parient_node_idx1,
           parient_node_idx2) {
    edge1 <-
      igraph::shortest_paths(
        graph = graph,
        from = parient_node_idx1,
        to = child_node_idx,
        mode = "out"
      )$vpath[[1]]
    
    edge2 <-
      igraph::shortest_paths(
        graph = graph,
        from = parient_node_idx2,
        to = child_node_idx,
        mode = "out"
      )$vpath[[1]]
    
    if (length(edge1) != 0 & length(edge2) != 0) {
      return(TRUE)
    } else{
      return(FALSE)
    }
    
  }


family_distance_color <-
  c(
    "1" = ggsci::pal_gsea()(n = 10)[10],
    "2" = ggsci::pal_gsea()(n = 10)[7],
    "3" = ggsci::pal_gsea()(n = 10)[6],
    "4" = ggsci::pal_gsea()(n = 10)[1]
  )

convert_mass_dataset2pmd <-
  function(object) {
    data <-
      as.matrix(massdataset::extract_expression_data(object = object))
    group <- massdataset::extract_sample_info(object) %>%
      dplyr::select(sample_id, group) %>%
      dplyr::rename(sample_name = sample_id,
                    sample_group = group) %>%
      as.data.frame()
    list(
      data = data,
      group = group,
      mz = as.numeric(massdataset::extract_variable_info(object)$mz),
      rt = as.numeric(massdataset::extract_variable_info(object)$rt)
    )
  }




polarity_color <-
  c(
    "positive" = "#FC4E07",
    "negative" = "#00AFBB",
    "mix" = "#E7B800"
  )



optimize_loess_span <-
  function(x, y, span_range = seq(0.2, 0.6, 0.1)) {
    span_rmse <-
      purrr::map(span_range, function(span) {
        # cat(span, " ")
        temp_data =
          data.frame(x, y)
        
        prediction <-
          purrr::map(
            2:(nrow(temp_data) - 1),
            .f = function(idx) {
              temp_result =
                loess(formula = y ~ x,
                      data = temp_data[-idx,],
                      span = span)
              prediction =
                try(predict(object = temp_result,
                            newdata = temp_data[idx, -2, drop = FALSE]))
              
              if (class(prediction) == "try-error") {
                data.frame(real = temp_data$y[idx],
                           prediction = NA)
              } else{
                data.frame(real = temp_data$y[idx],
                           prediction = as.numeric(prediction))
              }
            }
          ) %>%
          dplyr::bind_rows() %>%
          dplyr::filter(!is.na(prediction))
        
        if (all(is.na(prediction$prediction))) {
          temp_rmse = NA
        } else{
          temp_rmse = sqrt(sum((
            prediction$real - prediction$prediction
          ) ^ 2) / nrow(prediction))
        }
        
        data.frame(span = span, rmse = temp_rmse)
      }) %>%
      dplyr::bind_rows()
    
    # span_rmse
    plot <-
      data.frame(x, y) %>%
      ggplot(aes(x, y)) +
      geom_point(size = 5) +
      # geom_line() +
      base_theme
    
    span_rmse =
      span_rmse %>%
      dplyr::filter(!is.na(rmse))
    idx = which.min(span_rmse$rmse)
    # for(i in 1:nrow(span_rmse)){
    plot =
      plot +
      geom_smooth(
        se = FALSE,
        span = span_rmse$span[idx],
        color = ggsci::pal_lancet()(n = 9)[idx]
      )
    # }
    
    plot =
      plot +
      ggplot2::ggtitle(label = paste("Span: ", span_rmse$span[idx])) +
      theme(title = element_text(colour = ggsci::pal_lancet()(n = 9)[idx]))
    
    list(span_rmse, plot)
  }


base_theme =
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    strip.text = element_text(size = 12)
  )




get_mfuzz_center <- function(data,
                             c,
                             membership_cutoff = 0.5) {
  data <-
    data@assayData$exprs
  
  membership <- c$membership
  
  membership <-
    membership %>%
    as.data.frame() %>%
    purrr::map(function(x) {
      rownames(membership)[which(x >= membership_cutoff)]
    })
  
  centers <-
    membership %>%
    purrr::map(function(x) {
      apply(data[x, , drop = FALSE], 2, mean)
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  centers
  
}


lm_adjusted_cor <-
  function(data_set1,
           data_set2,
           sample_info,
           method = c("spearman", "pearson", "all"),
           threads = 5) {
    method = match.arg(method)
    library(future)
    library(furrr)
    # plan(strategy = multisession(workers = threads))
    
    cor =
      rownames(data_set1) %>%
      purrr::map(function(name) {
        cat(name, " ")
        x = as.numeric(data_set1[name,])
        
        temp_cor =
          purrr::map(
            as.data.frame(t(data_set2)),
            .f = function(y) {
              temp_data =
                data.frame(x = x,
                           y = y,
                           sample_info)
              temp_data <-
                temp_data %>%
                dplyr::filter(!is.na(y))
              temp_data$Gender[temp_data$Gender == 'F'] = 0
              temp_data$Gender[temp_data$Gender == 'M'] = 1
              temp_data$Gender = as.numeric(temp_data$Gender)
              
              temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
              temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
              temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
              temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
              temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
              
              # temp_data$SSPG = as.numeric(temp_data$SSPG)
              # temp_data$FPG = as.numeric(temp_data$FPG)
              
              ##only remain the subjects with less than 5 time points
              remain_subject_id =
                temp_data %>%
                dplyr::select(subject_id) %>%
                dplyr::group_by(subject_id) %>%
                dplyr::summarise(n = n()) %>%
                dplyr::filter(n >= 5) %>%
                pull(subject_id)
              
              temp_data =
                temp_data %>%
                dplyr::filter(subject_id %in% remain_subject_id)
              
              ##linear mixed model to adjusted the x and y
              adjusted_x =
                lme4::lmer(
                  formula = x ~ Gender + Adj.age + Ethnicity + (1 |
                                                                  subject_id),
                  data = temp_data
                ) %>%
                residuals()
              
              adjusted_y =
                lme4::lmer(
                  formula = y ~ Gender + Adj.age + Ethnicity + (1 |
                                                                  subject_id),
                  data = temp_data
                ) %>%
                residuals()
              
              if (method == "all") {
                cor_value1 =
                  cor.test(adjusted_x, adjusted_y, method = "spearman")
                
                cor_value2 =
                  cor.test(adjusted_x, adjusted_y, method = "pearson")
                
                result1 =
                  c(cor = unname(cor_value1$estimate),
                    p = unname(cor_value1$p.value))
                
                if (is.na(result1[1])) {
                  result1[1] = 0
                }
                
                if (is.na(result1[2])) {
                  result1[2] = 1
                }
                
                result2 =
                  c(cor = unname(cor_value2$estimate),
                    p = unname(cor_value2$p.value))
                
                if (is.na(result2[1])) {
                  result2[1] = 0
                }
                
                if (is.na(result2[2])) {
                  result2[2] = 1
                }
                
                list(result1 = result1, result2 = result2)
                
              } else {
                cor_value =
                  cor.test(adjusted_x, adjusted_y, method = method)
                
                result =
                  c(cor = unname(cor_value$estimate),
                    p = unname(cor_value$p.value))
                
                if (is.na(result[1])) {
                  result[1] = 0
                }
                
                if (is.na(result[2])) {
                  result[2] = 1
                }
                list(result)
              }
            }
          )
        # do.call(rbind, .) %>%
        # as.data.frame()
        
        if (method == "all") {
          temp_cor1 =
            temp_cor %>%
            purrr::map(function(x) {
              x[[1]]
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "metabolite") %>%
            dplyr::mutate(microbiome = name) %>%
            dplyr::select(microbiome, metabolite, dplyr::everything())
          
          temp_cor1$p_adjust = p.adjust(temp_cor1$p, method = "BH")
          
          temp_cor2 =
            temp_cor %>%
            purrr::map(function(x) {
              x[[2]]
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "metabolite") %>%
            dplyr::mutate(microbiome = name) %>%
            dplyr::select(microbiome, metabolite, dplyr::everything())
          
          temp_cor2$p_adjust = p.adjust(temp_cor2$p, method = "BH")
          
          list(temp_cor1 = temp_cor1,
               temp_cor2 = temp_cor2)
        } else{
          temp_cor =
            temp_cor %>%
            tibble::rownames_to_column(var = "metabolite") %>%
            dplyr::mutate(microbiome = name) %>%
            dplyr::select(microbiome, metabolite, dplyr::everything())
          temp_cor$p_adjust = p.adjust(temp_cor$p, method = "BH")
          list(temp_cor)
        }
      })
    # do.call(rbind, .) %>%
    # as.data.frame()
    
    if (method == "all") {
      cor1 =
        cor %>%
        purrr::map(function(x) {
          x[[1]]
        }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      cor2 =
        cor %>%
        purrr::map(function(x) {
          x[[2]]
        }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      cor1$p_adjust2 =
        p.adjust(cor1$p, method = "BH")
      
      cor2$p_adjust2 =
        p.adjust(cor2$p, method = "BH")
      list(spearman = cor1,
           pearson = cor2)
    } else{
      cor$p_adjust2 =
        p.adjust(cor$p, method = "BH")
      list(cor)
    }
  }


lm_adjust <-
  function(expression_data,
           sample_info,
           threads = 5) {
    library(future)
    library(furrr)
    # plan(strategy = multisession(workers = threads))
    new_expression_data <-
      rownames(expression_data) %>%
      purrr::map(function(name) {
        # cat(name, " ")
        x = as.numeric(expression_data[name,])
        temp_data =
          data.frame(x = x,
                     sample_info)
        
        
        temp_data$Gender[temp_data$Gender == 'Female'] = 0
        temp_data$Gender[temp_data$Gender == 'Male'] = 1
        temp_data$Gender = as.numeric(temp_data$Gender)
        
        temp_data$IRIS[temp_data$IRIS == 'IR'] = 1
        temp_data$IRIS[temp_data$IRIS == 'IS'] = 2
        temp_data$IRIS[temp_data$IRIS == 'Unknown'] = 0
        temp_data$IRIS = as.numeric(temp_data$IRIS)
        
        adjusted_x <-
          lm(x ~ Gender + BMI + IRIS, data = temp_data) %>%
          residuals()
        adjusted_x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(new_expression_data) <-
      colnames(expression_data)
    
    rownames(new_expression_data) <-
      rownames(expression_data)
    new_expression_data
  }


library(wesanderson)
wesanderson::wes_palette(name = "Zissou1",
                         n = 5,
                         type = "discrete")


sex_color <-
  c(
    "Female" = wesanderson::wes_palettes$Rushmore1[5],
    "Male" = wesanderson::wes_palettes$Rushmore1[1]
  )

ethnicity_color <-
  c(
    "Caucasian" = wesanderson::wes_palettes$Darjeeling2[1],
    "Asian" = wesanderson::wes_palettes$Darjeeling2[2],
    "Hispanics" = wesanderson::wes_palettes$Darjeeling2[3],
    "Black" = wesanderson::wes_palettes$Darjeeling2[4]
  )

iris_color <-
  c(
    "IS" = wesanderson::wes_palettes$Zissou1[5],
    "IR" = wesanderson::wes_palettes$Zissou1[1],
    "Unknown" = "grey"
  )


RColorBrewer::brewer.pal(n = 12, name = "Set3")
omics_color <-
  c(
    transcriptome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    proteomics = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12],
    metabolomics = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
    cytokine = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
    clinical_test = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
    lipidomics = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6],
    gut_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
    skin_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8],
    oral_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[10],
    nasal_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[11]
  )

sample_class_color <-
  c(
    blood = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
    gut_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
    skin_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8],
    oral_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[10],
    nasal_microbiome = RColorBrewer::brewer.pal(n = 12, name = "Set3")[11]
  )

optimize_hc_cluster_number <-
  function(hc,
           data,
           k_list = c(2:20)) {
    purrr::map(k_list, function(k) {
      hc_cluster <- cutree(hc, k = k)
      center <-
        unique(hc_cluster) %>%
        purrr::map(function(cluster) {
          data[which(hc_cluster == cluster),] %>%
            colMeans()
        }) %>%
        dplyr::bind_rows()
      sum(cor(t(center))[upper.tri(cor(t(center)))] > 0.8)
    }) %>%
      unlist()
  }



do_se_swan <-
  function(object,
           qt = "age",
           window_center,
           buckets_size) {
    expression_data <-
      massdataset::extract_expression_data(object)
    
    sample_info <-
      massdataset::extract_sample_info(object)
    
    variable_info <-
      massdataset::extract_variable_info(object)
    
    qt <-
      sample_info %>%
      dplyr::pull(qt)
    
    
    if (missing(window_center)) {
      window_center <-
        quantile(qt, probs = seq(.1, .9, .1))
    }
    
    if (missing(buckets_size)) {
      buckets_size <-
        (max(range(qt)) - min(range(qt))) / 2
    }
    
    p_value <-
      window_center %>%
      purrr::map(function(i) {
        cat(i, " ")
        left_range <-
          c(i - buckets_size / 2, i)
        right_range <-
          c(i, i + buckets_size / 2)
        
        left_idx <-
          which(qt >= left_range[1] &
                  qt < left_range[2])
        
        control_sample_id <-
          colnames(object)[left_idx]
        
        right_idx <-
          which(qt >= right_range[1] &
                  qt < right_range[2])
        
        case_sample_id <-
          colnames(object)[right_idx]
        
        p_value <-
          seq_len(nrow(expression_data)) %>%
          purrr::map(function(i) {
            wilcox.test(as.numeric(expression_data[i, control_sample_id]),
                        as.numeric(expression_data[i, case_sample_id]))$p.value
          }) %>%
          unlist()
        
        p_value <-
          data.frame(variable_id = variable_info$variable_id,
                     p_value) %>%
          dplyr::mutate(p_value_adjust = p.adjust(p_value, "BH"),
                        center = i)
        p_value
        
      }) %>%
      dplyr::bind_rows()
    
    p_value
    
  }



match_sample_id <-
  function(sample_info1,
           sample_info2,
           day_cutoff = 3) {
    sample_info1 <-
      sample_info1 %>%
      dplyr::select(sample_id, subject_id, collection_date)
    
    sample_info2 <-
      sample_info2 %>%
      dplyr::select(sample_id, subject_id, collection_date)
    
    sample_info1 %>%
      t() %>%
      as.data.frame() %>%
      purrr::map(function(x) {
        temp =
          sample_info2 %>%
          dplyr::filter(subject_id == x[2]) %>%
          dplyr::mutate(diff = abs(collection_date - as.Date(x[3]))) %>%
          dplyr::filter(diff == min(diff)) %>%
          dplyr::filter(diff <= day_cutoff) %>%
          head(1)
        if (nrow(temp) == 0) {
          return(NA)
        }
        temp$sample_id
      }) %>%
      unlist() %>%
      unname()
  }

get_go_result_sim <-
  function(result,
           sim.cutoff = 0.7) {
    if (nrow(result) == 0) {
      return(data.frame(
        name1 = character(),
        name2 = character(),
        sim = numeric()
      ))
    }
    
    bp_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "BP"],
                                        ont = "BP",
                                        measure = "Wang") %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)
    
    name <- apply(bp_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })
    
    bp_sim_matrix <-
      bp_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)
    
    mf_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "MF"],
                                        ont = "MF",
                                        measure = "Wang") %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)
    
    name <- apply(mf_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })
    
    mf_sim_matrix <-
      mf_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)
    
    cc_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "CC"],
                                        ont = "CC",
                                        measure = "Wang") %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > 0.7)
    
    name <- apply(cc_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })
    
    cc_sim_matrix <-
      cc_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)
    
    go_sim_matrix <-
      rbind(bp_sim_matrix, mf_sim_matrix, cc_sim_matrix) %>%
      as.data.frame()
    go_sim_matrix
  }


show_matrix_cluster =
  function(result,
           ont = NULL,
           measure = "Wang",
           remove_words,
           margin = 15,
           width = 14,
           height = 8,
           path = "GO_result",
           top = 20) {
    if (all(result$cluster == "Other")) {
      return(NULL)
    }
    
    if (!is.null(ont)) {
      mat =
        simplifyEnrichment::GO_similarity(
          go_id = result %>%
            dplyr::filter(cluster != 'Other') %>%
            dplyr::filter(ONTOLOGY == ont) %>%
            dplyr::pull(node),
          measure = measure,
          ont = ont
        )
      
      if (is.null(nrow(mat))) {
        return(NULL)
      }
    } else{
      mat =
        simplifyEnrichment::term_similarity_from_Reactome(
          term_id = result %>%
            dplyr::filter(cluster != 'Other') %>%
            dplyr::pull(node),
          method = "jaccard"
        )
    }
    
    if (!is.null(ont)) {
      cc_order <-
        result %>%
        dplyr::arrange(cluster) %>%
        dplyr::filter(cluster != 'Other') %>%
        dplyr::filter(ONTOLOGY == ont)
    } else{
      cc_order <-
        result %>%
        dplyr::arrange(cluster) %>%
        dplyr::filter(cluster != 'Other')
    }
    
    cc_order =
      cc_order[stringr::str_order(cc_order$cluster, numeric = TRUE),]
    
    ####only remain the top 20 cluster2 with smallest P values
    temp_p =
      cc_order %>%
      plyr::dlply(.variables = .(Direction)) %>%
      purrr::map(
        .f = function(x) {
          x %>%
            plyr::dlply(.variables = .(cluster)) %>%
            purrr::map(function(x) {
              min(x$p.adjust)
            }) %>%
            unlist() %>%
            sort() %>%
            head(top)
        }
      )
    
    temp_p =
      temp_p[c("UP", 'DOWN')] %>%
      unlist()
    
    names(temp_p) =
      names(temp_p) %>%
      stringr::str_replace_all("UP\\.", "") %>%
      stringr::str_replace_all("DOWN\\.", "")
    
    
    cc_order =
      cc_order %>%
      dplyr::filter(cluster %in% names(temp_p))
    
    mat = mat[cc_order$node, cc_order$node]
    
    keywords <-
      cc_order %>%
      plyr::dlply(.variables = .(cluster)) %>%
      purrr::map(
        .f = function(x) {
          x <-
            x %>%
            dplyr::select(Description, p.adjust) %>%
            dplyr::mutate(p.adjust = -log(as.numeric(p.adjust), 10)) %>%
            plyr::dlply(.variables = .(Description)) %>%
            purrr::map(
              .f = function(y) {
                data.frame(
                  word = stringr::str_split(y$Description, " ")[[1]],
                  p.adjust = y$p.adjust
                )
              }
            ) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            dplyr::arrange(word)
          
          rownames(x) <- NULL
          
          x <-
            x %>%
            dplyr::group_by(word) %>%
            dplyr::summarise(freq = sum(as.numeric(p.adjust))) %>%
            dplyr::ungroup() %>%
            dplyr::filter(!word %in% remove_words)
          
        }
      )
    
    keywords =
      keywords[unique(cc_order$cluster)]
    
    
    name <-
      lapply(keywords, function(x) {
        !is.null(x)
      }) %>%
      unlist() %>% which() %>%
      names()
    
    keywords =
      keywords[names(keywords) %in% name]
    
    library(ComplexHeatmap)
    library(circlize)
    col_fun = colorRamp2(c(0, 1), c("white", "red"))
    
    cluster_col <-
      unique(result$cluster)[unique(result$cluster) != 'Other']
    
    cluster_col <-
      colorRampPalette(colors = RColorBrewer::brewer.pal(n = 10, name = "Spectral")[c(7, 2, 6, 10, 9, 1, 5, 8, 4, 3)])(n = length(cluster_col))
    
    names(cluster_col) <-
      unique(result$cluster)[unique(result$cluster) != 'Other']
    
    direction_col = cluster_col[unique(cc_order$cluster)]
    direction_col[names(direction_col) %in% cc_order$cluster[cc_order$Direction == "UP"]] =
      ggsci::pal_aaas()(n = 10)[2]
    
    direction_col[names(direction_col) %in% cc_order$cluster[cc_order$Direction == "DOWN"]] =
      ggsci::pal_aaas()(n = 10)[1]
    
    ha1 =
      rowAnnotation(
        foo = anno_block(
          gp = gpar(fill = direction_col),
          labels = NULL,
          labels_gp = gpar(col = "white", fontsize = 10),
          width = unit(5, "mm")
        ),
        foo2 = anno_block(
          gp = gpar(fill = unique(cluster_col[cc_order$cluster])),
          labels = unique(stringr::str_replace(cc_order$cluster, "Cluster ", "")),
          labels_gp = gpar(col = "white", fontsize = 10),
          width = unit(5, "mm")
        )
      )
    
    ha2 =
      HeatmapAnnotation(foo = anno_block(
        gp = gpar(fill = unique(cluster_col[cc_order$cluster])),
        labels = unique(stringr::str_replace(cc_order$cluster, "Cluster ", "")),
        labels_gp = gpar(col = "white", fontsize = 10),
        height = unit(5, "mm")
      ))
    
    ht = Heatmap(
      mat,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      name = "GO similarity",
      col = col_fun,
      right_annotation = ha1,
      top_annotation = ha2,
      border = FALSE,
      row_split = as.numeric(factor(
        cc_order$cluster, levels = unique(cc_order$cluster)
      )),
      column_split = as.numeric(factor(
        cc_order$cluster, levels = unique(cc_order$cluster)
      )),
      column_title = NULL,
      row_title = NULL,
      row_gap = unit(0.5, "mm"),
      column_gap = unit(0.5, "mm")
    )
    
    cl =
      factor(cc_order$cluster, levels = unique(cc_order$cluster))
    
    align_to = split(seq_len(nrow(mat)), cl)
    align_to = align_to[names(align_to) %in% names(keywords)]
    align_to
    
    fontsize_range = c(8, 16)
    
    gbl = lapply(names(align_to), function(nm) {
      kw = keywords[[nm]]$word
      freq = keywords[[nm]]$freq
      fontsize = scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)
      word_cloud_grob(
        text = kw,
        fontsize = fontsize,
        col = colorRampPalette(colors = ggsci::pal_aaas()(n =
                                                            10))(length(kw))
      )
    })
    
    names(gbl) = names(align_to)
    
    # margin = unit(8, "pt")
    margin = unit(margin, "pt")
    gbl_h = lapply(gbl, function(x)
      convertHeight(grobHeight(x), "cm") + margin)
    
    gbl_h = do.call(unit.c, gbl_h)
    
    gbl_w = lapply(gbl, function(x)
      convertWidth(grobWidth(x), "cm"))
    
    gbl_w = do.call(unit.c, gbl_w)
    
    gbl_w = max(gbl_w) + margin
    
    panel_fun = function(index, nm) {
      # background
      grid.rect(gp = gpar(fill = "#DDDDDD", col = NA))
      # border
      grid.lines(
        c(0, 1, 1, 0),
        c(0, 0, 1, 1),
        gp = gpar(col = "#AAAAAA"),
        default.units = "npc"
      )
      gb = gbl[[nm]]
      # a viewport within the margins
      pushViewport(
        viewport(
          x = margin / 2,
          y = margin / 2,
          width = grobWidth(gb),
          height = grobHeight(gb),
          just = c("left", "bottom")
        )
      )
      grid.draw(gb)
      popViewport()
    }
    
    ht = ht + rowAnnotation(
      keywords = anno_link(
        align_to = align_to,
        which = "row",
        panel_fun = panel_fun,
        size = gbl_h,
        gap = unit(2, "mm"),
        width = gbl_w + unit(5, "mm"),
        # 5mm for the link
        link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"),
        internal_line = FALSE
      )
    ) # you can set it to TRUE to see what happens
    
    ht
    
    # if (!is.null(ont)) {
    #   height = case_when(ont == "BP" ~ 12,
    #                      ont == "MF" ~ 5,
    #                      ont == "CC" ~ 6)
    # } else{
    #   height = height
    # }
    
    # ht
    
    pdf(file = file.path(path, paste(ont, "_sim_matrix.pdf", sep = "")),
        width = width,
        height = height)
    draw(ht, ht_gap = unit(2, "pt"))
    dev.off()
  }


get_jaccard_index_for_three_databases <-
  function(result_go_cluster,
           result_kegg_cluster,
           result_reactome_cluster,
           variable_info) {
    met_data =
      rbind(result_go_cluster,
            result_kegg_cluster,
            result_reactome_cluster)
    
    temp_data <-
      met_data$geneID %>%
      stringr::str_split("/") %>%
      purrr::map(function(x) {
        if (stringr::str_detect(x[1], "ENSG")) {
          return(x)
        }
        
        if (stringr::str_detect(x[1], "[A-Za-z]")) {
          return(variable_info$ENSEMBL[match(x, variable_info$UNIPROT)])
        }
        
        return(variable_info$ENSEMBL[match(x, variable_info$ENTREZID)])
        
      })
    
    names(temp_data) =
      met_data$module_annotation
    
    ##calculate jaccard index
    jaccard_index <-
      purrr::map(
        1:(length(temp_data) - 1),
        .f = function(idx) {
          temp <-
            purrr::map(
              temp_data[(idx + 1):length(temp_data)],
              .f = function(y) {
                length(intersect(temp_data[[idx]], y)) / length(union(temp_data[[idx]], y))
              }
            ) %>%
            unlist() %>%
            data.frame(value = .) %>%
            tibble::rownames_to_column(var = "name2") %>%
            data.frame(name1 = names(temp_data)[idx], .)
          temp$name2 <-
            names(temp_data[(idx + 1):length(temp_data)])
          temp
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    return(jaccard_index)
    
  }

database_color =
  c(
    GO = ggsci::pal_d3()(n = 10)[1],
    KEGG = ggsci::pal_d3()(n = 10)[2],
    Reactome = ggsci::pal_d3()(n = 10)[3],
    HMDB = ggsci::pal_d3()(n = 10)[5]
  )

calculate_jaccard_index <-
  function(x, y) {
    x <- unique(x[!is.na(x)])
    y <- unique(y[!is.na(y)])
    
    length(intersect(x, y)) / length(unique(c(x, y)))
    
  }



cytokine_class_color <-
  c(
    "Proinflammatory" = ggsci::pal_nejm()(n = 9)[1],
    "Anti-inflammatory" = ggsci::pal_nejm()(n = 9)[2],
    "Proinflammatory/Anti-inflammatory" = ggsci::pal_nejm()(n = 9)[3]
  )


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

genus_color <-
  c(
    "Acidobacteria",
    "Actinobacteria",
    "Bacteroidetes",
    "Campilobacterota",
    "Cyanobacteria/Chloroplast",
    "Deinococcus-Thermus",
    "Firmicutes",
    "Fusobacteria",
    "Gemmatimonadetes",
    "Lentisphaerae",
    "Candidatus_Saccharibacteria",
    "Proteobacteria",
    "Plantae",
    "Spirochaetes",
    "SR1",
    "Synergistetes",
    "Unclassified_Bacteria",
    "Verrucomicrobia"
  )

names(genus_color) <-
  genus_color

genus_color <- gg_color_hue(length(genus_color))

combine_workds <- function(y) {
  y[y == "acids"] <- "acid"
  y[y == "Intrinsic"] <- "intrinsic"
  y[y == "Fatty"] <- "fatty"
  y[y == "Double"] <- "double"
  y[y == "Diseases"] <- "disease"
  y[y == "Activation"] <- "activation"
  y[y == "Antigen"] <- "antigen"
  y[y == "Assembly"] <- "assembly"
  y[y == "Autophagy"] <- "autophagy"
  y[y == "Cascades"] <- "cascade"
  y[y == "Cellular"] <- "cellular"
  y[y == "compounds"] <- "compound"
  y[y == "Damage"] <- "damaged"
  y[y == "genes"] <- "gene"
  y[y == "Genes"] <- "gene"
  y[y == "heme"] <- "Heme"
  y[y == "Immune"] <- "immune"
  y[y == "Infections"] <- "Infection"
  y[y == "Initiation"] <- "initiation"
  y[y == "innate"] <- "initiation"
  y[y == "iron"] <- "Iron"
  y[y == "mucosal"] <- "mucosa"
  y[y == "process"] <- "Processing"
  y[y == "regulate"] <- "regulation"
  y[y == "regulated"] <- "regulation"
  y[y == "regulates"] <- "regulation"
  y[y == "Regulates"] <- "regulation"
  y[y == "Regulation"] <- "regulation"
  y[y == "regulators"] <- "regulation"
  y[y == "Replication"] <- "replication"
  y[y == "signal"] <- "signaling"
  y[y == "Signaling"] <- "signaling"
  y[y == "Stress"] <- "stress"
  y[y == "Translation"] <- "translation"
  y[y == "Viral"] <- "virus"
  y[y == "viral"] <- "virus"
  y[stringr::str_detect(y, "ubiquinol|ubiquitin|Ubiquitin")] <-
    "ubiquitination"
  y[y == "vesicle"] <- "Vesicle"
  y[y == "transcription"] <- "Transcription"
  y[y == "Transcriptional"] <- "Transcription"
  y[y == "subunits"] <- "subunit"
  y[y == "Strand"] <- "strand"
  y[y == "Repair"] <- "repair"
  y[y == "Receptor"] <- "receptor"
  y[y == "receptor-beta"] <- "receptor"
  y[y == "receptors"] <- "receptor"
  y[y == "respiratory"] <- "respiration"
  y[y == "Respiratory"] <- "Respiratory"
  y[y == "Response"] <- "response"
  y[y == "oxidant"] <- "oxidation"
  y[y == "oxidative"] <- "oxidation"
  y[y == "Oxidative"] <- "oxidation"
  y[y == "oxidative"] <- "oxidation"
  y[y == "oxide"] <- "oxidation"
  y[y == "oxygen"] <- "oxidation"
  y[y == "translational"] <- "translation"
  y[y == "toll-like"] <- "Toll-like"
  y[y == "stimulates"] <- "stimulus"
  y[y == "modified"] <- "modification"
  y[y == "modulates"] <- "modification"
  y[y == "Mitotic"] <- "mitotic"
  y[y == "Metabolic"] <- "metabolism"
  y[y == "Metabolism"] <- "metabolism"
  y[y == "metabolic"] <- "metabolism"
  y[y == "metabolic"] <- "derivation"
  y[y == "derivative"] <- "derivation"
  y[y == "derivatives"] <- "derivation"
  y[y == "coupled"] <- "coupling"
  y[y == "Complex"] <- "complex"
  y[y == "activate"] <- "activation"
  y[y == "Activated"] <- "activation"
  y[y == "activity"] <- "activation"
  y <- y[!stringr::str_detect(y, "SARS")]
  y
}

######this is used for word cloud
remove_words <-
  c(
    "to",
    "of",
    "in",
    "type",
    "pathway",
    "IX",
    "part",
    "life",
    "control",
    "quality",
    "body",
    "late",
    "cell",
    "species",
    "cells",
    "or",
    "levels",
    "as",
    "on",
    "by",
    "small",
    "other",
    "involved",
    "alpha",
    "specific",
    "number",
    "through",
    "outer",
    "large",
    "rough",
    "early",
    "via",
    "smooth",
    "system",
    "into",
    "entry",
    "and",
    "T",
    "based",
    "within",
    "from",
    "built",
    "mediated",
    "-",
    "_",
    "animal",
    "the",
    "free",
    "a",
    "pool",
    "60S",
    "40S",
    "and",
    "chain",
    "Decay",
    "enhanced",
    "independent",
    "joining",
    "4",
    "2",
    "up",
    "take",
    "release",
    'Like',
    "presentation",
    "Class",
    "I",
    "mediated",
    "exchange",
    "&",
    "events",
    "B",
    "an",
    "",
    "at",
    "B",
    "Base",
    "c",
    "E",
    "during",
    "for",
    "Major",
    "NOTCH",
    "Of",
    "Opening",
    "Pathway",
    "processing",
    "free",
    letters,
    LETTERS,
    "family",
    "them",
    "ii",
    "class",
    1:7,
    "group",
    "phase",
    "ar",
    "orc",
    "new",
    "ap",
    "ends",
    "sars-cov-2",
    "upon",
    "ix",
    "major",
    "System",
    "with",
    "affected",
    "along",
    "AP",
    "AR",
    "associated",
    "Associated",
    "association",
    "containing",
    "down",
    "Ends",
    "II",
    "SARS-CoV",
    "ARS-CoV-2-host",
    "Network",
    "virus",
    "regulation",
    "Processing",
    "protein",
    "Protein",
    "Nonsense",
    "Nonsense-Mediated",
    "production",
    "pathways",
    "multiple",
    "scanning",
    "site",
    "The",
    "start",
    "pattern",
    "Processing",
    "Phase",
    "Packaging"
  )


pahtway_heatmap_all_cluster <-
  function(result_all_list,
           word_cloud_list,
           gene_marker,
           expression_data,
           top_pathway = 10,
           font_size = c(5, 15),
           word_count_breaks = c(1, 30),
           show_column_names = FALSE,
           add_text = FALSE) {
    library(plyr)
    library(tidygraph)
    library(ggraph)
    library(igraph)
    library(ComplexHeatmap)
    library(patchwork)
    library(sjPlot)
    
    temp_info <-
      seq_along(result_all_list) %>%
      purrr::map(function(i) {
        result_all <-
          result_all_list[[i]] %>%
          dplyr::mutate(cluster = names(result_all_list)[i]) %>%
          dplyr::arrange(p.adjust) %>%
          head(top_pathway)
        result_all <-
          result_all %>%
          dplyr::arrange(p.adjust) %>%
          plyr::dlply(.variables = .(module_annotation)) %>%
          purrr::map(
            .f = function(x) {
              gene_id <-
                stringr::str_split(x$geneID, pattern = "/")[[1]]
              
              x <-
                data.frame(x %>% dplyr::select(-geneID),
                           gene_id, stringsAsFactors = FALSE)
              
              gene_id1 =
                gene_marker$variable_id[match(x$gene_id, gene_marker$ENSEMBL)]
              gene_id2 =
                gene_marker$variable_id[match(x$gene_id, gene_marker$UNIPROT)]
              gene_id3 =
                gene_marker$variable_id[match(x$gene_id, gene_marker$ENTREZID)]
              
              gene_id =
                data.frame(gene_id1, gene_id2, gene_id3) %>%
                apply(1, function(x) {
                  x[!is.na(x)][1]
                })
              
              x$gene_id <-
                gene_id
              
              x$correlation <-
                gene_marker$correlation[match(x$gene_id, gene_marker$variable_id)]
              x$p_value <-
                gene_marker$p_value[match(x$gene_id, gene_marker$variable_id)]
              x$Count <- nrow(x)
              x
            }
          ) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::select(cluster,
                        module_annotation,
                        p.adjust,
                        correlation,
                        p_value,
                        Count,
                        gene_id)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    rownames(temp_info) <- NULL
    
    temp_info$id <-
      paste(temp_info$cluster,
            temp_info$module_annotation, sep = "_")
    
    ###quantitative pathways
    ##heatmap
    temp_data <-
      seq_along(result_all_list) %>%
      purrr::map(function(i) {
        result_all <-
          result_all_list[[i]] %>%
          dplyr::mutate(cluster = as.character(i)) %>%
          dplyr::arrange(p.adjust) %>%
          head(top_pathway)
        
        pathway <-
          result_all %>%
          plyr::dlply(.variables = .(module_annotation)) %>%
          purrr::map(
            .f = function(x) {
              name <- x$module_annotation[1]
              cluster <- names(result_all_list)[i]
              gene <- x$geneID %>%
                stringr::str_split("/") %>%
                `[[`(1)
              
              gene1 =
                gene_marker$variable_id[match(gene, gene_marker$ENSEMBL)]
              gene2 =
                gene_marker$variable_id[match(gene, gene_marker$ENTREZID)]
              gene3 =
                gene_marker$variable_id[match(gene, gene_marker$UNIPROT)]
              
              gene =
                data.frame(gene1, gene2, gene3) %>%
                apply(1, function(x) {
                  x[!is.na(x)][1]
                })
              
              data <- expression_data[gene,] %>%
                as.data.frame() %>%
                apply(2, mean) %>%
                as.data.frame() %>%
                t() %>%
                as.data.frame()
              
              data <-
                data.frame(
                  name,
                  cluster = cluster,
                  data,
                  stringsAsFactors = FALSE,
                  check.names = FALSE
                )
              
              data
            }
          ) %>%
          do.call(rbind, .) %>%
          as.data.frame()
        
        pathway_data <- pathway %>%
          dplyr::select(-c(name, cluster))
        
        range(pathway_data)
        
        if (abs(range(pathway_data)[1]) > range(pathway_data)[2]) {
          pathway_data[pathway_data < -range(pathway_data)[2]] <-
            -range(pathway_data)[2]
        } else{
          pathway_data[pathway_data > -range(pathway_data)[1]] <-
            -range(pathway_data)[1]
        }
        
        rownames(pathway) <-
          rownames(pathway_data) <-
          paste(pathway$cluster, pathway$name, sep = "_")
        
        list(pathway = pathway, pathway_data = pathway_data)
      })
    
    pathway <-
      temp_data %>%
      purrr::map(function(x) {
        x$pathway
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    pathway_data <-
      temp_data %>%
      purrr::map(function(x) {
        x$pathway_data
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    rownames(pathway) <-
      rownames(pathway_data) <-
      paste(pathway$cluster, pathway$name, sep = "_")
    
    lm_reg <-
      lm(formula = y ~ x,
         data = data.frame(y = font_size, x = word_count_breaks))
    
    text <-
      word_cloud_list %>%
      purrr::map(function(x) {
        z <-
          data.frame(text = x) %>%
          dplyr::group_by(text) %>%
          dplyr::summarise(fontsize = dplyr::n()) %>%
          as.data.frame()
        if (nrow(z) > 50) {
          z <- z %>%
            dplyr::filter(fontsize > 1)
        }
        z$fontsize <-
          predict(lm_reg, newdata = data.frame(x = z$fontsize))
        z %>%
          dplyr::arrange(dplyr::desc(fontsize)) %>%
          head(20)
      })
    
    names(text) <- unique(pathway$cluster)
    
    align_to = unique(pathway$cluster) %>%
      lapply(function(z) {
        which(pathway$cluster == z)
      })
    
    names(align_to) <-
      unique(pathway$cluster)
    
    col <-
      RColorBrewer::brewer.pal(n = 9, name = "BrBG")[1:length(unique(pathway$cluster))]
    
    names(col) <-
      unique(pathway$cluster)
    
    ###important genes heatmap
    library(circlize)
    col_fun <-
      colorRamp2(seq(min(as.matrix(pathway_data)),
                     max(as.matrix(pathway_data)), length.out = 9),
                 rev(RColorBrewer::brewer.pal(n = 9, name = "BrBG")))
    
    if (add_text) {
      ha_left <-
        rowAnnotation(cluster = anno_block(
          align_to = align_to,
          panel_fun = function(index, nm) {
            grid.rect(gp = gpar(fill = col[nm]))
            grid.text(stringr::str_replace(nm, "Cluster", ""),
                      0.5, 0.5)
          },
          width = unit(0.3, "cm")
        ))
      
      cor <-
        temp_info %>%
        dplyr::select(id, correlation, gene_id) %>%
        tidyr::pivot_wider(names_from = "gene_id",
                           values_from = "correlation") %>%
        as.data.frame()
      
      rownames(cor) <- cor$id
      cor <-
        cor %>%
        dplyr::select(-id) %>%
        as.matrix()
      
      mean_cor <-
        apply(as.data.frame(cor), 1, function(x) {
          mean(x, na.rm = TRUE)
        })
      
      cor_col <-
        circlize::colorRamp2(breaks = c(-1, 0, 1),
                             colors = c("darkblue", "white", "red"))
      
      ha_right <-
        rowAnnotation(
          Correlation = anno_boxplot(
            x = cor,
            width = unit(1, "cm"),
            outline = FALSE,
            gp = gpar(fill = cor_col(mean_cor))
          ),
          Count = anno_points(
            x = temp_info %>%
              dplyr::distinct(id, .keep_all = TRUE)
            %>% pull(Count),
            width = unit(1, "cm")
          ),
          textbox = anno_textbox(
            align_to = pathway$cluster,
            text = text,
            max_width = unit(40, "mm"),
            word_wrap = FALSE
          )
        )
      
      rownames(pathway_data) <-
        rownames(pathway_data) %>%
        stringr::str_replace("cluster", "")
      
      plot2 <-
        Heatmap(
          as.matrix(pathway_data),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(cex = 0.7),
          show_column_names = show_column_names,
          border = TRUE,
          col = col_fun,
          name = "z-score",
          clustering_method_rows = "ward.D",
          column_names_rot = 45,
          column_title = NULL,
          column_names_gp = gpar(cex = 0.5),
          row_split = pathway$cluster,
          row_title = NULL,
          left_annotation = ha_left,
          right_annotation = ha_right
        )
      
    } else{
      ha_right <-
        rowAnnotation(score = anno_points(x = score),
                      recover_score = anno_points(x = recover_score))
      
      ha_left <-
        rowAnnotation(cluster = anno_block(
          align_to = align_to,
          panel_fun = function(index, nm) {
            grid.rect(gp = gpar(fill = col[nm]))
            grid.text(nm, 0.5, 0.5)
          },
          width = unit(0.5, "cm")
        ))
      
      plot2 <-
        Heatmap(
          as.matrix(pathway_data),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_names = FALSE,
          show_column_names = show_column_names,
          border = TRUE,
          col = col_fun,
          name = "z-score",
          clustering_method_rows = "ward.D",
          column_names_rot = 45,
          column_split = c(rep(1, sum(
            colnames(pathway_data) != "PP"
          )),
          rep(2, 1)),
          column_title = NULL,
          row_split = pathway$cluster,
          row_title = NULL,
          right_annotation = ha_right,
          left_annotation = ha_left
        )
    }
    
    library(ggplotify)
    
    plot2 <- as.ggplot(plot2)
    plot2
  }


mci_test <- function(x, y, times = 1000) {
  mic <- minerva::mine(x, y)$MIC
  ####permutation to get the p values
  null_distributation <-
    purrr::map(1:times, function(i) {
      minerva::mine(x, sample(y, length(y), replace = FALSE))$MIC
    }) %>%
    unlist()
  ####get p value
  mean_value <-
    mean(null_distributation)
  sd_value <-
    sd(null_distributation)
  normal_distributation <-
    rnorm(n = 100000, mean = mean_value, sd = sd_value)
  p_value <-
    1 - sum(mic > normal_distributation) / 100000
  rm("normal_distributation")
  gc()
  data.frame(mic = mic,
             p_value = p_value)
}




network_class_color <-
  c(
    "function" = ggsci::pal_futurama()(n = 4)[1],
    "module" = ggsci::pal_futurama()(n = 4)[2],
    "pathway" = ggsci::pal_futurama()(n = 4)[3],
    "gene" = "#8DD3C7",
    "protein" = "#FFED6F",
    "metabolite" = "#BEBADA"
  )



new_coords <- function(range_left = 0,
                       range_right = 100,
                       old_x = c(5, 60, 10),
                       scale = FALSE) {
  if (scale) {
    temp = data.frame(
      x = c(min(old_x), max(old_x)),
      y = c(range_left, range_right)
    )
    lm_reg <- lm(y ~ x, data = temp)
    temp <-
    coefficients(lm_reg)[2] * old_x + coefficients(lm_reg)[1]
    return(temp)
  }
  if (length(old_x) == 1) {
    return(mean(c(range_left, range_right)))
  }
  
  if (length(old_x) == 2) {
    new_x <- data.frame(old_x,
                        rank = rank(old_x),
                        original_index = seq_along(old_x)) %>%
      dplyr::arrange(rank) %>%
      dplyr::mutate(new_x = c(range_left, range_right)) %>%
      dplyr::arrange(original_index) %>%
      pull(new_x)
    return(new_x)
  }
  
  if (length(old_x) > 2) {
    new_x <- data.frame(old_x,
                        rank = rank(old_x),
                        original_index = seq_along(old_x)) %>%
      dplyr::arrange(rank) %>%
      dplyr::mutate(new_x = seq(
        from = range_left,
        to = range_right,
        length.out = length(old_x)
      )) %>%
      dplyr::arrange(original_index) %>%
      pull(new_x)
    return(new_x)
  }
  
}
