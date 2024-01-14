no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)
library(microbiomedataset)

###load cluster in 2,4,6
cluster_info2 <-
  readxl::read_xlsx(
    "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_2/cluster2.xlsx"
  )

cluster_info4 <-
  readxl::read_xlsx(
    "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_4/cluster4.xlsx"
  )

cluster_info5 <-
  readxl::read_xlsx(
    "3-data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cluster_5/cluster5.xlsx"
  )

load("3-data_analysis/combined_omics/DE_SWAN/temp_data")

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
  dplyr::filter(p_value_adjust < 0.05)

variable_id1 <-
  cluster_info2$variable_id

variable_id_cluster2 <-
  cluster_info2$variable_id

variable_id_cluster4 <-
  cluster_info4$variable_id

variable_id_cluster5 <-
  cluster_info5$variable_id

variable_id_crest1 <-
  temp_data %>%
  dplyr::filter(center == 44) %>%
  pull(variable_id)

variable_id_crest2 <-
  temp_data %>%
  dplyr::filter(center == 60) %>%
  pull(variable_id)

library(ComplexHeatmap)

lt <-
  list(
    crest1 = variable_id_crest1,
    crest2 = variable_id_crest2,
    cluster2 = variable_id_cluster2,
    cluster4 = variable_id_cluster4,
    cluster5 = variable_id_cluster5
  )
m = make_comb_mat(lt)

plot <-
  UpSet(
    m,
    set_order = c("crest1", "crest2", "cluster2", "cluster4", "cluster5"),
    comb_col = c(ggsci::pal_jama()(n = 3))[comb_degree(m)],
    bg_col = "#F0F0FF",
    bg_pt_col = "grey"
  )

plot <- ggplotify::as.ggplot(plot)

plot
dir.create("3-data_analysis/combined_omics/crest_cluster_overlap")
setwd("3-data_analysis/combined_omics/crest_cluster_overlap")

ggsave(plot,
       filename = "all_overlap.pdf",
       width = 7,
       height = 4)


library(ggVennDiagram)
library(scales)
library(RColorBrewer)
plot <-
ggVennDiagram::ggVennDiagram(x = list(crest1 = variable_id_crest1,
                                      crest2 = variable_id_crest2,
                                      cluster2 = variable_id_cluster2)) +
  scale_fill_continuous(type = "gradient",
                        low=colorRampPalette(brewer.pal(6, "Reds"))(10)[1],
                        high = colorRampPalette(brewer.pal(6, "Reds"))(10)[10],
                        limits = c(0,3400))
plot
ggsave(plot,
       filename = "crest1_crest2_cluster2_overlap.pdf",
       width = 7,
       height = 7)


plot <-
  ggVennDiagram::ggVennDiagram(x = list(crest1 = variable_id_crest1,
                                        crest2 = variable_id_crest2,
                                        cluster4 = variable_id_cluster4)) +
  scale_fill_continuous(type = "gradient",
                        low=colorRampPalette(brewer.pal(6, "Reds"))(10)[1],
                        high = colorRampPalette(brewer.pal(6, "Reds"))(10)[10],
                        limits = c(0,3400))
plot
ggsave(plot,
       filename = "crest1_crest2_cluster4_overlap.pdf",
       width = 7,
       height = 7)

plot <-
  ggVennDiagram::ggVennDiagram(x = list(crest1 = variable_id_crest1,
                                        crest2 = variable_id_crest2,
                                        cluster5 = variable_id_cluster5)) +
  scale_fill_continuous(type = "gradient",
                        low=colorRampPalette(brewer.pal(6, "Reds"))(10)[1],
                        high = colorRampPalette(brewer.pal(6, "Reds"))(10)[10],
                        limits = c(0,3400))
plot
ggsave(plot,
       filename = "crest1_crest2_cluster5_overlap.pdf",
       width = 7,
       height = 7)

