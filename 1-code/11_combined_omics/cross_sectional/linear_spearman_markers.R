no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)
library(microbiomedataset)

###load("data)
load(
  "3-data_analysis/plasma_transcriptome/linear_spearman/cross_section/variable_info"
)

transcriptome_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = SYMBOL)

load("3-data_analysis/plasma_proteomics/linear_spearman/cross_section/variable_info")
proteomics_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = variable_id)

load(
  "3-data_analysis/plasma_metabolomics/linear_spearman/cross_section/variable_info"
)

metabolomics_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = Compound.name)

load("3-data_analysis/plasma_cytokine/linear_spearman/cross_section/variable_info")
cytokine_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = SYMBOL)

load("3-data_analysis/clinical_test/linear_spearman/cross_section/variable_info")
clinical_test_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = test_name)

load("3-data_analysis/plasma_lipidomics/linear_spearman/cross_section/variable_info")
lipidomics_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = variable_id)

load("3-data_analysis/gut_microbiome/linear_spearman/cross_section/variable_info")
gut_microbiome_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = Genus)

load("3-data_analysis/skin_microbiome/linear_spearman/cross_section/variable_info")
skin_microbiome_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = Genus)

load("3-data_analysis/oral_microbiome/linear_spearman/cross_section/variable_info")
oral_microbiome_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = Genus)

load("3-data_analysis/nasal_microbiome/linear_spearman/cross_section/variable_info")
nasal_microbiome_variable_info <-
  variable_info %>%
  dplyr::mutate(mol_name = Genus)

load("3-data_analysis/plasma_transcriptome/DEG/cross_section_loess/all_marker_name")
transcriptome_all_marker_name <- all_marker_name

load("3-data_analysis/plasma_proteomics/DEG/cross_section_loess/all_marker_name")
proteomics_all_marker_name <- all_marker_name

load("3-data_analysis/plasma_metabolomics/DEG/cross_section_loess/all_marker_name")
metabolomics_all_marker_name <- all_marker_name

load("3-data_analysis/plasma_metabolomics/DEG/cross_section_loess/all_marker_name")
metabolomics_all_marker_name <- all_marker_name

load("3-data_analysis/plasma_cytokine/DEG/cross_section_loess/all_marker_name")
cytokine_all_marker_name <- all_marker_name

load("3-data_analysis/clinical_test/DEG/cross_section_loess/all_marker_name")
clinical_test_all_marker_name <- all_marker_name

load("3-data_analysis/plasma_lipidomics/DEG/cross_section_loess/all_marker_name")
lipidomics_all_marker_name <- all_marker_name

load("3-data_analysis/gut_microbiome/DEG/cross_section_loess/all_marker_name")
gut_microbiome_all_marker_name <- all_marker_name

load("3-data_analysis/skin_microbiome/DEG/cross_section_loess/all_marker_name")
skin_microbiome_all_marker_name <- all_marker_name

load("3-data_analysis/oral_microbiome/DEG/cross_section_loess/all_marker_name")
oral_microbiome_all_marker_name <- all_marker_name

load("3-data_analysis/nasal_microbiome/DEG/cross_section_loess/all_marker_name")
nasal_microbiome_all_marker_name <- all_marker_name

load("3-data_analysis/plasma_transcriptome/DEG/cross_section_loess/all_marker_name")
transcriptome_all_marker_name_permutation <-
  all_marker_name

load(
  "3-data_analysis/plasma_proteomics/DEG/cross_section_loess/all_marker_name_permutation"
)
proteomics_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/plasma_metabolomics/DEG/cross_section_loess/all_marker_name_permutation"
)
metabolomics_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/plasma_metabolomics/DEG/cross_section_loess/all_marker_name_permutation"
)
metabolomics_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/plasma_cytokine/DEG/cross_section_loess/all_marker_name_permutation"
)
cytokine_all_marker_name_permutation <- all_marker_name_permutation

load(
  "3-data_analysis/clinical_test/DEG/cross_section_loess/all_marker_name_permutation"
)
clinical_test_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/plasma_lipidomics/DEG/cross_section_loess/all_marker_name_permutation"
)
lipidomics_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/gut_microbiome/DEG/cross_section_loess/all_marker_name_permutation"
)
gut_microbiome_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/skin_microbiome/DEG/cross_section_loess/all_marker_name_permutation"
)
skin_microbiome_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/oral_microbiome/DEG/cross_section_loess/all_marker_name_permutation"
)
oral_microbiome_all_marker_name_permutation <-
  all_marker_name_permutation

load(
  "3-data_analysis/nasal_microbiome/DEG/cross_section_loess/all_marker_name_permutation"
)
nasal_microbiome_all_marker_name_permutation <-
  all_marker_name_permutation



dir.create("3-data_analysis/combined_omics/linear_spearman_marker",
           recursive = TRUE)

setwd("3-data_analysis/combined_omics/linear_spearman_marker")

sum(transcriptome_variable_info$cor_p < 0.05)
sum(transcriptome_variable_info$cor_p < 0.05) / nrow(transcriptome_variable_info)
sum(transcriptome_variable_info$lm_p < 0.05)
length(transcriptome_all_marker_name)

sum(proteomics_variable_info$cor_p < 0.05)
sum(proteomics_variable_info$cor_p < 0.05) / nrow(proteomics_variable_info)
sum(proteomics_variable_info$lm_p < 0.05)
length(proteomics_all_marker_name)

sum(metabolomics_variable_info$cor_p < 0.05)
sum(metabolomics_variable_info$cor_p < 0.05) / nrow(metabolomics_variable_info)
sum(metabolomics_variable_info$lm_p < 0.05)
length(metabolomics_all_marker_name)

sum(cytokine_variable_info$cor_p < 0.05)
sum(cytokine_variable_info$cor_p < 0.05) / nrow(cytokine_variable_info)
sum(cytokine_variable_info$lm_p < 0.05)
length(cytokine_all_marker_name)

sum(clinical_test_variable_info$cor_p < 0.05)
sum(clinical_test_variable_info$cor_p < 0.05) / nrow(clinical_test_variable_info)
sum(clinical_test_variable_info$lm_p < 0.05)
length(clinical_test_all_marker_name)

sum(lipidomics_variable_info$cor_p < 0.05)
sum(lipidomics_variable_info$cor_p < 0.05) / nrow(lipidomics_variable_info)
sum(lipidomics_variable_info$cor_p < 0.05)
sum(lipidomics_variable_info$lm_p < 0.05)
length(lipidomics_all_marker_name)

sum(gut_microbiome_variable_info$cor_p < 0.05)
sum(gut_microbiome_variable_info$cor_p < 0.05) / nrow(gut_microbiome_variable_info)
sum(gut_microbiome_variable_info$lm_p < 0.05)
length(gut_microbiome_all_marker_name)

sum(skin_microbiome_variable_info$cor_p < 0.05)
sum(skin_microbiome_variable_info$cor_p < 0.05) / nrow(skin_microbiome_variable_info)
sum(skin_microbiome_variable_info$lm_p < 0.05)
length(skin_microbiome_all_marker_name)

sum(oral_microbiome_variable_info$cor_p < 0.05)
sum(oral_microbiome_variable_info$cor_p < 0.05) / nrow(oral_microbiome_variable_info)
sum(oral_microbiome_variable_info$lm_p < 0.05)
length(oral_microbiome_all_marker_name)

sum(nasal_microbiome_variable_info$cor_p < 0.05)
sum(nasal_microbiome_variable_info$cor_p < 0.05) / nrow(nasal_microbiome_variable_info)
sum(nasal_microbiome_variable_info$lm_p < 0.05)
length(nasal_microbiome_all_marker_name)

length(transcriptome_all_marker_name) +
  length(proteomics_all_marker_name) +
  length(metabolomics_all_marker_name) +
  length(lipidomics_all_marker_name) +
  length(cytokine_all_marker_name) +
  length(clinical_test_all_marker_name) +
  length(gut_microbiome_all_marker_name) +
  length(skin_microbiome_all_marker_name) +
  length(oral_microbiome_all_marker_name) +
  length(nasal_microbiome_all_marker_name)

sum(transcriptome_variable_info$permutated_p_value < 0.05) +
  sum(proteomics_variable_info$permutated_p_value < 0.05) +
  sum(metabolomics_variable_info$permutated_p_value < 0.05) +
  sum(cytokine_variable_info$permutated_p_value < 0.05) +
  sum(clinical_test_variable_info$permutated_p_value < 0.05) +
  sum(gut_microbiome_variable_info$permutated_p_value < 0.05) +
  sum(skin_microbiome_variable_info$permutated_p_value < 0.05) +
  sum(oral_microbiome_variable_info$permutated_p_value < 0.05) +
  sum(nasal_microbiome_variable_info$permutated_p_value < 0.05)

length(transcriptome_all_marker_name_permutation) +
  length(proteomics_all_marker_name_permutation) +
  length(metabolomics_all_marker_name_permutation) +
  length(lipidomics_all_marker_name_permutation) +
  length(cytokine_all_marker_name_permutation) +
  length(clinical_test_all_marker_name_permutation) +
  length(gut_microbiome_all_marker_name_permutation) +
  length(skin_microbiome_all_marker_name_permutation) +
  length(oral_microbiome_all_marker_name_permutation) +
  length(nasal_microbiome_all_marker_name_permutation)

nrow(transcriptome_variable_info) +
  nrow(proteomics_variable_info) +
  nrow(metabolomics_variable_info) +
  nrow(lipidomics_variable_info) +
  nrow(cytokine_variable_info) +
  nrow(clinical_test_variable_info) +
  nrow(gut_microbiome_variable_info) +
  nrow(skin_microbiome_variable_info) +
  nrow(oral_microbiome_variable_info) +
  nrow(nasal_microbiome_variable_info)

c(7106 / 8556,
  240 / 302,
  711 / 814,
  22 / 66,
  25 / 46,
  430 / 846,
  127 / 158,
  195 / 205,
  98 / 104,
  206 / 208)

######volcano plot
temp_data <-
  rbind(
    transcriptome_variable_info[, c("variable_id",
                                    "mol_name",
                                    "cor_p",
                                    "cor_p_adjust",
                                    "spearman_cor",
                                    "coefficient")] %>%
      dplyr::mutate(data_type = "transcriptome"),
    proteomics_variable_info[, c("variable_id",
                                 "mol_name",
                                 "cor_p",
                                 "cor_p_adjust",
                                 "spearman_cor",
                                 "coefficient")] %>%
      dplyr::mutate(data_type = "proteomics"),
    metabolomics_variable_info[, c("variable_id",
                                   "mol_name",
                                   "cor_p",
                                   "cor_p_adjust",
                                   "spearman_cor",
                                   "coefficient")] %>%
      dplyr::mutate(data_type = "metabolomics"),
    cytokine_variable_info[, c("variable_id",
                               "mol_name",
                               "cor_p",
                               "cor_p_adjust",
                               "spearman_cor",
                               "coefficient")] %>%
      dplyr::mutate(data_type = "cytokine"),
    clinical_test_variable_info[, c("variable_id",
                                    "mol_name",
                                    "cor_p",
                                    "cor_p_adjust",
                                    "spearman_cor",
                                    "coefficient")] %>%
      dplyr::mutate(data_type = "clinical_test"),
    lipidomics_variable_info[, c("variable_id",
                                 "mol_name",
                                 "cor_p",
                                 "cor_p_adjust",
                                 "spearman_cor",
                                 "coefficient")] %>%
      dplyr::mutate(data_type = "lipidomics"),
    gut_microbiome_variable_info[, c("variable_id",
                                     "mol_name",
                                     "cor_p",
                                     "cor_p_adjust",
                                     "spearman_cor",
                                     "coefficient")] %>%
      dplyr::mutate(data_type = "gut_microbiome"),
    skin_microbiome_variable_info[, c("variable_id",
                                      "mol_name",
                                      "cor_p",
                                      "cor_p_adjust",
                                      "spearman_cor",
                                      "coefficient")] %>%
      dplyr::mutate(data_type = "skin_microbiome"),
    oral_microbiome_variable_info[, c("variable_id",
                                      "mol_name",
                                      "cor_p",
                                      "cor_p_adjust",
                                      "spearman_cor",
                                      "coefficient")] %>%
      dplyr::mutate(data_type = "oral_microbiome"),
    nasal_microbiome_variable_info[, c("variable_id",
                                       "mol_name",
                                       "cor_p",
                                       "cor_p_adjust",
                                       "spearman_cor",
                                       "coefficient")] %>%
      dplyr::mutate(data_type = "nasal_microbiome")
  )

# dplyr::mutate(data_type = factor(
#   data_type,
#   levels = c(
#     "transcriptome",
#     "proteomics",
#     "metabolomics",
#     "cytokine",
#     "clinical_test",
#     "lipidomics",
#     "gut_microbiome",
#     "skin_microbiome",
#     "oral_microbiome",
#     "nasal_microbiome"
#   )
# ))

plot <-
  temp_data %>%
  dplyr::mutate(data_type = factor(data_type, levels = names(omics_color))) %>%
  ggplot(aes(x = spearman_cor, coefficient)) +
  geom_point(aes(color = data_type), show.legend = FALSE) +
  theme_base +
  labs(x = "Spearman correlation",
       y = "lm beta") +
  scale_color_manual(values = omics_color) +
  geom_smooth(color = "black",
              method = "lm",
              show.legend = FALSE) +
  facet_wrap(facets = vars(data_type),
             nrow = 2,
             scales = "free")

plot

# ggsave(plot, filename = "linear_spearman.pdf", width = 14, height = 5)


temp_data %>%
  dplyr::mutate(data_type = factor(data_type, levels = names(omics_color))) %>%
  dplyr::group_by(data_type) %>%
  dplyr::summarise(cor = cor(spearman_cor, coefficient))

temp_data %>%
  dplyr::mutate(data_type = factor(data_type, levels = names(omics_color))) %>%
  dplyr::group_by(data_type) %>%
  dplyr::summarise(n = sum(cor_p < 0.05))

###volcano plot
# top10_up_marker_name <-
#   temp_data %>%
#   dplyr::filter(cor_p < 0.05 & spearman_cor > 0) %>%
#   dplyr::arrange(desc(spearman_cor)) %>%
#   head(10) %>%
#   dplyr::pull(variable_id)
#
# top10_down_marker_name <-
#   variable_info %>%
#   dplyr::filter(lm_p < 0.05 & coefficient < 0) %>%
#   dplyr::arrange(desc(abs(coefficient))) %>%
#   head(10) %>%
#   dplyr::pull(SYMBOL)

volcano_plot <-
  temp_data %>%
  mutate(
    marker = case_when(
      cor_p < 0.05 & spearman_cor > 0 ~ data_type,
      cor_p < 0.05 &
        spearman_cor < 0 ~ data_type,
      TRUE ~ "No"
    )
  ) %>%
  ggplot(aes(spearman_cor, -log(cor_p, 10))) +
  geom_point(aes(size = -log(cor_p, 10),
                 color = marker),
             alpha = 1) +
  theme_base +
  scale_color_manual(values = c(omics_color, "No" = "grey")) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.05, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(p-values)")
# ggrepel::geom_text_repel(aes(label = ifelse(
#   SYMBOL %in% top10_up_marker_name,
#   SYMBOL, NA
# )), size = 3) +
# ggrepel::geom_text_repel(aes(label = ifelse(
#   SYMBOL %in% top10_down_marker_name,
#   SYMBOL, NA
# )), size = 3)

volcano_plot

# ggsave(volcano_plot,
#        filename = "volcano_plot.pdf",
#        width = 6,
#        height = 7)


library(plyr)
top10_up_marker_name <-
  temp_data %>%
  plyr::dlply(.variables = .(data_type)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::filter(cor_p_adjust < 0.2 & spearman_cor > 0) %>%
      dplyr::arrange(desc(spearman_cor)) %>%
      head(10) %>%
      dplyr::pull(variable_id)
  }) %>%
  unlist()

top10_down_marker_name <-
  temp_data %>%
  plyr::dlply(.variables = .(data_type)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::filter(cor_p_adjust < 0.2 & spearman_cor < 0) %>%
      dplyr::arrange(desc(spearman_cor)) %>%
      head(10) %>%
      dplyr::pull(variable_id)
  }) %>%
  unlist()


volcano_plot <-
  temp_data %>%
  dplyr::mutate(data_type = factor(data_type, levels = names(omics_color))) %>%
  mutate(
    marker = case_when(
      cor_p_adjust < 0.2 & spearman_cor > 0 ~ "Up",
      cor_p_adjust < 0.2 &
        spearman_cor < 0 ~ "Down",
      TRUE ~ "No"
    )
  ) %>%
  ggplot(aes(spearman_cor, -log(cor_p_adjust, 10))) +
  geom_point(aes(size = -log(cor_p_adjust, 10),
                 color = marker),
             alpha = 0.7) +
  theme_base +
  scale_color_manual(values = marker_color) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log(0.2, 10), linetype = 2) +
  labs(x = "Spearman correlation", y = "-log10(FDR p-values)") +
  facet_wrap(facets = vars(data_type),
             nrow = 2,
             scales = "free") +
  theme(legend.position = "bottom") +
  ggrepel::geom_text_repel(aes(label = ifelse(
    variable_id %in% top10_up_marker_name,
    mol_name, NA
  )), size = 3) +
  ggrepel::geom_text_repel(aes(label = ifelse(
    variable_id %in% top10_down_marker_name,
    mol_name, NA
  )), size = 3)

volcano_plot

# ggsave(volcano_plot,
#        filename = "volcano_plot_omics.pdf",
#        width = 14,
#        height = 6)

# library(ggpubr)
# labs <- paste0(df$group, " (", df$value, "%)")
# ggpubr::ggdonutchart(
#   df,
#   "value",
#   label = labs,
#   lab.pos = "in",
#   lab.font = "white",
#   fill = "group",
#   color = "white",
#   palette = c("#00AFBB", "#E7B800", "#FC4E07")
# )





###transcriptome
temp1 <-
  c(transcriptome_variable_info$variable_id[transcriptome_variable_info$cor_p < 0.05],
    transcriptome_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(transcriptome_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["transcriptome"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "transcriptome_donut.pdf",
#        width = 7,
#        height = 7)

# library(ggVennDiagram)
# library(venneuler)
# library(VennDiagram)
# ggVennDiagram::ggVennDiagram(x = list("Linear" = transcriptome_variable_info$variable_id[transcriptome_variable_info$cor_p < 0.05],
#                                       "Changed" = transcriptome_all_marker_name))
library(VennDiagram)
venn.diagram(
  x = list("Linear" = transcriptome_variable_info$variable_id[transcriptome_variable_info$cor_p < 0.05],
           "Changed" = transcriptome_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()

###proteomics
temp1 <-
  c(proteomics_variable_info$variable_id[proteomics_variable_info$cor_p < 0.05],
    proteomics_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(proteomics_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["proteomics"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "proteomics_donut.pdf",
#        width = 7,
#        height = 7)


library(VennDiagram)
venn.diagram(
  x = list("Linear" = proteomics_variable_info$variable_id[proteomics_variable_info$cor_p < 0.05],
           "Changed" = proteomics_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()




###metabolomics
temp1 <-
  c(metabolomics_variable_info$variable_id[metabolomics_variable_info$cor_p < 0.05],
    metabolomics_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(metabolomics_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["metabolomics"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "metabolomics_donut.pdf",
#        width = 7,
#        height = 7)

library(VennDiagram)
venn.diagram(
  x = list("Linear" = metabolomics_variable_info$variable_id[metabolomics_variable_info$cor_p < 0.05],
           "Changed" = metabolomics_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()



###cytokine
temp1 <-
  c(cytokine_variable_info$variable_id[cytokine_variable_info$cor_p < 0.05],
    cytokine_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(cytokine_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["cytokine"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "cytokine_donut.pdf",
#        width = 7,
#        height = 7)

library(VennDiagram)
venn.diagram(
  x = list("Linear" = cytokine_variable_info$variable_id[cytokine_variable_info$cor_p < 0.05],
           "Changed" = cytokine_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()


###clinical_test
temp1 <-
  c(clinical_test_variable_info$variable_id[clinical_test_variable_info$cor_p < 0.05],
    clinical_test_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(clinical_test_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["clinical_test"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "clinical_test_donut.pdf",
#        width = 7,
#        height = 7)

library(VennDiagram)
venn.diagram(
  x = list("Linear" = clinical_test_variable_info$variable_id[clinical_test_variable_info$cor_p < 0.05],
           "Changed" = clinical_test_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()


###lipidomics
temp1 <-
  c(lipidomics_variable_info$variable_id[lipidomics_variable_info$cor_p < 0.05],
    lipidomics_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(lipidomics_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["lipidomics"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "lipidomics_donut.pdf",
#        width = 7,
#        height = 7)


library(VennDiagram)
venn.diagram(
  x = list("Linear" = lipidomics_variable_info$variable_id[lipidomics_variable_info$cor_p < 0.05],
           "Changed" = lipidomics_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()

###gut_microbiome
temp1 <-
  c(gut_microbiome_variable_info$variable_id[gut_microbiome_variable_info$cor_p < 0.05],
    gut_microbiome_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(gut_microbiome_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["gut_microbiome"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "gut_microbiome_donut.pdf",
#        width = 7,
#        height = 7)

library(VennDiagram)
venn.diagram(
  x = list("Linear" = gut_microbiome_variable_info$variable_id[gut_microbiome_variable_info$cor_p < 0.05],
           "Changed" = gut_microbiome_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()


###skin_microbiome
temp1 <-
  c(skin_microbiome_variable_info$variable_id[skin_microbiome_variable_info$cor_p < 0.05],
    skin_microbiome_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(skin_microbiome_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["skin_microbiome"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "skin_microbiome_donut.pdf",
#        width = 7,
#        height = 7)

library(VennDiagram)
venn.diagram(
  x = list("Linear" = skin_microbiome_variable_info$variable_id[skin_microbiome_variable_info$cor_p < 0.05],
           "Changed" = skin_microbiome_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()


###oral_microbiome
temp1 <-
  c(oral_microbiome_variable_info$variable_id[oral_microbiome_variable_info$cor_p < 0.05],
    oral_microbiome_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(oral_microbiome_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["oral_microbiome"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "oral_microbiome_donut.pdf",
#        width = 7,
#        height = 7)


library(VennDiagram)
venn.diagram(
  x = list("Linear" = oral_microbiome_variable_info$variable_id[oral_microbiome_variable_info$cor_p < 0.05],
           "Changed" = oral_microbiome_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()



###nasal_microbiome
temp1 <-
  c(nasal_microbiome_variable_info$variable_id[nasal_microbiome_variable_info$cor_p < 0.05],
    nasal_microbiome_all_marker_name) %>%
  unique() %>%
  length()

temp2 <-
  sum(nasal_microbiome_variable_info$cor_p < 0.05)

temp_data <-
  data.frame(group = c("Other",
                       "Linear"),
             value = c(temp1 - temp2,
                       temp2)) %>%
  dplyr::mutate(value = round(value * 100 / sum(value), 2))


plot <-
  ggpubr::ggdonutchart(
    temp_data,
    "value",
    label = paste0(temp_data$group, " (", temp_data$value, "%)"),
    lab.pos = "in",
    lab.font = "white",
    fill = "group",
    color = "white",
    palette = c(unname(omics_color["nasal_microbiome"]), "grey")
  )
plot
# ggsave(plot,
#        filename = "nasal_microbiome_donut.pdf",
#        width = 7,
#        height = 7)

library(VennDiagram)
venn.diagram(
  x = list("Linear" = nasal_microbiome_variable_info$variable_id[nasal_microbiome_variable_info$cor_p < 0.05],
           "Changed" = nasal_microbiome_all_marker_name),
  category.names = c("Linear" , "Changed"),
  filename = NULL,
  fill = c("red", "darkblue")
) %>%
  grid.draw()
