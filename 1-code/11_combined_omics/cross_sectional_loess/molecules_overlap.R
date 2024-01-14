no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load spearman linear markers
###load("data)
{
  load(
    "3-data_analysis/plasma_transcriptome/linear_spearman/cross_section/variable_info"
  )
  
  transcriptome_variable_info <-
    variable_info %>%
    dplyr::mutate(mol_name = SYMBOL)
  
  load(
    "3-data_analysis/plasma_proteomics/linear_spearman/cross_section/variable_info"
  )
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
  
  load(
    "3-data_analysis/plasma_lipidomics/linear_spearman/cross_section/variable_info"
  )
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
}

{
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
}

####load crest marker
load("3-data_analysis/combined_omics/DE_SWAN/temp_data")

temp_data <-
  temp_data %>%
  dplyr::filter(p_value_adjust < 0.05)

dir.create("3-data_analysis/combined_omics/molecules_crest_linear_overlap",
           recursive = TRUE)

setwd("3-data_analysis/combined_omics/molecules_crest_linear_overlap")

###transcriptomics
transcriptome_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "transcriptome")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "transcriptome_", ""))

transcriptome_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "transcriptome")) %>%
  dplyr::filter(center == 61) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "transcriptome_", ""))

transcriptome_linear <-
  transcriptome_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = transcriptome_crest1$variable_id,
    crest2 = transcriptome_crest2$variable_id,
    linear = transcriptome_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

transcriptome_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["transcriptome"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["transcriptome"])))

library(ggplotify)
transcriptome_upset_plot <-
  as.ggplot(transcriptome_upset_plot)

transcriptome_upset_plot
save(transcriptome_upset_plot, file = "transcriptome_upset_plot")
ggsave(
  transcriptome_upset_plot,
  file = "transcriptome_upset.pdf",
  width = 7,
  height = 5
)


###proteomics
proteomics_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "proteomics")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "proteomics_", ""))

proteomics_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "proteomics")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "proteomics_", ""))

proteomics_linear <-
  proteomics_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = proteomics_crest1$variable_id,
    crest2 = proteomics_crest2$variable_id,
    linear = proteomics_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

proteomics_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["proteomics"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["proteomics"])))

library(ggplotify)
proteomics_upset_plot <-
  as.ggplot(proteomics_upset_plot)

proteomics_upset_plot
save(proteomics_upset_plot, file = "proteomics_upset_plot")
ggsave(
  proteomics_upset_plot,
  file = "proteomics_upset.pdf",
  width = 7,
  height = 5
)





###metabolomics
metabolomics_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "metabolomics")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "metabolomics_", ""))

metabolomics_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "metabolomics")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "metabolomics_", ""))

metabolomics_linear <-
  metabolomics_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = metabolomics_crest1$variable_id,
    crest2 = metabolomics_crest2$variable_id,
    linear = metabolomics_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

metabolomics_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["metabolomics"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["metabolomics"])))

library(ggplotify)
metabolomics_upset_plot <-
  as.ggplot(metabolomics_upset_plot)

metabolomics_upset_plot
save(metabolomics_upset_plot, file = "metabolomics_upset_plot")
ggsave(
  metabolomics_upset_plot,
  file = "metabolomics_upset.pdf",
  width = 7,
  height = 5
)








###cytokine
cytokine_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "cytokine")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "cytokine_", ""))

cytokine_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "cytokine")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "cytokine_", ""))

cytokine_linear <-
  cytokine_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = cytokine_crest1$variable_id,
    crest2 = cytokine_crest2$variable_id,
    linear = cytokine_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

cytokine_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["cytokine"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["cytokine"])))

library(ggplotify)
cytokine_upset_plot <-
  as.ggplot(cytokine_upset_plot)

cytokine_upset_plot
save(cytokine_upset_plot, file = "cytokine_upset_plot")
ggsave(
  cytokine_upset_plot,
  file = "cytokine_upset.pdf",
  width = 7,
  height = 5
)








###clinical_test
clinical_test_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "clinical_test")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "clinical_test_", ""))

clinical_test_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "clinical_test")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "clinical_test_", ""))

clinical_test_linear <-
  clinical_test_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = clinical_test_crest1$variable_id,
    crest2 = clinical_test_crest2$variable_id,
    linear = clinical_test_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

clinical_test_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["clinical_test"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["clinical_test"])))

library(ggplotify)
clinical_test_upset_plot <-
  as.ggplot(clinical_test_upset_plot)

clinical_test_upset_plot
save(clinical_test_upset_plot, file = "clinical_test_upset_plot")
ggsave(
  clinical_test_upset_plot,
  file = "clinical_test_upset.pdf",
  width = 7,
  height = 5
)








###lipidomics
lipidomics_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "lipidomics")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "lipidomics_", ""))

lipidomics_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "lipidomics")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "lipidomics_", ""))

lipidomics_linear <-
  lipidomics_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = lipidomics_crest1$variable_id,
    crest2 = lipidomics_crest2$variable_id,
    linear = lipidomics_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

lipidomics_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["lipidomics"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["lipidomics"])))

library(ggplotify)
lipidomics_upset_plot <-
  as.ggplot(lipidomics_upset_plot)

lipidomics_upset_plot
save(lipidomics_upset_plot, file = "lipidomics_upset_plot")
ggsave(
  lipidomics_upset_plot,
  file = "lipidomics_upset.pdf",
  width = 7,
  height = 5
)








###gut_microbiome
gut_microbiome_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "gut_microbiome")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "gut_microbiome_", ""))

gut_microbiome_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "gut_microbiome")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "gut_microbiome_", ""))

gut_microbiome_linear <-
  gut_microbiome_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = gut_microbiome_crest1$variable_id,
    crest2 = gut_microbiome_crest2$variable_id,
    linear = gut_microbiome_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

gut_microbiome_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["gut_microbiome"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["gut_microbiome"])))

library(ggplotify)
gut_microbiome_upset_plot <-
  as.ggplot(gut_microbiome_upset_plot)

gut_microbiome_upset_plot
save(gut_microbiome_upset_plot, file = "gut_microbiome_upset_plot")
ggsave(
  gut_microbiome_upset_plot,
  file = "gut_microbiome_upset.pdf",
  width = 7,
  height = 5
)








###skin_microbiome
skin_microbiome_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "skin_microbiome")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "skin_microbiome_", ""))

skin_microbiome_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "skin_microbiome")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "skin_microbiome_", ""))

skin_microbiome_linear <-
  skin_microbiome_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = skin_microbiome_crest1$variable_id,
    crest2 = skin_microbiome_crest2$variable_id,
    linear = skin_microbiome_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

skin_microbiome_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["skin_microbiome"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["skin_microbiome"])))

library(ggplotify)
skin_microbiome_upset_plot <-
  as.ggplot(skin_microbiome_upset_plot)

skin_microbiome_upset_plot
save(skin_microbiome_upset_plot, file = "skin_microbiome_upset_plot")
ggsave(
  skin_microbiome_upset_plot,
  file = "skin_microbiome_upset.pdf",
  width = 7,
  height = 5
)








###oral_microbiome
oral_microbiome_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "oral_microbiome")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "oral_microbiome_", ""))

oral_microbiome_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "oral_microbiome")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "oral_microbiome_", ""))

oral_microbiome_linear <-
  oral_microbiome_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = oral_microbiome_crest1$variable_id,
    crest2 = oral_microbiome_crest2$variable_id,
    linear = oral_microbiome_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

oral_microbiome_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["oral_microbiome"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["oral_microbiome"])))

library(ggplotify)
oral_microbiome_upset_plot <-
  as.ggplot(oral_microbiome_upset_plot)

oral_microbiome_upset_plot
save(oral_microbiome_upset_plot, file = "oral_microbiome_upset_plot")
ggsave(
  oral_microbiome_upset_plot,
  file = "oral_microbiome_upset.pdf",
  width = 7,
  height = 5
)







###nasal_microbiome
nasal_microbiome_crest1 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "nasal_microbiome")) %>%
  dplyr::filter(center == 44) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "nasal_microbiome_", ""))

nasal_microbiome_crest2 <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "nasal_microbiome")) %>%
  dplyr::filter(center == 57) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "nasal_microbiome_", ""))

nasal_microbiome_linear <-
  nasal_microbiome_variable_info %>%
  dplyr::filter(cor_p < 0.05)

library(ComplexHeatmap)

temp <-
  list(
    crest1 = nasal_microbiome_crest1$variable_id,
    crest2 = nasal_microbiome_crest2$variable_id,
    linear = nasal_microbiome_linear$variable_id
  )

m1 = make_comb_mat(temp)
m1

nasal_microbiome_upset_plot <-
  UpSet(m1,
        comb_col = omics_color["nasal_microbiome"],
        right_annotation = upset_right_annotation(m = m1,
                                                  gp = gpar(fill = omics_color["nasal_microbiome"])))

library(ggplotify)
nasal_microbiome_upset_plot <-
  as.ggplot(nasal_microbiome_upset_plot)

nasal_microbiome_upset_plot
save(nasal_microbiome_upset_plot, file = "nasal_microbiome_upset_plot")
ggsave(
  nasal_microbiome_upset_plot,
  file = "nasal_microbiome_upset.pdf",
  width = 7,
  height = 5
)


library(patchwork)

plot <-
transcriptome_upset_plot + proteomics_upset_plot + metabolomics_upset_plot + cytokine_upset_plot +
  clinical_test_upset_plot + lipidomics_upset_plot + gut_microbiome_upset_plot +
  skin_microbiome_upset_plot + oral_microbiome_upset_plot + nasal_microbiome_upset_plot +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = 'a')

ggsave(plot, filename = "upset_plot.pdf", width = 16, height = 6)

