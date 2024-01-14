no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)
library(microbiomedataset)

load("3-data_analysis/plasma_transcriptome/data_preparation/sample_info")

subject_info <-
  sample_info %>%
  dplyr::select(subject_id, adjusted_age) %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::distinct(subject_id, .keep_all = TRUE)

###load("data)
load("3-data_analysis/plasma_transcriptome/linear_spearman/cross_section/pca_object")
transcriptome_pca_object <-
  pca_object

load("3-data_analysis/plasma_proteomics/linear_spearman/cross_section/pca_object")
proteomics_pca_object <-
  pca_object

load("3-data_analysis/plasma_metabolomics/linear_spearman/cross_section/pca_object")
metabolomics_pca_object <-
  pca_object

load("3-data_analysis/plasma_cytokine/linear_spearman/cross_section/pca_object")
cytokine_pca_object <-
  pca_object

load("3-data_analysis/clinical_test/linear_spearman/cross_section/pca_object")
clinical_test_pca_object <-
  pca_object

load("3-data_analysis/plasma_lipidomics/linear_spearman/cross_section/pca_object")
lipidomics_pca_object <-
  pca_object

load("3-data_analysis/gut_microbiome/linear_spearman/cross_section/pca_object")
gut_microbiome_pca_object <-
  pca_object

load("3-data_analysis/skin_microbiome/linear_spearman/cross_section/pca_object")
skin_microbiome_pca_object <-
  pca_object

load("3-data_analysis/oral_microbiome/linear_spearman/cross_section/pca_object")
oral_microbiome_pca_object <-
  pca_object

load("3-data_analysis/nasal_microbiome/linear_spearman/cross_section/pca_object")
nasal_microbiome_pca_object <-
  pca_object

dir.create("3-data_analysis/combined_omics/PCA_marker")
setwd("3-data_analysis/combined_omics/PCA_marker")

temp_data <-
  rbind(
    data.frame(
      subject_id = rownames(transcriptome_pca_object$x),
      PC = apply(transcriptome_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "transcriptome"
    ),
    data.frame(
      subject_id = rownames(proteomics_pca_object$x),
      PC = apply(proteomics_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "proteomics"
    ),
    data.frame(
      subject_id = rownames(metabolomics_pca_object$x),
      PC = apply(metabolomics_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "metabolomics"
    ),
    data.frame(
      subject_id = rownames(cytokine_pca_object$x),
      PC = apply(cytokine_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "cytokine"
    ),
    data.frame(
      subject_id = rownames(clinical_test_pca_object$x),
      PC = apply(clinical_test_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "clinical_test"
    ),
    data.frame(
      subject_id = rownames(lipidomics_pca_object$x),
      PC = apply(lipidomics_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "lipidomics"
    ),
    data.frame(
      subject_id = rownames(gut_microbiome_pca_object$x),
      PC = apply(gut_microbiome_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "gut_microbiome"
    ),
    data.frame(
      subject_id = rownames(skin_microbiome_pca_object$x),
      PC = apply(skin_microbiome_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "skin_microbiome"
    ),
    data.frame(
      subject_id = rownames(oral_microbiome_pca_object$x),
      PC = apply(oral_microbiome_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "oral_microbiome"
    ),
    data.frame(
      subject_id = rownames(nasal_microbiome_pca_object$x),
      PC = apply(nasal_microbiome_pca_object$x[, 1, drop = FALSE], 1, mean),
      class = "nasal_microbiome"
    )
  ) %>%
  dplyr::left_join(subject_info, by = "subject_id") %>%
  dplyr::filter(!is.na(adjusted_age)) %>%
  dplyr::mutate(class = factor(class, levels = names(omics_color)))

plot <-
  temp_data %>%
  ggplot(aes(adjusted_age, PC)) +
  geom_point(aes(color = class), show.legend = FALSE) +
  theme_base +
  geom_smooth(aes(color = class), method = "lm",
              show.legend = FALSE) +
  facet_wrap(facets = vars(class),
             scales = "free",
             nrow = 2) +
  scale_color_manual(values = omics_color) +
  labs(x = "Age (years)", y = "PC1")
# scale_y_continuous(limits = c(-10, 15))
plot
# ggsave(plot,
#        filename = "pca_age.pdf",
#        width = 14,
#        height = 6)

temp_data %>%
  group_by(class) %>%
  dplyr::summarise(cor = cor(PC, adjusted_age, method = "spearman"))

library(plyr)
temp_data %>%
  plyr::dlply(.variables = .(class)) %>%
  purrr::map(function(x) {
    temp <- cor.test(x$adjusted_age, x$PC, method = "spearman")
    data.frame(cor = unname(abs(temp$estimate)),
               p = unname(abs(temp$p.value)),
               class = x$class[1])
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


####transcriptome vs proteomics
transcriptome_pc1 <-
  temp_data %>%
  dplyr::filter(class == "transcriptome") %>%
  dplyr::arrange(subject_id) %>% 
  dplyr::select(subject_id,
                transcriptome = PC)

proteomics_pc1 <-
  temp_data %>%
  dplyr::filter(class == "proteomics") %>%
  dplyr::arrange(subject_id) %>% 
  dplyr::select(subject_id,
                proteomics = PC)

metabolomics_pc1 <-
  temp_data %>%
  dplyr::filter(class == "metabolomics") %>%
  dplyr::arrange(subject_id) %>% 
  dplyr::select(subject_id,
                metabolomics = PC)

temp <-
transcriptome_pc1 %>% 
  dplyr::full_join(proteomics_pc1,
                   by = "subject_id") %>% 
  dplyr::full_join(metabolomics_pc1,
                   by = "subject_id")

cor.test(temp$transcriptome,
         temp$proteomics, method = "spearman")

cor.test(temp$transcriptome,
         temp$metabolomics, method = "spearman")

cor.test(temp$proteomics,
         temp$metabolomics, method = "spearman")
