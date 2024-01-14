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
load(
  "3-data_analysis/plasma_transcriptome/linear_spearman/cross_section/permanova_result"
)
transcriptome_permanova_result <-
  permanova_result

load(
  "3-data_analysis/plasma_proteomics/linear_spearman/cross_section/permanova_result"
)
proteomics_permanova_result <-
  permanova_result

load(
  "3-data_analysis/plasma_metabolomics/linear_spearman/cross_section/permanova_result"
)
metabolomics_permanova_result <-
  permanova_result

load("3-data_analysis/plasma_cytokine/linear_spearman/cross_section/permanova_result")
cytokine_permanova_result <-
  permanova_result

load("3-data_analysis/clinical_test/linear_spearman/cross_section/permanova_result")
clinical_test_permanova_result <-
  permanova_result

load(
  "3-data_analysis/plasma_lipidomics/linear_spearman/cross_section/permanova_result"
)
lipidomics_permanova_result <-
  permanova_result

load("3-data_analysis/gut_microbiome/linear_spearman/cross_section/permanova_result")
gut_microbiome_permanova_result <-
  permanova_result

load("3-data_analysis/skin_microbiome/linear_spearman/cross_section/permanova_result")
skin_microbiome_permanova_result <-
  permanova_result

load("3-data_analysis/oral_microbiome/linear_spearman/cross_section/permanova_result")
oral_microbiome_permanova_result <-
  permanova_result

load(
  "3-data_analysis/nasal_microbiome/linear_spearman/cross_section/permanova_result"
)
nasal_microbiome_permanova_result <-
  permanova_result


dir.create("3-data_analysis/combined_omics/age_linear_regression_r2")
setwd("3-data_analysis/combined_omics/age_linear_regression_r2")

transcriptome_permanova_result$R2
proteomics_permanova_result$R2
metabolomics_permanova_result$R2
cytokine_permanova_result$R2
clinical_test_permanova_result$R2
lipidomics_permanova_result$R2
gut_microbiome_permanova_result$R2
skin_microbiome_permanova_result$R2
oral_microbiome_permanova_result$R2
nasal_microbiome_permanova_result$R2


# ###load("data)
# load("3-data_analysis/plasma_transcriptome/linear_spearman/cross_section/lm_result")
# transcriptome_lm_result <-
#   lm_result
#
# load("3-data_analysis/plasma_proteomics/linear_spearman/cross_section/lm_result")
# proteomics_lm_result <-
#   lm_result
#
# load("3-data_analysis/plasma_metabolomics/linear_spearman/cross_section/lm_result")
# metabolomics_lm_result <-
#   lm_result
#
# load("3-data_analysis/plasma_cytokine/linear_spearman/cross_section/lm_result")
# cytokine_lm_result <-
#   lm_result
#
# load("3-data_analysis/clinical_test/linear_spearman/cross_section/lm_result")
# clinical_test_lm_result <-
#   lm_result
#
# load("3-data_analysis/plasma_lipidomics/linear_spearman/cross_section/lm_result")
# lipidomics_lm_result <-
#   lm_result
#
# load("3-data_analysis/gut_microbiome/linear_spearman/cross_section/lm_result")
# gut_microbiome_lm_result <-
#   lm_result
#
# load("3-data_analysis/skin_microbiome/linear_spearman/cross_section/lm_result")
# skin_microbiome_lm_result <-
#   lm_result
#
# load("3-data_analysis/oral_microbiome/linear_spearman/cross_section/lm_result")
# oral_microbiome_lm_result <-
#   lm_result
#
# load("3-data_analysis/nasal_microbiome/linear_spearman/cross_section/lm_result")
# nasal_microbiome_lm_result <-
#   lm_result
#
# dir.create("3-data_analysis/combined_omics/age_linear_regression_r2")
# setwd("3-data_analysis/combined_omics/age_linear_regression_r2")
#
# temp_data <-
#   rbind(
#     data.frame(
#       r2 = summary(transcriptome_lm_result)$r.squared,
#       class = "transcriptome"
#     ),
#     data.frame(
#       r2 = summary(proteomics_lm_result)$r.squared,
#       class = "proteomics"
#     ),
#     data.frame(
#       r2 = summary(metabolomics_lm_result)$r.squared,
#       class = "metabolomics"
#     ),
#     data.frame(
#       r2 = summary(cytokine_lm_result)$r.squared,
#       class = "cytokine"
#     ),
#     data.frame(
#       r2 = summary(clinical_test_lm_result)$r.squared,
#       class = "clinical_test"
#     ),
#     data.frame(
#       r2 = summary(lipidomics_lm_result)$r.squared,
#       class = "lipidomics"
#     ),
#     data.frame(
#       r2 = summary(gut_microbiome_lm_result)$r.squared,
#       class = "gut_microbiome"
#     ),
#     data.frame(
#       r2 = summary(skin_microbiome_lm_result)$r.squared,
#       class = "skin_microbiome"
#     ),
#     data.frame(
#       r2 = summary(oral_microbiome_lm_result)$r.squared,
#       class = "oral_microbiome"
#     ),
#     data.frame(
#       r2 = summary(nasal_microbiome_lm_result)$r.squared,
#       class = "nasal_microbiome"
#     )
#   )
#
# plot <-
#   temp_data %>%
#   ggplot(aes(adjusted_age, PC)) +
#   geom_point(aes(color = class), show.legend = FALSE) +
#   theme_base +
#   geom_smooth(aes(color = class), method = "lm",
#               show.legend = FALSE) +
#   facet_wrap(facets = vars(class),
#              scales = "free",
#              nrow = 2) +
#   scale_color_manual(values = omics_color) +
#   labs(x = "Age (years)", y = "PC1")
# # scale_y_continuous(limits = c(-10, 15))
# plot
# # ggsave(plot,
# #        filename = "pca_age.pdf",
# #        width = 14,
# #        height = 6)
#
# temp_data %>%
#   group_by(class) %>%
#   dplyr::summarise(cor = cor(PC, adjusted_age, method = "spearman"))
#
# library(plyr)
# temp_data %>%
#   plyr::dlply(.variables = .(class)) %>%
#   purrr::map(function(x) {
#     temp <- cor.test(x$adjusted_age, x$PC, method = "spearman")
#     data.frame(cor = unname(abs(temp$estimate)),
#                p = unname(abs(temp$p.value)),
#                class = x$class[1])
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()

setwd(masstools::get_project_wd())
###load("data)
load(
  "3-data_analysis/plasma_transcriptome/linear_spearman/cross_section/pls_performance"
)
transcriptome_pls_performance <-
  pls_performance

load(
  "3-data_analysis/plasma_proteomics/linear_spearman/cross_section/pls_performance"
)
proteomics_pls_performance <-
  pls_performance

load(
  "3-data_analysis/plasma_metabolomics/linear_spearman/cross_section/pls_performance"
)
metabolomics_pls_performance <-
  pls_performance

load("3-data_analysis/plasma_cytokine/linear_spearman/cross_section/pls_performance")
cytokine_pls_performance <-
  pls_performance

load("3-data_analysis/clinical_test/linear_spearman/cross_section/pls_performance")
clinical_test_pls_performance <-
  pls_performance

load(
  "3-data_analysis/plasma_lipidomics/linear_spearman/cross_section/pls_performance"
)
lipidomics_pls_performance <-
  pls_performance

load("3-data_analysis/gut_microbiome/linear_spearman/cross_section/pls_performance")
gut_microbiome_pls_performance <-
  pls_performance

load("3-data_analysis/skin_microbiome/linear_spearman/cross_section/pls_performance")
skin_microbiome_pls_performance <-
  pls_performance

load("3-data_analysis/oral_microbiome/linear_spearman/cross_section/pls_performance")
oral_microbiome_pls_performance <-
  pls_performance

load("3-data_analysis/nasal_microbiome/linear_spearman/cross_section/pls_performance")
nasal_microbiome_pls_performance <-
  pls_performance

dir.create("3-data_analysis/combined_omics/age_linear_regression_r2")
setwd("3-data_analysis/combined_omics/age_linear_regression_r2")

mean(transcriptome_pls_performance$measures$R2$values$value)
mean(proteomics_pls_performance$measures$R2$values$value)
mean(metabolomics_pls_performance$measures$R2$values$value)
mean(cytokine_pls_performance$measures$R2$values$value)
mean(clinical_test_pls_performance$measures$R2$values$value)
mean(lipidomics_pls_performance$measures$R2$values$value)

mean(gut_microbiome_pls_performance$measures$R2$values$value)
mean(skin_microbiome_pls_performance$measures$R2$values$value)
mean(oral_microbiome_pls_performance$measures$R2$values$value)
mean(nasal_microbiome_pls_performance$measures$R2$values$value)
