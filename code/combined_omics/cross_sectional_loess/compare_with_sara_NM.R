no_source()

rm(list = ls())

masstools::setwd_project()

sara_result <-
  readxl::read_xlsx(
    "data_analysis/sara_NM/41591_2019_719_MOESM2_ESM.xlsx",
    sheet = 1,
    skip = 3
  )

head(sara_result)

load(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/cor_data"
)

load(
  "data_analysis/combined_omics/fuzzy_c_means_clustering/cross_section_loess/final_cluster_info"
)

load(
  "data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

dir.create("data_analysis/sara_NM/compare")
setwd("data_analysis/sara_NM/compare")

colnames(sara_result)[1] <- "variable_id"

colnames(sara_result)[2:4] <- c("rho1", "p_value1", "fdr1")
colnames(sara_result)[5:7] <- c("rho2", "p_value2", "fdr2")

plot(sara_result$rho1, sara_result$rho2)

cor(sara_result$rho1, sara_result$rho2, use = "complete.obs")

sum(sara_result$fdr1 < 0.1, na.rm = TRUE)
sum(sara_result$fdr2 < 0.1, na.rm = TRUE)

sum(sara_result$fdr1 < 0.2, na.rm = TRUE)
sum(sara_result$fdr2 < 0.2, na.rm = TRUE)

cor_data$fdr <- p.adjust(cor_data$p_value, method = "fdr")

sara_result %>%
  dplyr::filter(fdr2 < 0.1) %>%
  dplyr::arrange(Molecule_Type) %>%
  dplyr::filter(!Molecule_Type %in% c("Metabolites")) %>%
  dplyr::arrange(variable_id)

cor_data %>%
  dplyr::filter(
    variable_id %in% c(
      "clinical_test_EGFR",
      "clinical_test_MONOAB",
      "clinical_test_MONO",
      "clinical_test_RDW",
      "clinical_test_A1C",
      "cytokine_PDGFBB",
      "cytokine_MCP3",
      "cytokine_MIG",
      "cytokine_VEGFD",
      "cytokine_IP10",
      "cytokine_GMCSF"
    )
  ) %>%
  dplyr::arrange(variable_id)

top10_pos <-
  sara_result %>%
  # dplyr::filter(fdr2 < 0.1 & rho2 > 0) %>%
  dplyr::arrange(Molecule_Type) %>%
  dplyr::filter(!Molecule_Type %in% c("Metabolites", 'Gut_MicrobialGenes', "Gut_Microbes")) %>%
  dplyr::arrange(variable_id) %>%
  dplyr::arrange(dplyr::desc(rho2)) %>%
  head(12)

top10_neg <-
  sara_result %>%
  # dplyr::filter(fdr2 < 0.2 & rho2 < 0) %>%
  dplyr::arrange(Molecule_Type) %>%
  dplyr::filter(
    !Molecule_Type %in% c(
      "Metabolites",
      'Gut_MicrobialGenes',
      "Gut_Microbes",
      "Nasal_Microbes"
    )
  ) %>%
  dplyr::arrange(variable_id) %>%
  dplyr::arrange(rho2) %>%
  head(20)

write.csv(top10_pos, "top10_pos.csv", row.names = FALSE)
write.csv(top10_neg, "top10_neg.csv", row.names = FALSE)


top10_pos

variable_info <-
  object_cross_section_loess@variable_info

final_cluster_info %>%
  dplyr::filter(
    variable_id %in% c(
      "clinical_test_A1C",
      "clinical_test_MONO",
      "clinical_test_MONOAB",
      "clinical_test_RDW",
      "cytokine_PDGFBB",
      "cytokine_MIG",
      "transcriptome_LIPM",
      "transcriptome_NOB1",
      "transcriptome_CD14",
      "transcriptome_ACAA2"
    )
  ) %>%
  dplyr::left_join(variable_info[, c("variable_id", "mol_name",
                                     "GENENAME",
                                     "test_name")],
                   by = "variable_id")

top10_neg %>%
  dplyr::arrange(Molecule_Type)


final_cluster_info %>%
  dplyr::filter(
    variable_id %in% c(
      "clinical_test_EGFR",
      "clinical_test_ALB",
      "proteomics_IGF2R",
      "proteomics_LV657",
      "transcriptome_POU5F1B",
      "transcriptome_LINC01011",
      "transcriptome_ORC3",
      "transcriptome_PPAPDC2",
      "transcriptome_WWC2-AS2",
      "transcriptome_LINC00936",
      "transcriptome_PTPRK",
      "transcriptome_PPP1CA",
      "transcriptome_FLJ32255",
      "transcriptome_PCP2",
      "transcriptome_BRICD5",
      "transcriptome_CAB39",
      "transcriptome_TCTEX1D4",
      "transcriptome_DCLRE1A"
    )
  ) %>%
  dplyr::left_join(variable_info[, c("variable_id", "mol_name",
                                     "GENENAME",
                                     "test_name")],
                   by = "variable_id")



up_10 <-
  cor_data %>%
  dplyr::filter(fdr < 0.05 & correlation > 0) %>%
  dplyr::arrange(desc(correlation)) %>%
  head(10)

down_10 <-
  cor_data %>%
  dplyr::filter(fdr < 0.05 & correlation < 0) %>%
  dplyr::arrange(correlation) %>%
  head(10)


final_cluster_info %>%
  dplyr::filter(variable_id %in% up_10$variable_id) %>%
  pull(cluster)

final_cluster_info %>%
  dplyr::filter(variable_id %in% down_10$variable_id) %>%
  pull(cluster)



