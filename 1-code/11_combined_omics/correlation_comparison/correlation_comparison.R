no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
###transcriptome
load(
  "3-data_analysis/plasma_transcriptome/correlation_with_age/cross_section/cor_data"
)
load(
  "3-data_analysis/plasma_transcriptome/correlation_with_age/cross_section_loess/cor_data_loess"
)
transcriptome_cor_data <-
  data.frame(cor_data, class = "transcriptome")
transcriptome_cor_data_loess <-
  data.frame(cor_data_loess, class = "transcriptome")

###proteomics
load("3-data_analysis/plasma_proteomics/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/plasma_proteomics/correlation_with_age/cross_section_loess/cor_data_loess"
)
proteomics_cor_data <-
  data.frame(cor_data, class = "proteomics")
proteomics_cor_data_loess <-
  data.frame(cor_data_loess, class = "proteomics")

###metabolomics
load(
  "3-data_analysis/plasma_metabolomics/correlation_with_age/cross_section/cor_data"
)
load(
  "3-data_analysis/plasma_metabolomics/correlation_with_age/cross_section_loess/cor_data_loess"
)
metabolomics_cor_data <-
  data.frame(cor_data, class = "metabolomics")
metabolomics_cor_data_loess <-
  data.frame(cor_data_loess, class = "metabolomics")

###cytokine
load("3-data_analysis/plasma_cytokine/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/plasma_cytokine/correlation_with_age/cross_section_loess/cor_data_loess"
)
cytokine_cor_data <-
  data.frame(cor_data, class = "cytokine")
cytokine_cor_data_loess <-
  data.frame(cor_data_loess, class = "cytokine")

###clinical_test
load("3-data_analysis/clinical_test/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/clinical_test/correlation_with_age/cross_section_loess/cor_data_loess"
)
clinical_test_cor_data <-
  data.frame(cor_data, class = "clinical_test")
clinical_test_cor_data_loess <-
  data.frame(cor_data_loess, class = "clinical_test")


###lipidomics
load("3-data_analysis/plasma_lipidomics/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/plasma_lipidomics/correlation_with_age/cross_section_loess/cor_data_loess"
)
lipidomics_cor_data <-
  data.frame(cor_data, class = "lipidomics")
lipidomics_cor_data_loess <-
  data.frame(cor_data_loess, class = "lipidomics")

###gut_microbiome
load("3-data_analysis/gut_microbiome/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/gut_microbiome/correlation_with_age/cross_section_loess/cor_data_loess"
)
gut_microbiome_cor_data <-
  data.frame(cor_data, class = "gut_microbiome")
gut_microbiome_cor_data_loess <-
  data.frame(cor_data_loess, class = "gut_microbiome")

###skin_microbiome
load("3-data_analysis/skin_microbiome/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/skin_microbiome/correlation_with_age/cross_section_loess/cor_data_loess"
)
skin_microbiome_cor_data <-
  data.frame(cor_data, class = "skin_microbiome")
skin_microbiome_cor_data_loess <-
  data.frame(cor_data_loess, class = "skin_microbiome")

###oral_microbiome
load("3-data_analysis/oral_microbiome/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/oral_microbiome/correlation_with_age/cross_section_loess/cor_data_loess"
)
oral_microbiome_cor_data <-
  data.frame(cor_data, class = "oral_microbiome")
oral_microbiome_cor_data_loess <-
  data.frame(cor_data_loess, class = "oral_microbiome")


###nasal_microbiome
load("3-data_analysis/nasal_microbiome/correlation_with_age/cross_section/cor_data")
load(
  "3-data_analysis/nasal_microbiome/correlation_with_age/cross_section_loess/cor_data_loess"
)
nasal_microbiome_cor_data <-
  data.frame(cor_data, class = "nasal_microbiome")
nasal_microbiome_cor_data_loess <-
  data.frame(cor_data_loess, class = "nasal_microbiome")


dir.create("3-data_analysis/combined_omics/correlation_comparison",
           recursive = TRUE)

setwd("3-data_analysis/combined_omics/correlation_comparison")

temp_data <-
  rbind(
    transcriptome_cor_data,
    proteomics_cor_data,
    metabolomics_cor_data,
    cytokine_cor_data,
    clinical_test_cor_data,
    lipidomics_cor_data,
    gut_microbiome_cor_data,
    skin_microbiome_cor_data,
    oral_microbiome_cor_data,
    nasal_microbiome_cor_data
  )

temp_data_loess <-
  rbind(
    transcriptome_cor_data_loess,
    proteomics_cor_data_loess,
    metabolomics_cor_data_loess,
    cytokine_cor_data_loess,
    clinical_test_cor_data_loess,
    lipidomics_cor_data_loess,
    gut_microbiome_cor_data_loess,
    skin_microbiome_cor_data_loess,
    oral_microbiome_cor_data_loess,
    nasal_microbiome_cor_data_loess
  )


temp_data <-
  data.frame(
    variable_id = temp_data$variable_id,
    correlation = temp_data$correlation,
    p_value = temp_data$p_value,
    correlation_loess = temp_data_loess$correlation,
    p_value_loess = temp_data_loess$p_value,
    class = temp_data$class
  )

plot <-
  temp_data %>%
  ggplot(aes(x = correlation, correlation_loess)) +
  geom_point(aes(color = class), alpha = 0.7) +
  geom_smooth(aes(color = class), method = "lm") +
  geom_smooth(color = "black", method = "lm") +
  scale_color_manual(values = omics_color) +
  base_theme +
  theme(legend.position = "bottom") +
  labs(x = "Correlation (raw data)",
       y = "Correlation (loess data)")

plot 

cor.test(temp_data$correlation, temp_data$correlation_loess)

# ggsave(plot, filename = "correlation_comparison.pdf", width = 7, height = 7)




