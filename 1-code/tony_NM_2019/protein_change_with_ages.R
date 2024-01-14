no_source()

rm(list = ls())

setwd(r4projects::get_project_wd())

source("1-code/100-tools.R")

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

load("3-data_analysis/plasma_proteomics/data_preparation/object")

object_cross_section_loess <-
  object_cross_section_loess %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(age = as.numeric(sample_id))

protein1 <-
  readxl::read_xlsx(
    "3-data_analysis/tony_NM_2019/41591_2019_673_MOESM3_ESM.xlsx",
    sheet = "ST1 Nomenclature 2,925 proteins",
    skip = 2
  )

protein2 <-
  readxl::read_xlsx(
    "3-data_analysis/tony_NM_2019/41591_2019_673_MOESM3_ESM.xlsx",
    sheet = "ST2 Nomenclature 1,305 proteins",
    skip = 2
  )

dir.create("3-data_analysis/tony_NM_2019")
setwd("3-data_analysis/tony_NM_2019")

proteomics_variable_info <-
  object@variable_info

intersect_protein_id1 <-
  intersect(protein1$UniProt, proteomics_variable_info$UNIPROT)
intersect_protein_id2 <-
  intersect(protein2$UniProt, proteomics_variable_info$UNIPROT)

intersect_protein_id <-
  unique(c(intersect_protein_id1, intersect_protein_id2))

intersect_protein_id <-
  intersect_protein_id[!is.na(intersect_protein_id)]

data1 <-
  protein1 %>%
  dplyr::filter(UniProt %in% intersect_protein_id) %>%
  dplyr::select(ID, UniProt)

data2 <-
  protein2 %>%
  dplyr::filter(UniProt %in% intersect_protein_id) %>%
  dplyr::select(R_ID, UniProt) %>%
  dplyr::rename(ID = R_ID)

data <-
  rbind(data1, data2) %>%
  dplyr::arrange(UniProt, ID) %>%
  dplyr::distinct(UniProt, .keep_all = TRUE)

dir.create("age_value_plot")

for (i in 1:nrow(data)) {
  cat(i, " ")
  uniprot_id <-
    data$UniProt[i]
  index <-
    which(proteomics_variable_info$variable_id == data$ID[i])
  
  if (length(index) == 0) {
    next()
  }
  
  plot <-
    object %>%
    ggplot_mass_dataset(direction = "variable",
                        variable_index = index) +
    geom_point(aes(x = adjusted_age)) +
    geom_smooth(aes(x = adjusted_age), method = "loess") +
    base_theme +
    labs(x = "Age (years)",
         y = "Z-score",
         title = proteomics_variable_info$variable_id[index])
  
  ggsave(
    plot,
    filename = file.path(
      "age_value_plot",
      paste0(proteomics_variable_info$variable_id[index], ".pdf")
    ),
    width = 7,
    height = 7
  )
}
