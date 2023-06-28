no_source()

rm(list = ls())

setwd(masstools::get_project_wd())

source("code/tools.R")

load("data_analysis/plasma_proteomics/data_preparation/object")

protein1 <-
  readxl::read_xlsx(
    "data_analysis/tony_NM_2019/41591_2019_673_MOESM3_ESM.xlsx",
    sheet = "ST1 Nomenclature 2,925 proteins",
    skip = 2
  )

protein2 <-
  readxl::read_xlsx(
    "data_analysis/tony_NM_2019/41591_2019_673_MOESM3_ESM.xlsx",
    sheet = "ST2 Nomenclature 1,305 proteins",
    skip = 2
  )

tony_protein <-
  readxl::read_xlsx("data_analysis/tony_NM/protein_waves.xlsx")

tony_protein <-
  tony_protein %>%
  dplyr::rename(variable_id = variable) %>%
  dplyr::left_join(protein1[, c("ID", "UniProt", "EntrezGeneID", "EntrezGeneSymbol")],
                   by = c("variable_id" = "ID"))

variable_info <-
  object@variable_info

######overwhole interaction
ipop_protein_id <-
  variable_info$UNIPROT[!is.na(variable_info$UNIPROT)] %>%
  unique()

tony_protein_id <-
  tony_protein$UniProt[!is.na(tony_protein$UniProt)] %>%
  stringr::str_replace_all(" ", ",") %>%
  stringr::str_split(",") %>%
  unlist() %>%
  stringr::str_trim(side = "both")

tony_protein_id <-
  tony_protein_id[tony_protein_id != ""]

length(ipop_protein_id)
length(tony_protein_id)

intersect_protein_id <-
  intersect(ipop_protein_id,
            tony_protein_id)

intersect_protein_id

length(intersect_protein_id)

load("data_analysis/combined_omics/DE_SWAN/temp_data")

dir.create(
  "data_analysis/plasma_proteomics/protein_waves/compared_with_tony_paper",
  recursive = TRUE
)
setwd("data_analysis/plasma_proteomics/protein_waves/compared_with_tony_paper")

ipop_protein <-
  temp_data %>%
  dplyr::filter(stringr::str_detect(variable_id, "proteomics")) %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "proteomics_", ""))

variable_info <-
  object@variable_info

ipop_protein <-
  ipop_protein %>%
  dplyr::left_join(variable_info, by = "variable_id")

###crest 1
ipop_crest1 <-
  ipop_protein %>%
  dplyr::filter(center >= 41 & center <= 45)

tony_crest1 <-
  tony_protein %>%
  dplyr::filter(qvalue.34 < 0.05) %>%
  dplyr::select(variable_id, UniProt, EntrezGeneID, EntrezGeneSymbol)

# dim(ipop_crest1)
# dim(tony_crest1)

ipop_protein_id1 <-
  ipop_crest1$UNIPROT[!is.na(ipop_crest1$UNIPROT)]

tony_protein_id1 <-
  tony_crest1$UniProt[!is.na(tony_crest1$UniProt)] %>%
  stringr::str_replace_all(" ", ",") %>%
  stringr::str_split(",") %>%
  unlist() %>%
  stringr::str_trim(side = "both")

tony_protein_id1 <-
  tony_protein_id1[tony_protein_id1 != ""]

# length(ipop_protein_id1)
# length(tony_protein_id1)

intersect(ipop_protein_id1,
          tony_protein_id1) %>%
  length()

###crest 2
ipop_crest2 <-
  ipop_protein %>%
  dplyr::filter(center >= 54 & center <= 60)
# dplyr::filter(p_value_adjust < 0.05)

tony_crest2 <-
  tony_protein %>%
  dplyr::filter(qvalue.60 < 0.05) %>%
  dplyr::select(variable_id, UniProt, EntrezGeneID, EntrezGeneSymbol)

# dim(ipop_crest2)
# dim(tony_crest2)

ipop_protein_id2 <-
  ipop_crest2$UNIPROT[!is.na(ipop_crest2$UNIPROT)]

tony_protein_id2 <-
  tony_crest2$UniProt[!is.na(tony_crest2$UniProt)] %>%
  stringr::str_replace_all(" ", ",") %>%
  stringr::str_split(",") %>%
  unlist() %>%
  stringr::str_trim(side = "both")

tony_protein_id2 <-
  tony_protein_id2[tony_protein_id2 != ""]

# length(ipop_protein_id2)
# length(tony_protein_id2)

intersect(ipop_protein_id2,
          tony_protein_id2) %>%
  length()

intersect(ipop_protein_id,
          tony_protein_id) %>%
  length()

intersect(ipop_protein_id1,
          tony_protein_id1) %>%
  length()

intersect(ipop_protein_id2,
          tony_protein_id2) %>%
  length()


length(tony_protein_id1)
length(tony_protein_id2)
intersect(tony_protein_id1,
          tony_protein_id2)

library(ggVennDiagram)
ggVennDiagram(x = list(tony_protein_id1,
                       tony_protein_id2))
