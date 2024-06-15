###
no_function()

tinyTools::setwd_project()
library(tidyverse)
rm(list = ls())
####load raw data
load("data/from_xin/Revision_MultiOmes_0509.RData")
setwd("data_analysis/proteome/data_preparation")
expression_data = 
swathprot.PCR.df

colnames(expression_data)

sample_info = 
  expression_data %>% 
  dplyr::select(SampleID, SubjectID:CL4)

expression_data =
expression_data %>% 
  dplyr::select(-c(SampleID, SubjectID:CL4)) %>% 
  t() %>%
  as.data.frame()

colnames(expression_data) = sample_info$SampleID

variable_info =
  data.frame(variable_id = rownames(expression_data))

dim(expression_data)
dim(sample_info)
dim(variable_info)

colnames(expression_data) == sample_info$SampleID
rownames(expression_data)  == variable_info$variable_id

sample_info = 
sample_info %>% 
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID)


###add information to variable_info
variable_info$variable_id

library(clusterProfiler)
library(org.Hs.eg.db)
library(plyr)
library(UniprotR)


UNIPROT = 
  clusterProfiler::bitr(
    geneID = variable_info$variable_id,
    fromType = "SYMBOL",
    toType = "UNIPROT",
    OrgDb = org.Hs.eg.db, 
    drop = FALSE
  )
    
# protein_info =
#   UNIPROT %>%
#   plyr::dlply(.variables = .(SYMBOL)) %>%
#   purrr::map(function(x) {
#     if (nrow(x) == 1) {
#       return(x)
#     }
#     temp =
#       UniprotR::GetNamesTaxa(ProteinAccList = x$UNIPROT)
#     
#     temp
#   })
# 
# 
# save(protein_info, file = "protein_info")
load("protein_info")



UNIPROT = 
  UNIPROT %>% 
  plyr::dlply(.variables = .(SYMBOL))

UNIPROT2 = 
purrr::map2(UNIPROT, protein_info,
            function(x, y) {
              if (nrow(x) == 1) {
                return(x)
              }
              
              entry = stringr::str_replace(y$Entry.name, "_HUMAN", "")
              gene = y$Gene.names %>% 
                stringr::str_split(" ") %>% 
                purrr::map(function(z){
                  z[1]
                }) %>% 
                unlist()
              
              temp_idx = 
              data.frame(entry, gene) %>% 
                apply(1, function(z){
                  sum(z == x$SYMBOL[1])
                })
              temp_idx[is.na(temp_idx)] = 0
              
              final_idx = 
              which(temp_idx == temp_idx[which.max(temp_idx)])
              
              
              if(length(final_idx) > 0){
                return(x[final_idx,])
              }else{
                  return(x)
                }
            })


idx = 
UNIPROT2 %>% lapply(nrow) %>% unlist() %>% `>`(1) %>% which()


UNIPROT2[[idx[2]]]

variable_info %>% 
  dplyr::left_join(UNIPROT, by = c("variable_id" = "SYMBOL"))




save(sample_info, file = "sample_info")
save(expression_data, file = "expression_data")
save(variable_info, file = "variable_info")

