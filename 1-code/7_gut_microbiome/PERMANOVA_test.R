no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3-data_analysis/gut_microbiome/data_preparation/object_cross_section")

dir.create("3-data_analysis/gut_microbiome/linear_spearman/cross_section/",
           recursive = TRUE)

setwd("3-data_analysis/gut_microbiome/linear_spearman/cross_section")

object_cross_section

dim(object_cross_section)


####only remain the genus level
library(microbiomedataset)

object_cross_section <-
  object_cross_section %>%
  summarize_variables(what = "sum_intensity", group_by = "Genus")

##only remain the genus at least 10% samples
dim(object_cross_section)

non_zero_per <-
  apply(object_cross_section, 1, function(x) {
    sum(x != 0) / ncol(object_cross_section)
  })

idx <-
  which(non_zero_per > 0.1)

object_cross_section <-
  object_cross_section[idx, ]


object_cross_section <-
  object_cross_section %>%
  transform2relative_intensity()

###find marker which are change according to aging

###linear regression
library(tidyverse)
library(ggpubr)
library(rstatix)

expression_data <-
  extract_expression_data(object_cross_section) %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

sample_info <-
  object_cross_section@sample_info

#######adjust BMI, sex, and IRIS, ethnicity
expression_data <-
  lm_adjust(expression_data = expression_data,
            sample_info = sample_info,
            threads = 3)

temp_object <- object_cross_section
temp_object@expression_data <- expression_data

####Permanova test
library(vegan)
temp_data <-
  t(expression_data)

# Run PERMANOVA test
permanova_result <- 
  adonis2(temp_data ~ sample_info$adjusted_age, 
          permutations = 999)

permanova_result

save(permanova_result, file = "permanova_result")




####linear regression data
dim(temp_data)

temp_data2 <-
  data.frame(age = sample_info$adjusted_age,
             temp_data)


lm_result <- lm(age ~ ., data = temp_data2)

summary(lm_result)$r.squared

save(lm_result, file = "lm_result")






####PLS
library(mixOmics)

pls_result <-
  pls(X = temp_data,
      Y = sample_info$adjusted_age,
      ncomp = 10) # assuming 10 components, but adjust as needed
summary(pls_result)


# Evaluate the performance of the PLS model
pls_performance <-
  perf(pls_result, validation = "Mfold", folds = 7)

save(pls_performance, file = "pls_performance")

mean(pls_performance$measures$R2$summary$mean)

