no_source()

rm(list = ls())
setwd(masstools::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load("3-data_analysis/clinical_test/data_preparation/object_cross_section")

dir.create("3-data_analysis/clinical_test/linear_spearman/cross_section/",
           recursive = TRUE)

setwd("3-data_analysis/clinical_test/linear_spearman/cross_section")

object_cross_section

dim(object_cross_section)

sum(is.na(object_cross_section))
sum(object_cross_section < 0)

show_variable_missing_values(object_cross_section, percentage = TRUE)

####remove variables has lots of NA
object_cross_section <-
  object_cross_section %>%
  mutate_variable_na_freq()

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(na_freq < 0.5)

show_sample_missing_values(object_cross_section, percentage = TRUE)

####remove samples has lots of NA
object_cross_section <-
  object_cross_section %>%
  mutate_sample_na_freq()

object_cross_section <-
  object_cross_section %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(na_freq < 0.5)

object_cross_section <-
  object_cross_section %>%
  impute_mv(method = "knn")

###find marker which are change according to aging
###linear regression
library(tidyverse)
library(ggpubr)
library(rstatix)

expression_data <-
  extract_expression_data(object_cross_section) %>%
  apply(1, function(x) {
    (x - mean(x)) /sd(x)
  }) %>%
  t() %>%
  as.data.frame()

library(plyr)

sample_info <-
  object_cross_section@sample_info

plot(sample_info$adjusted_age, as.numeric(expression_data[16,]))
cor.test(sample_info$adjusted_age, as.numeric(expression_data[16,]), method = "spearman")

#######adjust BMI, sex, and IRIS, ethnicity
expression_data <-
  lm_adjust(expression_data = expression_data,
            sample_info = sample_info,
            threads = 3)

plot(sample_info$adjusted_age, as.numeric(expression_data[16,]))
cor.test(sample_info$adjusted_age, as.numeric(expression_data[16,]), method = "spearman")


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

