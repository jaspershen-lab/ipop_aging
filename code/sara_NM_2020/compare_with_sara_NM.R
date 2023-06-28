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

marker_fdr0.2 <-
  sara_result %>%
  dplyr::filter(fdr2 < 0.2)


openxlsx::write.xlsx(
  marker_fdr0.2,
  file = "marker_fdr0.2.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)
