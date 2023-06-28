no_source()

masstools::setwd_project()
rm(list = ls())
source("code/tools.R")

load("data_analysis/phenotype/data_preparation/phenotype_data")

load("data_analysis/plasma_cytokine/data_preparation/object")
cytokine_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

cytokine_object <-
  object

load("data_analysis/plasma_lipidomics/data_preparation/object")
lipidomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

lipidomics_object <-
  object

load("data_analysis/plasma_metabolomics/data_preparation/metabolite/object")
metabolomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

metabolomics_object <-
  object

load("data_analysis/plasma_proteomics/data_preparation/object")
proteomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

proteomics_object <-
  object

load("data_analysis/plasma_transcriptome/data_preparation/object")
transcriptome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

transcriptome_object <-
  object

load("data_analysis/clinical_test/data_preparation/object")
clinical_test_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

clinical_test_object <-
  object

load("data_analysis/gut_microbiome/data_preparation/object")
gut_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

gut_microbiome_object <-
  object

load("data_analysis/nasal_microbiome/data_preparation/object")
nasal_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

nasal_microbiome_object <-
  object

load("data_analysis/skin_microbiome/data_preparation/object")
skin_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

skin_microbiome_object <-
  object

load("data_analysis/oral_microbiome/data_preparation/object")
oral_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

oral_microbiome_object <-
  object

dir.create("data_analysis/summary_info")
setwd("data_analysis/summary_info")

dim(transcriptome_object)
dim(proteomics_object)
dim(metabolomics_object)
dim(cytokine_object)
dim(clinical_test_object)
dim(lipidomics_object)
dim(gut_microbiome_object)
unique(gut_microbiome_object@variable_info$Genus) %>%
  length

unique(skin_microbiome_object@variable_info$Genus) %>%
  length

unique(oral_microbiome_object@variable_info$Genus) %>%
  length

unique(nasal_microbiome_object@variable_info$Genus) %>%
  length

all_subject_id <-
  c(
    transcriptome_sample_info$subject_id,
    proteomics_sample_info$subject_id,
    metabolomics_sample_info$subject_id,
    cytokine_sample_info$subject_id,
    clinical_test_sample_info$subject_id,
    lipidomics_sample_info$subject_id,
    gut_microbiome_sample_info$subject_id,
    skin_microbiome_sample_info$subject_id,
    oral_microbiome_sample_info$subject_id,
    nasal_microbiome_sample_info$subject_id
  ) %>%
  unique()

library(plyr)

sample_number <-
  purrr::map2(.x = list(
    cytokine_sample_info,
    lipidomics_sample_info,
    metabolomics_sample_info,
    proteomics_sample_info,
    transcriptome_sample_info,
    clinical_test_sample_info,
    gut_microbiome_sample_info,
    nasal_microbiome_sample_info,
    skin_microbiome_sample_info,
    oral_microbiome_sample_info
  ),
  .y = c(
    "cytokine",
    "lipidomics",
    "metabolomics",
    "proteomics",
    "transcriptome",
    "clinical_test",
    "gut_microbiome",
    "nasal_microbiome",
    "skin_microbiome",
    "oral_microbiome"
  ),
  function(x, y) {
    number = x %>%
      dplyr::count(subject_id,
                   name = paste0(y, "_number"))
    
    days <-
      x %>%
      plyr::dlply(.variables = .(subject_id)) %>%
      purrr::map(function(z) {
        z <-
          z %>%
          arrange(collection_date)
        data.frame(
          subject_id = z$subject_id[1],
          day_range = as.numeric(z$collection_date[nrow(z)] - z$collection_date[1])
        )
      }) %>%
      dplyr::bind_rows() %>%
      as.data.frame()
    
    colnames(days)[2] <-
      paste0(y, "_day_range")
    
    number %>%
      dplyr::left_join(days, by = "subject_id")
  })

sample_number <-
  sample_number[[1]] %>%
  dplyr::full_join(sample_number[[2]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[3]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[4]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[5]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[6]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[7]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[8]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[9]], by = c("subject_id")) %>%
  dplyr::full_join(sample_number[[10]], by = c("subject_id"))

transcriptome_subject_info <-
  transcriptome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

proteomics_subject_info <-
  proteomics_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

metabolomics_subject_info <-
  metabolomics_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

cytokine_subject_info <-
  cytokine_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

lipidomics_subject_info <-
  lipidomics_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

clinical_test_subject_info <-
  clinical_test_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

gut_microbiome_subject_info <-
  gut_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

skin_microbiome_subject_info <-
  skin_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

oral_microbiome_subject_info <-
  oral_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

nasal_microbiome_subject_info <-
  nasal_microbiome_sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(
    subject_id,
    subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

subject_info <-
  transcriptome_subject_info %>%
  dplyr::full_join(proteomics_subject_info, by = "subject_id") %>%
  dplyr::full_join(metabolomics_subject_info, by = "subject_id") %>%
  dplyr::full_join(lipidomics_subject_info, by = "subject_id") %>%
  dplyr::full_join(clinical_test_subject_info, by = "subject_id") %>%
  dplyr::full_join(cytokine_subject_info, by = "subject_id") %>%
  dplyr::full_join(gut_microbiome_subject_info, by = "subject_id") %>%
  dplyr::full_join(oral_microbiome_subject_info, by = "subject_id") %>%
  dplyr::full_join(nasal_microbiome_subject_info, by = "subject_id") %>%
  dplyr::full_join(skin_microbiome_subject_info, by = "subject_id")

subject_id_random <-
  subject_info %>%
  dplyr::select(matches("subject_id_random")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

consented <-
  subject_info %>%
  dplyr::select(matches("consented")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

IRIS <-
  subject_info %>%
  dplyr::select(matches("IRIS")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

SSPG <-
  subject_info %>%
  dplyr::select(matches("SSPG")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

FPG <-
  subject_info %>%
  dplyr::select(matches("FPG")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

diabetes_class <-
  subject_info %>%
  dplyr::select(matches("diabetes_class")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

Gender <-
  subject_info %>%
  dplyr::select(matches("Gender")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

Ethnicity <-
  subject_info %>%
  dplyr::select(matches("Ethnicity")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  })

adjusted_age <-
  subject_info %>%
  dplyr::select(matches("adjusted_age")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

BMI <-
  subject_info %>%
  dplyr::select(matches("BMI")) %>%
  apply(1, function(x) {
    x <-
      as.character(x)
    x <- x[!is.na(x)]
    x[1]
  }) %>%
  as.numeric()

subject_info <-
  data.frame(
    subject_id = subject_info$subject_id,
    subject_id_random = subject_id_random,
    consented,
    IRIS,
    SSPG,
    FPG,
    diabetes_class,
    Gender,
    Ethnicity,
    adjusted_age,
    BMI
  )

subject_info <-
  subject_info %>%
  dplyr::left_join(sample_number, by = "subject_id")


subject_info$subject_id_random <-
  masstools::name_duplicated(subject_info$subject_id_random)

#####circos plot
library(circlize)

df <-
  data.frame(
    factors = subject_info$subject_id_random,
    x = 1,
    y = 1,
    subject_info,
    stringsAsFactors = TRUE
  ) %>%
  dplyr::arrange(adjusted_age) %>%
  dplyr::mutate(factors = factor(factors, levels = factors))

circos.par(
  "track.height" = 0.2,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0, nrow(df) - 1), 90),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.5, 1.5))

##age
range(df$adjusted_age, na.rm = TRUE)
temp_value <- df$adjusted_age

circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = c(0.8 * min(temp_value), 1.1 * max(temp_value, na.rm = TRUE)),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = c(0.8 * min(temp_value),
             round((
               min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)
             ) / 2, 2),
             round(max(
               temp_value, na.rm = TRUE
             ), 2)),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.lines(
      x = mean(xlim, na.rm = TRUE),
      y =  temp_value[i],
      pch = 16,
      cex = 8,
      type = "h",
      col = ggsci::pal_aaas()(n = 10)[4],
      lwd = 2
    )
    
    #plot country labels
    circos.text(
      x = 1,
      y = 105,
      labels = name,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.5
      # adj = aa
    )
    
    circos.points(
      x = mean(xlim),
      y =  temp_value[i],
      pch = 16,
      cex = 0.8,
      col = ggsci::pal_aaas()(n = 10)[4]
    )
  }
)

##BMI
range(df$BMI, na.rm = TRUE)
temp_value <- df$BMI

circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = c(
    0.8 * min(temp_value, na.rm = TRUE),
    1.1 * max(temp_value, na.rm = TRUE)
  ),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = c(
        0.8 * min(temp_value, na.rm = TRUE),
        round((
          min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)
        ) / 2, 2),
        round(max(temp_value, na.rm = TRUE), 2)
      ),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.lines(
      x = mean(xlim, na.rm = TRUE),
      y =  temp_value[i],
      pch = 16,
      cex = 8,
      type = "h",
      col = ggsci::pal_tron()(n = 10)[1],
      lwd = 2
    )
    
    circos.points(
      x = mean(xlim),
      y =  temp_value[i],
      pch = 16,
      cex = 0.8,
      col = ggsci::pal_tron()(n = 10)[1]
    )
  }
)

## sex
temp_sex <- df$Gender
temp_sex[is.na(temp_sex)] <- "grey"
temp_sex[temp_sex == "Female"] <- sex_color["Female"]
temp_sex[temp_sex == "Male"] <- sex_color["Male"]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_sex[i],
      bg.border = "black"
    )
  }
)

## Ethnicity
temp_ethnicity <- df$Ethnicity
temp_ethnicity[is.na(temp_ethnicity)] <- "grey"
temp_ethnicity[temp_ethnicity == "Caucasian"] <-
  ethnicity_color["Caucasian"]
temp_ethnicity[temp_ethnicity == "Asian"] <-
  ethnicity_color["Asian"]
temp_ethnicity[temp_ethnicity == "Hispanics"] <-
  ethnicity_color["Hispanics"]
temp_ethnicity[temp_ethnicity == "Black"] <-
  ethnicity_color["Black"]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    #plot country labels
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_ethnicity[i],
      bg.border = "black"
    )
  }
)

## IRIS
temp_iris <- df$IRIS
temp_iris[is.na(temp_iris)] <- "grey"
temp_iris[temp_iris == "IR"] <-
  iris_color["IR"]
temp_iris[temp_iris == "Unknown"] <-
  iris_color["Unknown"]
temp_iris[temp_iris == "IS"] <-
  iris_color["IS"]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    #plot country labels
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_iris[i],
      bg.border = "black"
    )
  }
)


###age
#####age
age <-
  df$adjusted_age
library(gghalves)
plot_age <-
  age %>%
  data.frame(class = "class", value = .) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(
    binaxis = "y",
    color = ggsci::pal_aaas()(n = 10)[4],
    fill = ggsci::pal_aaas()(n = 10)[4],
    shape = 16,
    binwidth = 1,
    stackdir = "center"
  ) +
  theme_bw() +
  labs(x = "", y = "") +
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
plot_age

ggsave(plot_age,
       filename = "plot_age.pdf",
       width = 3,
       height = 10)

###BMI

bmi <-
  df$BMI

plot_bmi <-
  bmi %>%
  data.frame(class = "class", value = .) %>%
  ggplot(aes(x = class, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(
    binaxis = "y",
    color = ggsci::pal_tron()(n = 10)[1],
    fill = ggsci::pal_tron()(n = 10)[1],
    shape = 16,
    binwidth = 0.6,
    stackdir = "center"
  ) +
  theme_bw() +
  labs(x = "", y = "") +
  # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
plot_bmi

ggsave(plot_bmi,
       filename = "plot_bmi.pdf",
       width = 3,
       height = 10)

##sex

sex <-
  df$Gender

plot_sex <-
  sex %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("Female", "Male"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = sex_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_sex

# ggsave(plot_sex,
#        filename = "plot_sex.pdf",
#        width = 1.5,
#        height = 10)


##ethnicity
ethnicity <-
  df$Ethnicity

plot_ethnicity <-
  ethnicity %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c(
    "Asian", "Black", "Caucasian", "Hispanics"
  ))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = ethnicity_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_ethnicity

# ggsave(
#   plot_ethnicity,
#   filename = "plot_ethnicity.pdf",
#   width = 1.5,
#   height = 10
# )


##IRIS
iris <-
  df$IRIS

plot_iris <-
  iris %>%
  data.frame(class = "class", value = .) %>%
  dplyr::mutate(value = factor(value, levels = c("IR", "IS", "Unknown"))) %>%
  ggplot(aes(x = class)) +
  geom_bar(
    aes(fill = value),
    color = "black",
    position = "stack",
    show.legend = FALSE,
    width = 2
  ) +
  scale_fill_manual(values = iris_color) +
  theme_bw() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_iris

# ggsave(plot_iris,
#        filename = "plot_iris.pdf",
#        width = 1.5,
#        height = 10)

library(plyr)

temp_data <-
  subject_info %>%
  dplyr::arrange(desc(adjusted_age)) %>%
  dplyr::mutate(subject_id_random = factor(subject_id_random, subject_id_random)) %>%
  dplyr::select(subject_id_random, matches("_number")) %>%
  dplyr::rename_all(
    .funs = function(x) {
      stringr::str_replace(x, "_number", "")
    }
  ) %>%
  tidyr::pivot_longer(
    cols = -subject_id_random,
    names_to = "class",
    values_to = "number"
  )

###max sample number
max_sample_number <-
  temp_data %>%
  plyr::dlply(.variables = .(subject_id_random)) %>%
  purrr::map(function(x) {
    sum(x$number, na.rm = TRUE)
  }) %>%
  unlist() %>%
  max()

temp_data$number[is.na(temp_data$number)] <- 0

temp_data <-
  temp_data %>%
  plyr::dlply(.variables = .(subject_id_random)) %>%
  purrr::map(function(x) {
    rbind(
      x,
      data.frame(
        subject_id_random = x$subject_id_random[1],
        class = "non",
        number = max_sample_number - sum(x$number)
      )
    )
  }) %>%
  dplyr::bind_rows()

temp_data <-
  temp_data %>%
  dplyr::mutate(class = factor(
    class,
    levels = c(
      "transcriptome",
      "proteomics",
      "metabolomics",
      "cytokine",
      "clinical_test",
      "lipidomics",
      "gut_microbiome",
      "skin_microbiome",
      "oral_microbiome",
      "nasal_microbiome",
      "non"
    ) %>% rev()
  ))

subject_sample_number <-
  temp_data %>%
  plyr::dlply(.variables = .(subject_id_random)) %>%
  purrr::map(function(x) {
    sum(x$number[x$class != "non"])
  }) %>%
  unlist()

temp_data2 <-
  subject_info %>%
  dplyr::select(subject_id_random, matches("_day_range")) %>%
  dplyr::rename_all(
    .funs = function(x) {
      stringr::str_replace(x, "_day_range", "")
    }
  ) %>%
  tidyr::pivot_longer(
    cols = -subject_id_random,
    names_to = "class",
    values_to = "days"
  ) %>%
  plyr::dlply(.variables = .(subject_id_random)) %>%
  purrr::map(function(x) {
    data.frame(
      subject_id_random = x$subject_id_random[1],
      days = max(x$days, na.rm = TRUE)
    )
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(subject_id_random = factor(subject_id_random, levels = levels(temp_data$subject_id_random))) %>%
  dplyr::mutate(scaled_days = days * max_sample_number / max(days))

plot <-
  temp_data %>%
  ggplot(aes(y = subject_id_random, x = number)) +
  geom_bar(aes(fill = class), stat = "identity") +
  geom_line(
    aes(x = scaled_days,
        y = subject_id_random),
    group = 1,
    data = temp_data2,
    color = "red"
  ) +
  geom_point(
    aes(x = scaled_days,
        y = subject_id_random),
    shape = 16,
    color = "red",
    data = temp_data2
  ) +
  theme_base +
  theme(
    axis.text.y = element_text(size = 5),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c(omics_color, "non" = "grey")) +
  labs(y = "") +
  scale_y_discrete(label = subject_sample_number) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0)),
    sec.axis = dup_axis(
      name = "Collection range (days)",
      breaks = c(c(0, 500, 1000, 1500, 2000) * max_sample_number /
                   max(2471)),
      labels = c(0, 500, 1000, 1500, 2000)
    )
  )

plot

# ggsave(plot,
#        filename = "sample_number.pdf",
#        width = 4,
#        height = 8)


dim(transcriptome_sample_info)
dim(transcriptome_subject_info)

dim(proteomics_sample_info)
dim(proteomics_subject_info)

dim(metabolomics_sample_info)
dim(metabolomics_subject_info)

dim(cytokine_sample_info)
dim(cytokine_subject_info)

dim(clinical_test_sample_info)
dim(clinical_test_subject_info)

dim(lipidomics_sample_info)
dim(lipidomics_subject_info)

dim(gut_microbiome_sample_info)
dim(gut_microbiome_subject_info)

dim(skin_microbiome_sample_info)
dim(skin_microbiome_subject_info)

dim(oral_microbiome_sample_info)
dim(oral_microbiome_subject_info)

dim(nasal_microbiome_sample_info)
dim(nasal_microbiome_subject_info)


###age and sex distributation
plot <-
  df %>%
  ggplot() +
  geom_histogram(aes(x = adjusted_age, fill = Gender),
                 binwidth = 2,
                 color = "black") +
  scale_fill_manual(values = sex_color) +
  base_theme +
  labs(x = "Age (yreas)", y = "Count") +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1))
plot
ggsave(plot,
       filename = "age_sex.pdf",
       width = 7,
       height = 7)

###age and ethnicity distributation
plot <-
  df %>%
  ggplot() +
  geom_histogram(aes(x = adjusted_age, fill = Ethnicity),
                 binwidth = 2,
                 color = "black") +
  scale_fill_manual(values = ethnicity_color) +
  base_theme +
  labs(x = "Age (yreas)", y = "Count") +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1))

plot
ggsave(plot,
       filename = "age_ethnicity.pdf",
       width = 7,
       height = 7)

###age and IRIS distributation
plot <-
  df %>%
  ggplot() +
  geom_histogram(aes(x = adjusted_age, fill = IRIS),
                 binwidth = 2,
                 color = "black") +
  scale_fill_manual(values = iris_color) +
  base_theme +
  labs(x = "Age (yreas)", y = "Count") +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1))

plot
ggsave(plot,
       filename = "age_iris.pdf",
       width = 7,
       height = 7)


###age and BMI distributation
plot <-
  df %>%
  ggplot(aes(adjusted_age, BMI)) +
  geom_point(aes(adjusted_age, BMI),
             size = 5) +
  base_theme +
  geom_smooth(method = "lm", color = "red") +
  labs(x = "Age (yreas)", y = "BMI") +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1))

plot

ggsave(plot,
       filename = "age_bmi.pdf",
       width = 7,
       height = 7)
