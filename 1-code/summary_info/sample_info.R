no_source()

setwd(r4projects::get_project_wd())

load(
  "3-data_analysis/combined_omics/data_preparation/cross_section_loess/object_cross_section_loess"
)

load("3-data_analysis/plasma_cytokine/data_preparation/object")
cytokine_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/plasma_lipidomics/data_preparation/object")
lipidomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/plasma_metabolomics/data_preparation/metabolite/object")
metabolomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/plasma_proteomics/data_preparation/object")
proteomics_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/plasma_transcriptome/data_preparation/object")
transcriptome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/clinical_test/data_preparation/object")
clinical_test_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/gut_microbiome/data_preparation/object")
gut_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/nasal_microbiome/data_preparation/object")
nasal_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/skin_microbiome/data_preparation/object")
skin_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

load("3-data_analysis/oral_microbiome/data_preparation/object")
oral_microbiome_sample_info <-
  object@sample_info %>%
  dplyr::filter(!is.na(adjusted_age))

age_index <-
  data.frame(from = c(25, 40, 45, 50, 55, 60, 65),
             to =   c(40, 45, 50, 55, 60, 65, 75))

list(
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
) %>%
  purrr::map(function(sample_info) {
    apply(age_index, 1, function(x) {
      sum(sample_info$adjusted_age > x[1] &
            sample_info$adjusted_age <= x[2])
    })
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

list(
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
) %>%
  purrr::map(function(sample_info) {
    range(sample_info$adjusted_age)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

c(
  cytokine_sample_info$subject_id,
  lipidomics_sample_info$subject_id,
  metabolomics_sample_info$subject_id,
  proteomics_sample_info$subject_id,
  transcriptome_sample_info$subject_id,
  clinical_test_sample_info$subject_id,
  gut_microbiome_sample_info$subject_id,
  nasal_microbiome_sample_info$subject_id,
  skin_microbiome_sample_info$subject_id,
  oral_microbiome_sample_info$subject_id
) %>%
  unique()

dir.create("3-data_analysis/sample_info")
setwd("3-data_analysis/sample_info")

temp_data <-
  rbind(
    transcriptome_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "blood"),
    proteomics_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "blood"),
    metabolomics_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "blood"),
    cytokine_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "blood"),
    clinical_test_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "clinical_test"),
    lipidomics_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "blood"),
    gut_microbiome_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "gut_microbiome"),
    skin_microbiome_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "skin_microbiome"),
    oral_microbiome_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "oral_microbiome"),
    nasal_microbiome_sample_info %>%
      dplyr::select(subject_id_random, collection_date) %>%
      dplyr::mutate(class = "nasal_microbiome")
  ) %>%
  dplyr::filter(!is.na(collection_date)) %>%
  dplyr::distinct(subject_id_random, collection_date, class, .keep_all = TRUE)

library(plyr)
temp_20230607 <-
  temp_data %>%
  plyr::dlply(.variables = .(subject_id_random)) %>%
  purrr::map(function(x) {
    number = nrow(x)
    day = range(x$collection_date)[2] -
      range(x$collection_date)[1]
    data.frame(day = as.numeric(day), number)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

range(temp_20230607$day)
mean(temp_20230607$day)
median(temp_20230607$day)
mean(temp_20230607$number)

table(temp_data$class)


plot <-
  temp_data %>%
  ggplot(aes(collection_date, subject_id_random)) +
  geom_point(aes(color = class, shape = class), alpha = 0.5) +
  scale_color_manual(values = sample_class_color) +
  labs() +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Collection date", y = "")
plot

# ggsave(plot,
#        filename = "sample_collection.pdf",
#        width = 7,
#        height = 7)


transcriptome_sample_info$class <- "transcriptome"
proteomics_sample_info$class <- "proteomics"
metabolomics_sample_info$class <- "metabolomics"
cytokine_sample_info$class <- "cytokine"
clinical_test_sample_info$class <- "clinical_test"
lipidomics_sample_info$class <- "lipidomics"
gut_microbiome_sample_info$class <- "gut_microbiome"
skin_microbiome_sample_info$class <- "skin_microbiome"
oral_microbiome_sample_info$class <- "oral_microbiome"
nasal_microbiome_sample_info$class <- "nasal_microbiome"


temp_data <-
  list(
    transcriptome_sample_info,
    proteomics_sample_info,
    metabolomics_sample_info,
    cytokine_sample_info,
    clinical_test_sample_info,
    lipidomics_sample_info,
    gut_microbiome_sample_info,
    skin_microbiome_sample_info,
    oral_microbiome_sample_info,
    nasal_microbiome_sample_info
  ) %>%
  purrr::map(function(x) {
    x$collection_date <-
      as.character(x$collection_date)
    
    x$collection_date[is.na(x$collection_date)] <-
      x$collection_date[is.na(x$collection_date)] %>%
      purrr::map(function(i) {
        sample(letters, 20) %>%
          paste(collapse = "")
      }) %>%
      unlist()
    
    x <-
      x %>%
      dplyr::mutate(
        sample_id = paste(subject_id, collection_date, sep = "_"),
        class = x$class[1]
      ) %>%
      dplyr::select(sample_id, class) %>%
      dplyr::distinct(sample_id, .keep_all = TRUE)
    
    x
    
  })





temp_data <-
  temp_data %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  tidyr::pivot_wider(names_from = "class", values_from = "class") %>%
  dplyr::mutate(
    transcriptome = case_when(is.na(transcriptome) ~ FALSE,
                              TRUE ~ TRUE),
    proteomics = case_when(is.na(proteomics) ~ FALSE,
                           TRUE ~ TRUE),
    metabolomics = case_when(is.na(metabolomics) ~ FALSE,
                             TRUE ~ TRUE),
    cytokine = case_when(is.na(cytokine) ~ FALSE,
                         TRUE ~ TRUE),
    clinical_test = case_when(is.na(clinical_test) ~ FALSE,
                              TRUE ~ TRUE),
    lipidomics = case_when(is.na(lipidomics) ~ FALSE,
                           TRUE ~ TRUE),
    gut_microbiome = case_when(is.na(gut_microbiome) ~ FALSE,
                               TRUE ~ TRUE),
    skin_microbiome = case_when(is.na(skin_microbiome) ~ FALSE,
                                TRUE ~ TRUE),
    oral_microbiome = case_when(is.na(oral_microbiome) ~ FALSE,
                                TRUE ~ TRUE),
    nasal_microbiome = case_when(is.na(nasal_microbiome) ~ FALSE,
                                 TRUE ~ TRUE)
  )

library(ComplexUpset)

plot =
  upset(
    data = temp_data,
    name = "",
    intersect = colnames(temp_data)[-1],
    # set_sizes = TRUE,
    min_degree = 4,
    # group_by='sets',
    # min_size = 3,
    sort_sets = FALSE,
    sort_intersections = "ascending",
    sort_intersections_by = c('degree', "cardinality"),
    stripes = alpha(unname(omics_color[colnames(temp_data)[-1]]), 0.8),
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        mode = "intersect",
        text = list(angle = 90, size = 2),
        mapping = aes(fill = 'bars_color')
      ) +
        scale_fill_manual(values = c('bars_color' = 'black'),
                          guide = 'none') +
        theme_classic() +
        labs(x = "", y = "Sample number") +
        scale_y_continuous(expand = expansion(mult = c(0, 0))) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    ),
    matrix = (
      intersection_matrix(
        geom = geom_point(shape = 16,
                          size = 1),
        segment = geom_segment(linetype = 1),
        outline_color = list(active = 'white',
                             inactive = 'white')
      )
      # scale_color_manual(
      #   values=c('TRUE' = 'black', 'FALSE' = 'grey')
      #   # labels=c('TRUE'='yes', 'FALSE'='no'),
      #   # breaks=c('TRUE', 'FALSE'),
      #   # name=''
      # )
    )
    # queries=list(
    #   upset_query(
    #     intersect=c('RNA'),
    #     color=omics_color["RNA"],
    #     fill=omics_color["RNA"],
    #     only_components=c('intersections_matrix', 'Intersection size')
    #   )
    # )
  )

plot
# ggsave(plot,
#        filename = "upset_plot_of_all_samples.pdf",
#        width = 12,
#        height = 5)
