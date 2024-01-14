no_source()

rm(list = ls())
setwd(r4projects::get_project_wd())
source("1-code/100-tools.R")

library(tidyverse)
library(tidymass)

###load("data)
load(
  "3-data_analysis/plasma_transcriptome/data_preparation/object_cross_section_loess"
)

dir.create(
  "3-data_analysis/plasma_transcriptome/fuzzy_c_means_clustering/cross_section_loess/",
  recursive = TRUE
)

setwd(
  "3-data_analysis/plasma_transcriptome/fuzzy_c_means_clustering/cross_section_loess"
)

object_cross_section_loess

####clustering
library(Mfuzz)

object_cross_section_loess[1, , drop = TRUE] %>%
  unlist() %>%
  density() %>%
  plot()

######cluster all variables using the age range
object_cross_section_loess@sample_info$adjusted_age %>% range

object_cross_section_loess@sample_info$adjusted_age %>% plot

object_cross_section_loess@sample_info %>%
  ggplot(aes(adjusted_age)) +
  geom_histogram(binwidth = 5, color = "black")

range(object_cross_section_loess@sample_info$adjusted_age)

age_index <-
  data.frame(from = c(25, 40, 45, 50, 55, 60, 65),
             to =   c(40, 45, 50, 55, 60, 65, 75.5))

apply(age_index, 1, function(x) {
  sum(
    object_cross_section_loess@sample_info$adjusted_age > x[1] &
      object_cross_section_loess@sample_info$adjusted_age <= x[2]
  )
})

# temp_data <-
#   log(object_cross_section_loess@expression_data + 1, 2)

temp_data <-
  object_cross_section_loess@expression_data

temp_data_mean <-
  seq_len(nrow(age_index)) %>%
  purrr::map(function(i) {
    idx <- as.numeric(age_index[i, ])
    temp_data[, which(
      object_cross_section_loess@sample_info$adjusted_age > idx[1] &
        object_cross_section_loess@sample_info$adjusted_age <= idx[2]
    )] %>%
      apply(1, mean)
  }) %>%
  dplyr::bind_cols() %>%
  as.data.frame()

colnames(temp_data_mean) <-
  paste(age_index$from, age_index$to, sep = "_")

rownames(temp_data_mean) <-
  object_cross_section_loess@variable_info$variable_id

# temp_data_mean <-
#   temp_data_mean %>%
#   apply(1, function(x) {
#     x - mean(x)
#   }) %>%
#   t() %>%
#   as.data.frame()

idx = 3

object_cross_section_loess[idx, , drop = TRUE] %>%
  unlist() %>%
  density() %>%
  plot()

temp_data[idx,] %>%
  as.numeric() %>%
  density() %>%
  plot()

time <- colnames(temp_data_mean)

temp_data <- rbind(time, temp_data_mean)

row.names(temp_data)[1] <- "time"
rownames(temp_data)

write.table(
  temp_data,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")
data.s <- standardise(data)
# data.s <- data
m1 <- mestimate(data.s)
m1

plot <-
  Dmin(
    data.s,
    m = m1,
    crange = seq(2, 40, 2),
    repeats = 3,
    visu = TRUE
  )

plot <-
  plot %>%
  data.frame(distance = plot,
             k = seq(2, 40, 2)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21, size = 4, fill = "black") +
  # geom_smooth() +
  geom_segment(aes(
    x = k,
    y = 0,
    xend = k,
    yend = distance
  )) +
  theme_bw() +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "Cluster number",
       y = "Min. centroid distance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

plot

ggsave(plot,
       filename = "distance_k_number.pdf",
       width = 7,
       height = 7)

cluster_number <- 11

c <- mfuzz(data.s, c = cluster_number, m = m1)

# ####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
# center <- c$centers

membership_cutoff <- 0.5

center <-
  get_mfuzz_center(data = data.s,
                   c = c,
                   membership_cutoff = 0.5)
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

mfuzz.plot(
  eset = data.s,
  # min.mem = 0.6,
  cl = c,
  mfrow = c(3, 4),
  time.labels = time,
  new.window = FALSE
)

library(ComplexHeatmap)

###
cluster_color <-
  ggsci::pal_jama()(n = 7)[1:cluster_number]

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

####plot for each cluster
idx <- 1

temp_data <-
  data.s %>% as.data.frame() %>%
  t() %>% as.data.frame()

for (idx in 1:cluster_number) {
  cat(idx, " ")
  
  cluster_data <-
    cluster_info %>%
    # dplyr::filter(cluster == idx) %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > membership_cutoff)
  
  path <- paste("cluster", idx, sep = "_")
  dir.create(path)
  
  openxlsx::write.xlsx(
    cluster_data,
    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
    asTable = TRUE,
    overwrite = TRUE
  )
  
  temp_center <-
    center[idx, , drop = TRUE] %>%
    unlist() %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE) %>%
    dplyr::mutate(time = factor(time, levels = time)) %>%
    dplyr::mutate(time_point = time)
  
  temp_center$time <-
    temp_center$time %>%
    as.character() %>%
    stringr::str_split(pattern = "_") %>%
    lapply(function(x) {
      mean(as.numeric(x))
    }) %>%
    unlist()
  
  temp <-
    temp_data[cluster_data$variable_id, ] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "time",
      values_to = "value"
    ) %>%
    dplyr::mutate(time = factor(time, levels = unique(time))) %>%
    dplyr::mutate(time_point = time)
  
  temp$time <-
    temp$time %>%
    as.character() %>%
    stringr::str_split(pattern = "_") %>%
    lapply(function(x) {
      mean(as.numeric(x))
    }) %>%
    unlist()
  
  plot <-
    temp %>%
    dplyr::arrange(membership, variable_id) %>%
    dplyr::mutate(variable_id = factor(variable_id, levels = unique(variable_id))) %>%
    ggplot(aes(time, value, group = variable_id)) +
    geom_line(aes(color = membership), alpha = 0.7) +
    theme_bw() +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 12
      ),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score",
      title = paste(
        "Cluster ",
        idx,
        " (",
        nrow(cluster_data),
        " metabolic features)",
        sep = ""
      )
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    ) +
    geom_hline(yintercept = 0) +
    viridis::scale_color_viridis() +
    scale_x_continuous(breaks = c(temp_center$time),
                       labels = temp_center$time_point)
  
  plot
  
  ggsave(
    plot,
    filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
  
}

dim(cluster_data)

table(cluster_info$cluster)

cluster_info <-
  unique(cluster_info$cluster) %>%
  purrr::map(function(x) {
    temp <-
      cluster_info %>%
      # dplyr::filter(cluster == x) %>%
      dplyr::select(variable_id, paste0("X", x), cluster)
    colnames(temp)[2] <- "membership"
    temp <-
      temp %>%
      dplyr::filter(membership >= membership_cutoff)
    temp <-
      temp %>%
      dplyr::mutate(cluster_raw = cluster) %>%
      dplyr::mutate(cluster = x)
    temp
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

cluster_info %>%
  dplyr::count(cluster)

cluster_info %>%
  dplyr::filter(membership > 0.5) %>%
  dplyr::count(cluster)

final_cluster_info <-
  cluster_info

save(final_cluster_info, file = "final_cluster_info")

openxlsx::write.xlsx(
  final_cluster_info,
  file = "final_cluster_info.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)
