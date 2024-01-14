library(tidyverse)
library(tidymass)

setwd(masstools::get_project_wd())
setwd("3-data_analysis/other/")

rm(list = ls())

load("object_integration_ms2_POS")

object

load("hmdb_ms2.rda")
load("massbank_ms2.rda")
load("metlin_ms2.rda")
load("mona_ms2.rda")
load("mpsnyder_hilic_ms2.rda")
load("mpsnyder_rplc_ms2.rda")
load("nist_ms2.rda")

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "positive",
    threads = 8,
    database = hmdb_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "positive",
    threads = 8,
    database = massbank_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "positive",
    threads = 8,
    database = metlin_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "positive",
    threads = 8,
    database = mona_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "positive",
    threads = 8,
    database = mpsnyder_hilic_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "positive",
    threads = 8,
    database = mpsnyder_rplc_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "positive",
    threads = 8,
    database = nist_ms2
  )


object_pos <-
  object

save(object_pos, file = "object_pos")









load("object_integration_ms2_NEG")

object

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "negative",
    threads = 8,
    database = hmdb_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "negative",
    threads = 8,
    database = massbank_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "negative",
    threads = 8,
    database = metlin_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "negative",
    threads = 8,
    database = mona_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "negative",
    threads = 8,
    database = mpsnyder_hilic_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "negative",
    threads = 8,
    database = mpsnyder_rplc_ms2
  )

object <-
  annotate_metabolites_mass_dataset(
    object = object,
    ms1.match.ppm = 5,
    rt.match.tol = 1000000,
    polarity = "negative",
    threads = 8,
    database = nist_ms2
  )

object_neg <-
  object

save(object_neg, file = "object_neg")
