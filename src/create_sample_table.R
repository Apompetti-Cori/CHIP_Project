library(data.table)
library(tidyverse)
setwd("/home/apompetti/data/14Dec2023_CHIP_Exploration/")
here::i_am(".here")
vaf_files <- list.files(here::here("data/vaf_anno_files"), 
                        pattern = "CUH.*vaf_anno.tsv", full.names = TRUE, recursive = TRUE)

x <- vroom::vroom(vaf_files, id = "filepath")

y <- x %>%
  mutate(
    sample_id = str_replace_all(filepath,
                                ".+/|_vaf_anno.tsv|\\._vaf_anno.tsv",
                                ""
    ),
    batch = str_extract(filepath, 
                        "/mnt/data/data_ap/14Dec2023_CHIP_Exploration/data/vaf_anno_files/(.+)/", 
                        group = 1)
  ) %>%
  mutate(record_id = str_replace_all(sample_id,"CUH0+", "")) %>%
  mutate(is_rep = str_detect(sample_id, "R")) %>%
  mutate(sample_id = str_replace_all(sample_id,"CUH0+|R", "")) %>%
  mutate(record_id = case_when(nchar(sample_id) == 4 ~ paste0("CUH0", record_id),
                               nchar(sample_id) == 3 ~ paste0("CUH00", record_id),
                               nchar(sample_id) == 2 ~ paste0("CUH000", record_id))) %>%
  mutate(sample_id = case_when(nchar(sample_id) == 4 ~ paste0("CUH0", sample_id),
                               nchar(sample_id) == 3 ~ paste0("CUH00", sample_id),
                               nchar(sample_id) == 2 ~ paste0("CUH000", sample_id))) %>%
  mutate(
    chip_relevant = case_when(
      vaf >= 0.02 & vaf <= 0.4 & dp > 99 ~ TRUE,
      TRUE ~ FALSE
    ),
    impact = factor(impact, levels = c("MODIFIER", "LOW", "MODERATE", "HIGH")),
    chip_relevant = factor(chip_relevant, levels = c(TRUE, FALSE))
  )

md <- read_tsv(here::here("data/metadata/2023-12-14_CUH_methylation_ages.tsv"))

y <- left_join(y, md, by = "sample_id") %>% 
  mutate(
    abserror = abs(error),
    error_classifier = case_when(
      error > 0 & error <= 1 ~ "minimal_age_increase",
      error > 1 & error <= 5 ~ "slight_age_increase",
      error > 5 & error <= 10 ~ "moderate_age_increase",
      error > 10 & error <= 20 ~ "high_age_increase",
      error > 20  ~ "severe_age_increase",
      error < 0 & error >= -1  ~ "minimal_age_decrease",
      error < -1 & error >= -5 ~ "slight_age_decrease",
      error < -5 & error >= -10 ~ "moderate_age_decrease",
      error < -10 & error >= -20 ~ "high_age_decrease",
      error < -20  ~ "severe_age_decrease",
      TRUE ~ "no_change"
    ),
    abserror_classifier = case_when(
      abserror > 0 & abserror <= 1 ~ "minimal_age_change",
      abserror > 1 & abserror <= 5 ~ "slight_age_change",
      abserror > 5 & abserror <= 10 ~ "moderate_age_change",
      abserror > 10 & abserror <= 20 ~ "high_age_change",
      abserror > 20  ~ "severe_age_change",
      TRUE ~ "no_change"
    )
  ) %>%
  mutate(error_classifier = factor(error_classifier, levels = c("minimal_age_decrease", "minimal_age_increase",
                                                                "slight_age_decrease", "slight_age_increase",
                                                                "moderate_age_decrease", "moderate_age_increase",
                                                                "high_age_decrease", "high_age_increase",
                                                                "severe_age_decrease", "severe_age_increase")),
         abserror_classifier = factor(abserror_classifier, levels = c("minimal_age_change", "slight_age_change",
                                                                      "moderate_age_change", "high_age_change",
                                                                      "severe_age_change"))) %>% 
  filter(!is.na(error))

gnomad <- read_tsv(here::here("data/metadata/gnomad_chip_loci.tsv")) %>%
  rename(
    chrom = Chromosome,
    pos = Position
  )

y <- left_join(y, gnomad, by = c("chrom", "pos"))

fwrite(y, file = here::here("results/tables/CHIP_sample_table.tsv"), 
       sep = "\t")
saveRDS(y, here::here("results/rds/CHIP_sample_table.rds"))