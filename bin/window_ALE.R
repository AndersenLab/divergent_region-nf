#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

ale <- data.table::fread(args[1], col.names = c("CHROM", "START_BIN", "END_BIN", "temp_chrom","temp_start","temp_end","DEPTH", "ALE")) %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(DEPTH = mean(DEPTH),
                   ALE = mean(ALE)) %>%
  dplyr::mutate(SAMPLE = args[2])


write.table(ale, file = glue::glue("{args[2]}_window_ale.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")