#!/usr/bin/env Rscript

library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)

fin_genome_table <- args[1]
missingness_cutoff <- args[2]
redundancy_cutoff <- args[3]
fout_genome_list <- args[4]

fin_genome_table %>%
  read_csv() %>%
  filter(missingness < !! missingness_cutoff) %>%
  filter(redundancy < !! redundancy_cutoff) %>%
  select(genome) %>%
  write_csv(fout_genome_list, col_names = F)
