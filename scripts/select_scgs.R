#!/usr/bin/env Rscript

library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)

fin_score_table <- args[1]
fin_candidate_scg_table <- args[2]
candidate_scg_cutoff <- args[3]
fout_scg_list <- args[4]
fout_genome_table <- args[5]

scores <- read_csv(fin_score_table)
candidate_scgs <- read_csv(fin_candidate_scg_table)

n_genomes <-
  scores$genome %>%
  unique() %>%
  length()

scgs <-
  candidate_scgs %>%
  filter(single_copy_presence / !! n_genomes >= !! candidate_scg_cutoff) %>%
  select(profile)

write_csv(scgs, fout_scg_list, col_names = F)

n_scgs <- nrow(scgs)

genomes <-
  scores %>%
  right_join(scgs) %>%
  left_join(candidate_scgs) %>%
  filter(score >= gene_cutoff) %>%
  count(genome, profile) %>%
  group_by(genome) %>%
  summarize(n_present = n(), n_redundant = sum(n > 1)) %>%
  mutate(
    completeness = n_present / !! n_scgs,
    redundancy = n_redundant / !! n_scgs
  ) %>%
  select(genome, completeness, redundancy)

write_csv(genomes, fout_genome_table)
