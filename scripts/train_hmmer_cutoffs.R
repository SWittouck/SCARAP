#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(ROCR)

args = commandArgs(trailingOnly=TRUE)

fin_score_table <- args[1]
fout_scg_table <- args[2]

genes <- read_csv(fin_score_table)

cutoff <- function(score, is_best_copy) {

  if(all(is_best_copy)) return(min(score))

  pred <- prediction(score, is_best_copy)
  perf <- performance(pred, "f")

  tibble(
    cutoff = perf@x.values[[1]],
    f = perf@y.values[[1]]
  ) %>%
    filter(f == max(f, na.rm = T)) %>%
    slice(1) %>%
    pull(cutoff)

}

# add variable is_best_copy (of the genome)
genes <- genes %>%
  group_by(profile) %>%
  arrange(- score) %>%
  mutate(is_best_copy = ! duplicated(genome)) %>%
  ungroup()

# create profiles table with score cutoff per profile
profiles <- genes %>%
  group_by(profile) %>%
  summarize(gene_cutoff = cutoff(score, is_best_copy))

# determine total number of genomes
n_genomes <- genes$genome %>%
  unique() %>%
  length()

# count genomes where profile is present
profiles_presence <- genes %>%
  left_join(profiles) %>%
  filter(score >= gene_cutoff) %>%
  count(profile, genome) %>%
  group_by(profile) %>%
  summarize(
    single_copy_presence = sum(n == 1),
    multi_copy_presence = sum(n > 1)
  )

# add genome presence to profiles table
profiles <- left_join(profiles, profiles_presence)

# write profiles table to fout_scg_table
write_csv(profiles, path = fout_scg_table, col_names = T)
