#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

fin_score_table <- args[1]
fin_candidate_scg_table <- args[2]
fin_genome_list <- args[3]
fin_scg_list <- args[4]
fout_scg_matrix <- args[5]

scores <- read_csv(fin_score_table)
candidate_scgs <- read_csv(fin_candidate_scg_table)
genomes <- read_lines(fin_genome_list)
scgs <- read_lines(fin_scg_list)

genes <- 
  scores %>%
  filter(genome %in% genomes) %>%
  filter(profile %in% scgs) %>%
  left_join(candidate_scgs) %>%
  filter(score >= gene_cutoff) %>%
  add_count(genome, profile) %>%
  filter(n == 1) %>%
  select(gene, genome, profile) 

scg_matrix <- spread(genes, value = gene, key = profile)

write_csv(scg_matrix, fout_scg_matrix)
  
