# progenomics

Toolkit for prokaryotic comparative genomics.

## Dependencies

* Biopython version 1.67
* OrthoFinder version 2.1.2
* blast version 2.6.0
* MCL version 14-137
* HMMER version 3.1b2,
* R version 3.5.1
* R packages:
    * ROCR version 1.0.7
    * tidyverse version 1.2.1

## Workflows

**Extraction of single-copy core genes**

In this workflow, we start from a set of genomes (up to ~ 3000 if you want to run overnight on a decent desktop) and we want to extract a complete set of single-copy core genes (SCGs) of those genomes. As input data, we need to have a set of predicted protein sequences (.faa file) for each genome. We need to supply the paths to those .faa files to progenomics as a single txt file. If all .faa files are in the same folder, we can create such a file as follows:

    ls genomes/*.faa.gz > genomepaths.txt

Next, we extract candidate SCGs from the genomes. These are gene families that are present in a single copy in at least K genomes out of N randomly chosen seed genomes. We could set K to 25 and N to 30, for example:

    progenomics prepare_candidate_scgs \
     --fin_genomepaths genomepaths.txt \
     --n_seed_genomes 30 \
     --min_presence_in_seeds 25 \
     --dout cand_scgs \
     --threads 8

We now select the "real" SCGs from the candidates by requiring that they are present in a single copy in P% of the total number of genomes. For P = 95, this gives us:

    progenomics select_scgs \
     --fin_score_table cand_scgs/score_table.csv \
     --fin_candidate_scg_table cand_scgs/candidate_scg_table.csv \
     --candidate_scg_cutoff 0.95 \
     --fout_scg_list scg_list.txt \
     --fout_genome_table genome_table.csv

As a result, we get two output files: a list with SCG names and a table with for each genome, the percentage "missingness" and "redundancy"; those two measures can be used for genome quality control. For this demonstration, we keep all of our genomes and save their names in a txt file:

    cut -f 1 genome_table.csv > selected_genomes.txt

Finally, we can construct a "SCG matrix" where the rows are genomes, the columns are SCGs and the cells contain the actual names of individual genes:

    progenomics construct_scg_matrix \
     --fin_score_table cand_scgs/score_table.csv \
     --fin_candidate_scg_table cand_scgs/candidate_scg_table.csv \
     --fin_genome_list selected_genomes.txt \
     --fin_scg_list scg_list.txt \
     --fout_scg_matrix scg_matrix.csv
