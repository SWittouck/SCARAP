# Progenomics: toolkit for prokaryotic comparative genomics.

Progenomics is a general toolkit-under-construction for comparative genomics of prokaryotes. It should be able to handle large genome datasets of small to medium sequence divergence (i.e., genomes from the same species, genus, family and possibly order). What is currently implemented is a pipeline to get the __core genome__ for up to thousands of genomes overnight on a decent desktop computer. A __pangenome pipeline__ is planned for the near future.

Progenomics depends on [OrthoFinder](https://github.com/davidemms/OrthoFinder) for gene family inference. 

## Dependencies

* [Python3](https://www.python.org/) version >= 3.6.7
* Python libraries:
    * [Biopython](https://biopython.org/) version >= 1.67
    * [pandas](https://pandas.pydata.org/) version >= 0.24.1
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) version >= 2.6.0
* [MCL](https://www.micans.org/mcl/index.html?sec_software) version >= 14-137
* [OrthoFinder](https://github.com/davidemms/OrthoFinder) version >= 2.1.2
* [mafft](https://mafft.cbrc.jp/alignment/software/) version >= 7.407
* [HMMER](http://hmmer.org/) version >= 3.1b2,
* [R](https://www.r-project.org/) version >= 3.5.1
* R packages:
    * ROCR version >= 1.0.7
    * tidyverse version >= 1.2.1

## Core genome pipeline

**How it works**

Progenomics works in four stages to be able to rapidly determine single-copy core genes (SCGs) for large genome datasets:

1. Gene family inference on a small, random subset of N (e.g. 30) seed genomes, using OrthoFinder. 
2. Selection of candidate SCGs by requiring single-copy presence in at least K (e.g. 25) out of N seed genomes. 
3. Search for the candidate SCGs in all genomes, using HMMER, and training of SCG-specific score cutoffs.  
4. Selection of the definitive SCGs by enforing that a candidate SCG is present in a single copy in P% (e.g. 95%) of all genomes.

**Tutorial**

In this workflow, we start from a set of genomes (up to ~ 3000 if you want to run overnight on a decent desktop computer) and we want to extract a complete set of single-copy core genes (SCGs) of those genomes. As input data, we need to have a set of predicted protein sequences (.faa file) for each genome. We need to supply the __paths to these .faa files__ to progenomics as a single txt file. If all .faa files are in the same folder, we can create such a file as follows:

    ls genomes/*.faa.gz > genomepaths.txt

Next, we extract candidate SCGs from the genomes and search for them in the complete genome dataset (__steps 1 - 3__). Candidate SCGs are gene families that are present in a single copy in at least K genomes out of N randomly chosen seed genomes. We could set K to 25 and N to 30, for example:

    progenomics prepare_candidate_scgs \
     --fin_genomepaths genomepaths.txt \
     --n_seed_genomes 30 \
     --min_presence_in_seeds 25 \
     --dout cand_scgs \
     --threads 8

We now select the definitive SCGs from the candidates by requiring that they are present in a single copy in P% of the total number of genomes (__step 4__). For P = 95, this gives us:

    progenomics select_scgs \
     --fin_score_table cand_scgs/score_table.csv \
     --fin_candidate_scg_table cand_scgs/candidate_scg_table.csv \
     --candidate_scg_cutoff 0.95 \
     --fout_scg_list scg_list.txt \
     --fout_genome_table genome_table.csv

As a result, we get two output files: a list with SCG names and a table with for each genome, the percentage "completeness" and "redundancy"; those two measures can be used for __genome quality filtering__. For this demonstration, we keep all of our genomes and save their names in a txt file:

    cut -f 1 genome_table.csv > selected_genomes.txt

Finally, we can construct a __SCG matrix__ where the rows are genomes, the columns are SCGs and the cells contain the actual names of individual genes:

    progenomics construct_scg_matrix \
     --fin_score_table cand_scgs/score_table.csv \
     --fin_candidate_scg_table cand_scgs/candidate_scg_table.csv \
     --fin_genome_list selected_genomes.txt \
     --fin_scg_list scg_list.txt \
     --fout_scg_matrix scg_matrix.csv

## Pangenome pipeline

Coming soon!

## License

Progenomics is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/Roary/blob/master/GPL-LICENSE).

## Feedback

All feedback and suggestions very welcome at stijn.wittouck[at]uantwerpen.be. You are of course also welcome to file [issues](https://github.com/SWittouck/progenomics/issues). 

## Citation

If you use progenomics in your publication, please try to cite the following preprint:

[Wittouck, Stijn, Sander Wuyts, Conor J Meehan, Vera van Noort, and Sarah Lebeer. 2019. “A Genome-Based Species Taxonomy of the Lactobacillus Genus Complex.” BioRxiv, January, 537084. doi:10.1101/537084.](https://www.biorxiv.org/content/10.1101/537084v1) 

If citing preprints is not allowed for your journal, don't worry about it. Hopefully, a peer-reviewed publication will be available soon! 
