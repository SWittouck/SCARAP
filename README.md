# Progenomics: toolkit for prokaryotic comparative genomics

Progenomics is a toolkit for fast and scalable comparative genomics. It has been designed for prokaryotes but should work for eukaryotic genomes as well. Progenomics can handle large genome datasets on a range of taxonomic levels; it has been tested on datasets up until the order level. Its most useful features are fast pangenome inference, sensitive search of query sequences in a pangenome, rapid core genome inference using a heuristic strategy and the construction of concatenated core gene alignments ("supermatrices"). 

## Dependencies

Essential dependencies: 

* [Python3](https://www.python.org/) version >= 3.6.7
* Python libraries:
    * [biopython](https://biopython.org/) version >= 1.67
    * [ete3](http://etetoolkit.org/) version >= 3.1.1
    * [scipy](https://www.scipy.org/) version >= 1.4.1
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) version >= 7.407
* [MMseqs2](https://github.com/soedinglab/MMseqs2) version >= d36dea2

Dependencies for the search module and core pipeline:

* [HMMER](http://hmmer.org/) version >= 3.1b2

Dependencies when using [OrthoFinder](https://github.com/davidemms/OrthoFinder) for pangenome inference:

* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) version >= 2.6.0
* [MCL](https://www.micans.org/mcl/index.html?sec_software) version >= 14-137
* [OrthoFinder](https://github.com/davidemms/OrthoFinder) version >= 2.1.2

## Usage

Progenomics is able to perform a number of specific tasks related to prokaryotic core and pangenomes (see also `progenomics -h`):

* `pan`: infer a pangenome from a set of faa files
* `build`: build a profile HMM database for a core/pangenome
* `search`: search query genes in a core/pangenome database
* `checkgenomes`: assess the quality of genomes in a core genome
* `checkgroups`: assess the quality of orthogroups in a core genome
* `filter`: filter the genomes/orthogroups in a pangenome
* `supermatrix`: construct a concatenated core orthogroup alignment from a core genome

A full core and pangenome pipeline are also implemented:

* `pan-pipeline`: infer a pangenome, build a profile HMM database and train score cutoffs from a set of faa files
* `core-pipeline`: infer a core genome, build a profile HMM database and train score cutoffs from a set of faa files

### Inferring a pangenome

If you want to infer the pangenome of a set of genomes, you only need their faa files (fasta files with protein sequences) as input. If the faa files are stored in a folder `faas`, you can infer the pangenome using 16 threads by running: 

    ls faas/*.faa > faapaths.txt
    progenomics pan faapaths.txt pan -t 16
    
The pangenome will be stored in `pan/pangenome.tsv`. 

The above example will use the builtin "FH" strategy to infer the pangenome; it is fast and scales more or less linearly with the number of input genomes. If you prefer to use OrthoFinder for pangenome inference, you can run:

    ls faas/*.faa > faapaths.txt
    progenomics pan faapaths.txt pan -d O-B -t 16
    
This will be a bit slower though.

### Searching a pangenome

Disclaimer: many aspects of this pipeline can still change, especially the way that hmmer score cutoffs are trained.

**Quick version**

If you want to infer a pangenome of your genomes as well as build a pangenome database that you can later query with one or more genes of interest, you can run:

    ls faas/*.faa > faapaths.txt
    progenomics pan-pipeline faapaths.txt pan -t 16

If you then want to identify whether some genes of interest (let's say in a file called `querygenes.fasta`) are present in the pangenome database, you can run:

    echo querygenes.fasta > querypath.txt
    progenomics search querypath.txt pangenome/db hits

This will produce a `hits` output folder with the file `hits.tsv`.

**Detailed version**

The pangenome pipeline can also be performed by running individual tasks:

    ls faas/*.faa > faapaths.txt
    progenomics pan faapaths.txt pan
    progenomics build faapaths.txt pan/pangenome.tsv db
    progenomics search faapaths.txt db pan2 -s pan -p pan/pangenome.tsv

The final `search` step is required because it will train a hmmer score cutoff for each profile HMM in the pangenome database and add these cutoffs to the database. In addition, it produces an orthogroup assignment for each protein in the set of input genomes (`pan2/hits.tsv`). Importantly, these assignments are not always the same as the orthogroup assignments listed in `pan/pangenome.tsv` because they are produced by a hmmer search with orthogroup-specific cutoffs, while the original orthogroup assignments have been produced by the pangenome inference process. A comparison between these two strategies of orthogroup assignment could be interesting.

### Inferring a core genome only

Let's say we want to infer the core genome for a set of genomes and we have one faa file (amino acid sequences of predicted genes) per genome in the folder `faas`. This can be done using the core pipeline of progenomics, which can be a lot faster than full pangenome inference. 

**Quick version**

To get the core genome, we can simply run the following commands:

    ls faas/*.faa > faapaths.txt
    progenomics core-pipeline faapaths.txt core -t 16

This will create the output folder `core`, run progenomics with 16 threads and produce the file `coregenome.tsv`. This output tsv file contains the core orthogroups and has the columns gene, genome and orthogroup.

If we now want to construct a supermatrix (concatenated alignment) of these core orthogroups, we could do it as follows:

    progenomics supermatrix faapaths.txt core/coregenome.tsv supermatrix

This will create a `supermatrix` output folder, with in it a file supermatrix.fasta.

And that's it! Three lines of code to get from the faa files to the supermatrix fasta file, ready to start constructing your phylogenetic tree.

**Detailed version**

If we want more fine-grained control, we could achieve the same result by running individual progenomics tasks. These individual tasks also give insight in how the core genome pipeline actually works.

**Step 1:** infer the pangenome of a random subset of seed genomes (e.g. 30).

    mkdir seeds cands
    ls faas/*.faa > faapaths.txt
    shuf -n 30 faapaths.txt > seeds/faapaths.txt
    progenomics pan seeds/faapaths.txt seeds/pan

**Step 2:** build a profile HMM database of "candidate core orthogroups" that are present in at least M seed genomes (e.g. 25).

    progenomics build seeds/faapaths.txt seeds/pan/pangenome.tsv cands/db -m 25

**Step 3:** identify the candidate core genes in the full set of genomes by searching all proteins of all genomes against the database of candidate core genes.

    progenomics search faapaths.txt cands/db cands/core -y core

**Step 4:** identify the core genes from the candidates by imposing a minimum percentage presence cutoff (e.g. 98%) in the full set of genomes.

    progenomics checkgroups cands/core/coregenome.tsv cands/groups
    awk '{ if ($2 > 0.98) { print $1 } }' cands/groups/orthogroups.tsv \
      > orthogroups.txt
    progenomics filter cands/core/coregenome.tsv core -o orthogroups.txt

The output folder `core` will now contain the file coregenome.tsv.

## License

Progenomics is free software, licensed under [GPLv3](https://github.com/SWittouck/progenomics/blob/master/LICENSE).

## Feedback

All feedback and suggestions very welcome at stijn.wittouck[at]uantwerpen.be. You are of course also welcome to file [issues](https://github.com/SWittouck/progenomics/issues).

## Citation

When you use progenomics for your publication, please cite:

[Wittouck, S., Wuyts, S., Meehan, C. J., van Noort, V., & Lebeer, S. (2019). A Genome-Based Species Taxonomy of the Lactobacillus Genus Complex. mSystems, 4(5), e00264–19.](https://doi.org/10.1128/mSystems.00264-19)
