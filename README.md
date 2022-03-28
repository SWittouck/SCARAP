# SCARAP: toolkit for comparative genomics of prokaryotes

SCARAP is a toolkit with modules for various tasks related to comparative genomics of prokaryotes. SCARAP has been designed to be fast and scalable. Its main feature is pangenome inference, but it also has modules for direct core genome inference (without inferring the full pangenome), subsampling representatives from a (large) set of genomes and constructing a concatenated core gene alignment ("supermatrix") that can later be used for phylogeny inference.  SCARAP has been designed for prokaryotes but should work for eukaryotic genomes as well. It can handle large genome datasets on a range of taxonomic levels; it has been tested on datasets with prokaryotic genomes from the species to the order level. 

## Dependencies

Essential dependencies: 

* [Python3](https://www.python.org/) version >= 3.6.7
* Python libraries:
    * [biopython](https://biopython.org/) version >= 1.67
    * [ete3](http://etetoolkit.org/) version >= 3.1.1
    * [scipy](https://www.scipy.org/) version >= 1.4.1
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) version >= 7.407
* [MMseqs2](https://github.com/soedinglab/MMseqs2) release 11, 12 or 13

## Usage

SCARAP is able to perform a number of specific tasks related to prokaryotic comparative genomics (see also `scarap -h`). 

The most useful modules of SCARAP are the following: 

* `pan`: infer a pangenome from a set of faa files
* `core`: infer a core genome from a set of faa files
* `sample`: sample a subset of representative genomes

Modules for other useful tasks are also available: 

* `build`: build a profile database for a core/pangenome
* `search`: search query genes in a profile database
* `checkgenomes`: assess the genomes in a core genome
* `checkgroups`: assess the orthogroups in a core genome
* `filter`: filter the genomes/orthogroups in a pangenome
* `concat`: construct a concatenated core orthogroup alignment from a core genome
* `fetch`: fetch sequences and store in fasta per orthogroup

### Inferring a pangenome

If you want to infer the pangenome of a set of genomes, you only need their faa files (fasta files with protein sequences) as input. If the faa files are stored in a folder `faas`, you can infer the pangenome using 16 threads by running: 

      scarap pan ./faas ./pan -t 16
    
The pangenome will be stored in `pan/pangenome.tsv`. 

The pangenome is stored in a "long format": a table with the columns gene, genome and orthogroup. 

### Inferring a core genome 

If you want to infer the core genome of a set of genomes directly, without inferring the full pangenome first, you can also do this with SCARAP. The reason you might want to do this, is because it is faster and because you sometimes don't need more than the core genome (e.g. when you are planning to infer a phylogeny). 

You can infer the core genome, given a set of faa files in a folder `faas`, in the following way:

      scarap core ./faas ./core -t 16
      
The core genome will be stored in `core/genes.tsv`. 

### Subsampling a set of genomes 

If you have a (large) dataset of genomes that you wish to subsample in a representative way, you can do this using the `sample` module. You will need to precompute the pangenome or core genome to do this; SCARAP calculates average nucleotide identity (ANI) or core nucleotide identity (cANI) values in the subsampling process, and it uses the single-copy orthogroups from a pan- or core genome to do this. 

For example, if you want to sample 100 genomes given a set of faa files in a folder `faas`: 

      scarap core ./faas ./core -t 16
      scarap sample ./faas ./core/genes.tsv ./representatives -m 100 -t 16
      
The representative genomes will be stored in `representatives/seeds.txt`. 

### Building a "supermatrix" for a set of genomes

You can build a concatenated alignment of core genes ("supermatrix") for a set of genomes using the `concat` module. 

Let's say you want to build a supermatrix of 100 core genes for a set of genomes, with faa files given in a folder `faas`: 

      scarap core ./faas ./core -m 100 -t 16
      scarap concat ./faas ./core/genes.tsv ./supermatrix -t 16
      
The amino acid supermatrix will be saved in `supermatrix/supermatrix_aas.fasta`. 
      
If you want to produce a nucleotide-level supermatrix, this can be achieved by giving a folder with ffn files (nucleotide sequences of predicted genes) as an additional argument: 

      scarap concat ./faas ./core/genes.tsv ./supermatrix -n ./ffns -t 16
      
The nucleotide-level supermatrix will be saved in `supermatrix/supermatrix_nucs.fasta`. 

## License

SCARAP is free software, licensed under [GPLv3](https://github.com/SWittouck/scarap/blob/master/LICENSE).

## Feedback

All feedback and suggestions very welcome at stijn.wittouck[at]uantwerpen.be. You are of course also welcome to file [issues](https://github.com/SWittouck/scarap/issues).

## Citation

When you use SCARAP for your publication, please cite:

[Wittouck, S., Wuyts, S., Meehan, C. J., van Noort, V., & Lebeer, S. (2019). A Genome-Based Species Taxonomy of the Lactobacillus Genus Complex. mSystems, 4(5), e00264–19.](https://doi.org/10.1128/mSystems.00264-19)
