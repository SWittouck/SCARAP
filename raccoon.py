#!/usr/bin/env python3

import argparse
import sys
import os
import glob
import subprocess
import concurrent.futures
import shutil
import re
import random
from Bio.Align.Applications import MafftCommandline
import pandas as pd
import gzip
import pathlib

"""
Class to store relations and metadata of a gene
"""
class Gene:

    def __init__(self, name, genome, orthogroup):
        self.name = name
        self.genome = genome
        self.orthogroup = orthogroup

    def __str__(self):
        return(self.name)

    @property
    def genome(self):
        return(self._genome)

    @genome.setter
    def genome(self, genome):
        self._genome = genome
        genome.genes.append(self)

    @property
    def orthogroup(self):
        return(self._orthogroup)

    @orthogroup.setter
    def orthogroup(self, orthogroup):
        self._orthogroup = orthogroup
        orthogroup.genes.append(self)


"""
Class to store relations and metadata of a genome
"""
class Genome:

    def __init__(self, name):
        self.name = name
        self.genes = []

    def __str__(self):
        return(self.name)

    def add_gene(self, gene):
        gene.genome = self


"""
Class to store relations and metadata of an orthogroup
"""
class Orthogroup:

    def __init__(self, name):
        self.name = name
        self.genes = []

    def __str__(self):
        return(self.name)

    def add_gene(self, gene):
        gene.orthogroup = self

    def get_n_genomes(self):
        genomes = set()
        for gene in self.genes:
            genomes.add(gene.genome)
        return(len(genomes))


""" 
Class to store the list of genomes and list of orthogroups
"""
class Taxon:

    def __init__(self):
        self.genomes = []
        self.orthogroups = []
        
# low-level functions

def read_paths(fin_paths):
    with open(fin_paths) as hin_paths:
        paths = [p.strip() for p in hin_paths.readlines()]
    return(paths)
    
"""
Function to parse an "orthogroups.csv" file
Returns an object of class Taxon  
"""
def read_orthogroups(din_orthofinder):
    taxon = Taxon()
    hin_orthofinder_csv = open(din_orthofinder + "/Orthogroups.csv", "r")
    genome_names = next(hin_orthofinder_csv)
    genome_names = genome_names.strip().split("\t")
    for genome_name in genome_names:
        genome = Genome(genome_name)
        taxon.genomes.append(genome)
    for line in hin_orthofinder_csv:
        elements = line.split("\t")
        orthogroup_name = elements[0]
        del(elements[0])
        orthogroup = Orthogroup(orthogroup_name)
        taxon.orthogroups.append(orthogroup)
        for ix, genome in enumerate(taxon.genomes):
            element = elements[ix].strip()
            if element == "":
                continue
            for gene_name in element.split(","):
                gene = Gene(gene_name.strip(), genome, orthogroup)
    print("found " + str(len(taxon.orthogroups)) + " orthogroups")
    return(taxon)
    
def select_orthogroups(taxon, min_genomes):
    orthogroups = taxon.orthogroups
    orthogroups_selected = []
    for orthogroup in orthogroups:
        if orthogroup.get_n_genomes() >= min_genomes:
            orthogroups_selected.append(orthogroup)
    print("selected " + str(len(orthogroups_selected)) + " orthogroups")
    return(orthogroups_selected)
    
def read_domtbl(fin_domtbl):
    domtbl = pd.read_csv(
        fin_domtbl, 
        delim_whitespace = True, 
        comment = '#', header = None,
        usecols = [0, 3, 13, 15, 16],
        names = ['gene', 'profile', 'score', 'hmm_from', 'hmm_to'],
        dtype = {'gene': str, 'profile': str, 'score': float, 'hmm_from': int, 'hmm_to': int}
    )
    return(domtbl)

"""
Function to perform initial processing of a hmmer domtbl file. 
- Takes the best profile hit per gene. 
- Looks up the genome of each gene and adds it to the score table. 
Returns a score table in the form of a data frame. 
Can handle zipped genomes. 
"""
def process_domtbl(domtbl, fins_genomes):
  
    scores_to_keep = (domtbl
        .loc[:, ['gene', 'score']]
        .groupby('gene')
        .idxmax()
        .score
    )
    
    scores = (domtbl
        .iloc[scores_to_keep, :]
        .set_index("gene")
    )
    
    genes = pd.DataFrame()
    
    for fin_genome in fins_genomes:
        genome = (re.compile("([^/]+).faa.gz")
            .search(fin_genome)
            .group(1)
        )
        with gzip.open(fin_genome, mode = "rt") as hin:
            fasta = hin.read()
        genes_genome = (re.compile(">([^ ]+)")
            .findall(fasta)
        )
        genes_genome = pd.DataFrame({"gene": genes_genome})
        try: 
            genes_genome.loc[:, "genome"] = genome
            genes = genes.append(genes_genome)
        except ValueError:
            pass
    
    genes = genes.set_index("gene")
    
    scores = scores.merge(genes, how = 'left', left_index = True, right_index = True)

    return(scores)

def construct_scg_table(fin_score_table, fout_scg_table):
    with open(os.devnull, 'w') as devnull:
        subprocess.run(["Rscript", "scripts/construct_scg_table.R", fin_score_table, 
            fout_scg_table], stdout = devnull, stderr = devnull)

# low-level functions calling external tools

def run_orthofinder(fins_genomes, dout_orthofinder):
    os.makedirs(dout_orthofinder, exist_ok = True)
    for fin_genome in fins_genomes:
        genome = pathlib.PurePosixPath(fin_genome).name
        fout_genome = dout_orthofinder + "/" + genome
        shutil.copyfile(fin_genome, fout_genome)
        subprocess.run(["gunzip", fout_genome])
    subprocess.run(["orthofinder", "-M", "msa", "-os", "-t", "8",
        "-f", dout_orthofinder])

def run_mafft(fin_aa_seqs, fout_aa_seqs_aligned):
    print("started aligning " + pathlib.PurePosixPath(fin_aa_seqs).stem)
    mafft_cline = MafftCommandline(input = fin_aa_seqs)
    stdout, stderr = mafft_cline()
    with open(fout_aa_seqs_aligned, "w") as hout_aa_seqs_aligned:
        hout_aa_seqs_aligned.write(stdout)
    print("finished aligning " + pathlib.PurePosixPath(fin_aa_seqs).stem)
  
def run_mafft_parallel(fins_aa_seqs, fouts_aa_seqs_aligned):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_mafft,
            fins_aa_seqs, fouts_aa_seqs_aligned
        )

def run_hmmbuild(fin_aa_seqs_aligned, fout_hmm):
    subprocess.run(
        ["hmmbuild", fout_hmm, fin_aa_seqs_aligned],
        stdout = subprocess.PIPE,
    )
    print("constructed profile of " + pathlib.PurePosixPath(fin_aa_seqs_aligned).stem)

def run_hmmbuild_parallel(fins_aa_seqs_aligned, fouts_hmms):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_hmmbuild, 
            fins_aa_seqs_aligned, fouts_hmms
        )

def run_hmmpress(fins_hmms, dout_hmm_db):
    print("pressing hmm database for " + str(len(fins_hmms)) + " gene families")
    fout_hmm_db = dout_hmm_db + "/hmm_db"
    subprocess.run(
        " ".join(["cat"] + fins_hmms + [">", fout_hmm_db]),
        shell = True
    )
    subprocess.run(
        ["hmmpress", fout_hmm_db],
        stdout = subprocess.PIPE
    )

def run_hmmsearch(din_hmm_db, fins_genomes, fout_domtbl):
    fin_hmm_db = din_hmm_db + "/hmm_db"
    fout_temp_genomes = din_hmm_db + "/genomes.temp"
    subprocess.run(
        " ".join(["cat"] + fins_genomes + [">", fout_temp_genomes]),
        shell = True
    )
    subprocess.run(["hmmsearch", "-o", "/dev/null", "--domtblout", fout_domtbl, 
        fin_hmm_db, fout_temp_genomes])
    subprocess.run(["rm", fout_temp_genomes])
  
def run_gunzip(fin):
    subprocess.run(["gunzip", fin])
    
def run_gunzip_parallel(fins):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_gunzip, fins)
        
def run_gzip(fin):
    subprocess.run(["gzip", fin])
    
def run_gzip_parallel(fins):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_gzip, fins)
        
# high-level stuff

def construct_profile_db(fins_aa_seqs, dout):
  
    # define variables for output files and directories
    dout_alignments = dout + "/alignments"
    dout_profiles = dout + "/profiles"
    dout_hmm_db = dout + "/hmm_db"
    
    # create directories
    douts_to_create = [dout_alignments, dout_profiles, dout_hmm_db]
    for dout_to_create in douts_to_create:
        os.makedirs(dout_to_create, exist_ok = True)
        
    # define file lists
    families = [pathlib.PurePosixPath(fin).stem for fin in fins_aa_seqs]
    fouts_alignments = [dout_alignments + "/" + fam + ".aln"
        for fam in families]
    fouts_profiles = [dout_profiles + "/" + fam + ".hmm"
        for fam in families]
    
    # align families, build profiles, build hmm database
    run_mafft_parallel(fins_aa_seqs, fouts_alignments)
    run_hmmbuild_parallel(fouts_alignments, fouts_profiles)
    run_hmmpress(fouts_profiles, dout_hmm_db)
    
def construct_genome_table(fin_score_table, fin_candidate_scg_table,
    candidate_scg_cutoff, fout_genome_table):
    with open(os.devnull, 'w') as devnull:
        subprocess.run(["Rscript", "scripts/construct_genome_table", 
            fin_score_table, fin_candidate_scg_table, candidate_scg_cutoff, 
            fout_genome_table], stdout = devnull, stderr = devnull)
            
def construct_scg_matrix(fin_score_table, fin_genome_table,
    fin_genome_list, candidate_scg_cutoff, fout_scg_matrix):
    with open(os.devnull, 'w') as devnull:
        subprocess.run(["Rscript", "scripts/construct_scg_matrix", 
            fin_score_table, fin_genome_table, fin_genome_list, 
            candidate_scg_cutoff, fout_scg_matrix], stdout = devnull, 
            stderr = devnull)

def parse_arguments():
  
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "task",
        help = "the task you want to perform",
        choices = [
            "construct_profile_db",
            "prepare_candidate_scgs",
            "construct_genome_table",
            "construct_scg_matrix",
            "construct_supermatrix",
            "calculate_scnis"
        ]
    )
    task = parser.parse_args(sys.argv[1:2]).task
    sys.argv[0] = sys.argv[0] + " " + task
    parser = argparse.ArgumentParser()
    
    if task == "construct_profile_db":
        parser.add_argument(
            "--fin_aa_seqs_paths",
            help = "a txt file with the paths of the gene families to put in the database",
            required = True
        )
        parser.add_argument(
            "--dout",
            help = "the output directory for the alignments, profiles hmms and hmmer database",
            required = True
        )
        
    elif task == "prepare_candidate_scgs":
        parser.add_argument(
            "--fin_genomepaths",
            help = "a txt file with paths to the genomes as fastas of amino acid sequences",
            required = True
        )
        parser.add_argument(
            "--n_seed_genomes",
            help = "the number of seed genomes to use",
            required = True
        )
        parser.add_argument(
            "--min_presence_in_seeds",
            help = "the minimum number of seed genomes a candidate SCG needs to be present in",
            required = True
        )
        parser.add_argument(
            "--dout",
            help = "the output directory for results",
            required = True
        )
        
    elif task == "construct_genome_table":
        parser.add_argument(
            "--fin_score_table",
            help = "a table with information on hmmer scores",
            required = True
        )
        parser.add_argument(
            "--fin_candidate_scg_table",
            help = "a table with information on candidate scgs",
            required = True
        )
        parser.add_argument(
            "--candidate_scg_cutoff",
            help = "the minimum percentage single copy presence for a candidate scg",
            required = True
        )
        parser.add_argument(
            "--fout_genome_table",
            help = "path to store the genome table",
            required = True
        )
    
    elif task == "construct_scg_matrix":
        parser.add_argument(
            "--fin_score_table",
            help = "a table with information on hmmer scores",
            required = True
        )
        parser.add_argument(
            "--fin_candidate_scg_table",
            help = "a table with information on candidate scgs",
            required = True
        )
        parser.add_argument(
            "--fin_genome_list",
            help = "a list of genomes to include",
            required = True
        )
        parser.add_argument(
            "--candidate_scg_cutoff",
            help = "the minimum percentage single copy presence for a candidate scg, in the genome list",
            required = True
        )
        parser.add_argument(
            "--fout_scg_matrix",
            help = "path to store the sgc matrix",
            required = True
        )
        
    else:
        print("Argparse should have dropped an 'unrecognized option' error")
    
    args = parser.parse_args(sys.argv[2:])
    args.task = task
    
    return(args)
    
if __name__ == "__main__":

    print("")
    print("hi, this is raccoon")
    print("")
    
    args = parse_arguments()

    if args.task == "construct_profile_db":
    
        fins_aa_seqs = read_paths(args.fin_aa_seqs_paths)
        construct_profile_db(fins_aa_seqs, args.dout)

    elif args.task == "prepare_candidate_scgs":
      
        # define variables for output directories
        dout_orthofinder = args.dout + "/orthofinder"
        dout_cand_scgs = args.dout + "/candidate_scgs"
        dout_hmm_db = dout_cand_scgs + "/hmm_db"
        
        # define variables for output files
        fout_cand_scg_domtbl = args.dout + "/candidate_scg_domtbl.txt"
        fout_score_table = args.dout + "/score_table.csv"
        fout_scg_table = args.dout + "/scg_table.csv"
        
        # create directories
        os.makedirs(dout_cand_scgs, exist_ok = True)
        
        # 1) select seeds and construct gene families
        fins_genomes = read_paths(args.fin_genomepaths)
        fins_seeds = random.sample(fins_genomes, int(args.n_seed_genomes))
        run_orthofinder(fins_seeds, dout_orthofinder)
        
        # 2) identify candidate scgs
        dout_orthofinder_results = glob.glob(dout_orthofinder + "/Results_*")[0]
        taxon = read_orthogroups(dout_orthofinder_results)
        orthogroups = select_orthogroups(taxon, int(args.min_presence_in_seeds))
        
        # define file lists for candidate scgs
        dout_orthofinder_ogs = glob.glob(dout_orthofinder_results + "/Orthologues_*/Sequences")[0]
        fouts_cand_scg_seqs = [dout_orthofinder_ogs + "/" + orthogroup.name + ".fa"
            for orthogroup in orthogroups]
            
        # 3) scan orthogroups in all genomes
        construct_profile_db(fouts_cand_scg_seqs, dout_cand_scgs)
        print("unzipping all genomes")
        run_gunzip_parallel(fins_genomes)
        print("performing hmmer search")
        fins_genomes_unzipped = [str(pathlib.PurePosixPath(fin).with_suffix("")) 
            for fin in fins_genomes]
        run_hmmsearch(dout_hmm_db, fins_genomes_unzipped, fout_cand_scg_domtbl)
        print("rezipping all genomes")
        run_gzip_parallel(fins_genomes_unzipped)
        print("reading hmmer domtbl")
        domtbl = read_domtbl(fout_cand_scg_domtbl)
        print("creating score table")
        score_table = process_domtbl(domtbl, fins_genomes)
        score_table.to_csv(fout_score_table)
        print("constructing scg table")
        construct_scg_table(fout_score_table, fout_scg_table)
        
    elif task == "construct_genome_table":
        construct_genome_table(args.fin_score_table, args.fin_candidate_scg_table,
            args.candidate_scg_cutoff, args.fout_genome_table)
            
    elif task == "construct_scg_matrix":
        construct_scg_matrix(args.fin_score_table, args.fin_genome_table,
            args.fin_genome_list, args.candidate_scg_cutoff, 
            args.fout_scg_matrix)

    else:
    
        print("this task is not yet implemented")
        
    print("raccoon out")
    print("")
