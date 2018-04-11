#!/usr/bin/env python3

import argparse
import sys
import os
import glob
import subprocess
import concurrent.futures
import shutil
import re
from random import randrange
from Bio.Align.Applications import MafftCommandline

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


"""
Method to parse an "orthogroups.csv" file
Returns an object of class Taxon  
"""
def parse_orthogroups(din_orthofinder):
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
    return(taxon)

def select_orthogroups(orthogroups, min_genomes):
    orthogroups_selected = []
    for orthogroup in orthogroups:
        if orthogroup.get_n_genomes() >= min_genomes:
            orthogroups_selected.append(orthogroup)
    return(orthogroups_selected)

def align(fin_sequences, fout_alignment):
    mafft_cline = MafftCommandline(input = fin_sequences)
    stdout, stderr = mafft_cline()
    with open(fout_alignment, "w") as hout_alignment:
        hout_alignment.write(stdout)

def align_orthogroups(orthogroups, din_orthogroups, dout_alignments):
    os.makedirs(dout_alignments, exist_ok = True)
    n = len(orthogroups)
    fins_sequences = [din_orthogroups + "/" + orthogroup.name + ".fa"
        for orthogroup in orthogroups]
    fouts_alignments = [dout_alignments + "/" + orthogroup.name + ".fasta"
        for orthogroup in orthogroups]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(align, 
            fins_sequences, fouts_alignments
        )

def construct_profile(fin_alignment, fout_profile):
   subprocess.run(
        ["hmmbuild", fout_profile, fin_alignment],
        stdout = subprocess.PIPE,
    )
 
def construct_profile_db(orthogroups, din_alignments, dout_profiles, dout_profile_db):
    fins_alignments = [din_alignments + "/" + orthogroup.name + ".fasta"
        for orthogroup in orthogroups]
    fouts_profiles = [dout_profiles + "/" + orthogroup.name + ".hmm"
        for orthogroup in orthogroups]
    os.makedirs(dout_profiles, exist_ok = True)
    os.makedirs(dout_profile_db, exist_ok = True)
    n = len(orthogroups)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(construct_profile, 
            fins_alignments, fouts_profiles
        )
    fout_profile_db = dout_profile_db + "/profiles"
    subprocess.run(
        " ".join(["cat", dout_profiles + "/*.hmm", ">", fout_profile_db]), 
        shell = True
    )
    subprocess.run(
        ["hmmpress", fout_profile_db],
        stdout = subprocess.PIPE
    )

def run_hmmscan(fin_sequences, fin_profile_db, fout_scores):
    subprocess.run(
        ["hmmscan", "--domtblout", fout_scores, 
            fin_profile_db, fin_sequences],
        stdout = subprocess.PIPE
    )
   
def run_hmmscan_parallel(fins_sequences, fin_hmm_db, fout_scores):
    dout_scores_temp = os.path.dirname(fout_scores) + "/scores_temp"
    os.makedirs(dout_scores_temp, exist_ok = True)
    n = len(fins_sequences)
    fouts_scores_temp = ["%s/%i.tsv" % (dout_scores_temp, i) for i in range(n)]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_hmmscan, 
            fins_sequences, 
            [fin_hmm_db] * n, 
            fouts_scores_temp
        )
    for fout_scores_temp in fouts_scores_temp:
        subprocess.run(
            " ".join(["cat", fout_scores_temp, ">>", fout_scores]),
            shell = True
        )
        subprocess.run(" ".join(["rm", fout_scores_temp]), shell = True)
    subprocess.run(" ".join(["rm -r", dout_scores_temp]), shell = True)

def select_genomes_random(fins_genomes, n):
    fins_seeds = []
    for _ in range(n):
        pos = randrange(len(fins_genomes))
        fins_seeds.append(fins_genomes[pos])
        fins_genomes[pos] = fins_genomes[-1]
        del fins_genomes[-1]
    return(fins_seeds)

def run_orthofinder(fins_genomes, dout_orthofinder):
    os.makedirs(dout_orthofinder, exist_ok = True)
    for fin_genome in fins_genomes:
        genome = re.search("[^/]+$", fin_genome).group(0)
        shutil.copyfile(fin_genome, dout_orthofinder + "/" + genome)
        subprocess.run(["gunzip", dout_orthofinder + "/" + genome]) 
    command = ["orthofinder", "-M", "msa", "-os", "-t", "8",
        "-f", dout_orthofinder]
    subprocess.run(command)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "task",
        help = "the task you want to perform",
        choices = [
            "construct_profile_db",
            "get_core_genes"
        ]
    )
    parser.add_argument(
        "--din_orthofinder",
        help = "the input directory with orthofinder files"
    )
    parser.add_argument(
        "--dout_profile_db",
        help = "the output directory for the hmm profiles of the orthogroups"
    )
    parser.add_argument(
        "--min_genomes",
        help = "minimum number of genomes an orthogroup has to be present in"
    )
    args = parser.parse_args()
    return(args)

if __name__ == "__main__":

    print("\nthis is raccoon, version unknown")
    args = parse_arguments()

    if args.task == "construct_profile_db":
        if args.din_orthofinder is None:
            print("argument --din_orthofinder is required")
        elif args.dout_profile_db is None:
            print("argument --dout_profile_db is required")
        elif args.min_genomes is None:
            print("argument --min_genomes is required")
        else:
            print("starting profile construction on orthofinder output")
            taxon = parse_orthogroups(args.din_orthofinder)
            orthogroups = select_orthogroups(taxon.orthogroups, int(args.min_genomes))
            print("found %i orthogroups" % (len(orthogroups)))
            d_alignments = args.dout_profile_db + "/alignments"
            d_profiles = args.dout_profile_db + "/profiles"
            d_profile_db = args.dout_profile_db + "/hmmer_db"
            din_orthofinder = args.din_orthofinder + "/Orthologues_*/Sequences"
            din_orthofinder = glob.glob(din_orthofinder)[0]
            align_orthogroups(orthogroups, din_orthofinder, d_alignments)
            construct_profile_db(orthogroups, d_alignments, d_profiles, d_profile_db)

    elif args.task == "get_core_genes":
        print("this task is not yet implemented")

    else:
        print("this task is not yet implemented")

    print("")
