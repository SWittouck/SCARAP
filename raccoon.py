#!/usr/bin/env python3

import argparse
import sys
import os
import glob
import subprocess
import concurrent.futures
from Bio.Align.Applications import MafftCommandline

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


class Genome:

    def __init__(self, name):
        self.name = name
        self.genes = []

    def __str__(self):
        return(self.name)

    def add_gene(self, gene):
        gene.genome = self


class Orthogroup:

    def __init__(self, name):
        self.name = name
        self.genes = []

    def __str__(self):
        return(self.name)

    def add_gene(self, gene):
        gene.orthogroup = self


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
def parse_orthofinder(din_orthofinder):
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

def select_orthogroups(orthogroups, n_genes_min):
    orthogroups_selected = []
    for orthogroup in orthogroups:
        if len(orthogroup.genes) >= n_genes_min:
            orthogroups_selected.append(orthogroup)
    return(orthogroups_selected)

def align_orthogroup(orthogroup, din_orthogroups, dout_alignments):
    fin_orthogroup = din_orthogroups + "/" + orthogroup.name + ".fa"
    fout_alignment = dout_alignments + "/" + orthogroup.name + ".fasta"
    mafft_cline = MafftCommandline(input = fin_orthogroup)
    stdout, stderr = mafft_cline()
    with open(fout_alignment, "w") as hout_alignment:
        hout_alignment.write(stdout)

def align_orthogroups(orthogroups, din_orthogroups, dout_alignments):
    os.makedirs(dout_alignments, exist_ok = True)
    n = len(orthogroups)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(align_orthogroup, 
            orthogroups, [din_orthogroups] * n, [dout_alignments] * n
        )

def construct_profile(orthogroup, din_alignments, dout_profiles):
    fin_alignment = din_alignments + "/" + orthogroup.name + ".fasta"
    fout_profile = dout_profiles + "/" + orthogroup.name + ".hmm"
    subprocess.run(
        ["hmmbuild", fout_profile, fin_alignment],
        stdout = subprocess.PIPE,
    )
 
def construct_hmmer_db(orthogroups, din_alignments, dout_profiles, dout_hmmer_db):
    os.makedirs(dout_profiles, exist_ok = True)
    os.makedirs(dout_hmmer_db, exist_ok = True)
    n = len(orthogroups)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(construct_profile, 
            orthogroups, [din_alignments] * n, [dout_profiles] * n
        )
    fout_hmmer_db = dout_hmmer_db + "/profiles"
    subprocess.run(
        " ".join(["cat", dout_profiles + "/*.hmm", ">", fout_hmmer_db]), 
        shell = True
    )
    subprocess.run(
        ["hmmpress", fout_hmmer_db],
        stdout = subprocess.PIPE
    )

def run_hmmer(fin_sequences, fin_hmm_db, fout_scores):
    subprocess.run(
        ["hmmscan", "--domtblout", fout_scores, 
            fin_hmm_db, fin_sequences],
        stdout = subprocess.PIPE
    )
   
def perform_hmmscan(fins_sequences, fin_hmm_db, fout_scores):
    dout_scores_temp = os.path.dirname(fout_scores) + "/scores_temp"
    os.mkdir(dout_scores_temp)
    n = len(fins_sequences)
    fouts_scores_temp = ["%s/%i.tsv" % (dout_scores_temp, i) for i in range(n)]
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_hmmer, 
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
    # to do: savely remove dout_scores_temp

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "task",
        help = "the task you want to perform",
        choices = [
            "select_seeds",
            "construct_orthogroups_seeds",
            "construct_profiles_seeds",
            "scan_seeds",
            "scan_genomes",
            "collect_core_fastas"
        ]
    )
    parser.add_argument(
        "--din_orthofinder_seeds",
        help = "the input directory with orthofinder files"
    )
    parser.add_argument(
        "--dout_profiles_seeds",
        help = "the output directory for the hmm profiles of the seed orthogroups"
    )
    parser.add_argument(
        "--din_queries",
        help = "the input directory with fasta file(s) containing query sequences"
    )
    parser.add_argument(
        "--din_hmmer_db",
        help = "the input directory with hmmer database files (.h3f, .h3i, .h3m and .h3p)"
    )
    parser.add_argument(
        "--fout_scores",
        help = "the output file to store hmmer scores"
    )
    args = parser.parse_args()
    return(args)

if __name__ == "__main__":
    print("\nthis is raccoon, version unknown")
    args = parse_arguments()
    if args.task == "select_seeds":
        print("seed selection is not yet implemented")
    elif args.task == "construct_orthogroups_seeds":
        print("construction of seed orthogroups is not yet implemented")
    elif args.task == "construct_profiles_seeds":
        if args.din_orthofinder_seeds is None:
            print("argument --din_orthofinder_seeds is required")
        elif args.dout_profiles_seeds is None:
            print("argument --dout_profiles_seeds is required")
        else:
            print("starting profile construction on seeds")
            taxon = parse_orthofinder(args.din_orthofinder_seeds)
            orthogroups = select_orthogroups(taxon.orthogroups, 10)
            print("found %i orthogroups" % (len(orthogroups)))
            d_alignments = args.dout_profiles_seeds + "/alignments"
            d_profiles = args.dout_profiles_seeds + "/profiles"
            d_hmmer_db = args.dout_profiles_seeds + "/hmmer_db"
            din_orthogroups = args.din_orthofinder_seeds + "/Orthologues_*/Sequences"
            din_orthogroups = glob.glob(din_orthogroups)[0]
            align_orthogroups(orthogroups, din_orthogroups, d_alignments)
            construct_hmmer_db(orthogroups, d_alignments, d_profiles, d_hmmer_db)
            fins_orthofinder_seeds = [din_orthogroups + "/" + 
                orthogroup.name + ".fa" for orthogroup in orthogroups]
            # f_seed_scores = args.dout_profiles_seeds + "/seed_scores.tsv"
            # perform_hmmscan(fins_orthofinder_seeds, d_hmmer_db + "/profiles", f_seed_scores)
    elif args.task == "scan_genomes":
        if args.din_queries is None:
            print("argument --din_queries is required")
        elif args.din_hmmer_db is None:
            print("argument --din_hmmer_db is required")
        elif args.fout_scores is None:
            print("argument --fout_scores is required")
        else:
            fins_queries = [os.path.join(args.din_queries, f) 
                for f in os.listdir(args.din_queries)]
            perform_hmmscan(fins_queries, args.din_hmmer_db, args.fout_scores)
    else:
        print("this task is not yet implemented")
    print("")
