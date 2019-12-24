import glob
import logging
import os
import pandas as pd
import re
import shutil

from Bio import SeqIO, SeqRecord
from pathlib import Path
from utils import *

def gather_orthogroup_sequences(pangenome, faapaths, dout_orthogroups,
    min_genomes = 1):

    # filter orthogroups based on number of genomes they occur in
    if min_genomes > 1:
        pangenome = pangenome.groupby("orthogroup").\
            filter(lambda x: len(set(x["genome"])) >= min_genomes)

    # make dictionary to store faapaths of genomes
    genome_to_faapath = {Path(faapath).stem: faapath for faapath in
        faapaths}

    # gather sequences and store in file per orthogroup
    for genome, genes in pangenome.groupby("genome"):
        faapath = genome_to_faapath[genome]
        gene_to_orthogroup = dict(zip(genes.gene, genes.orthogroup))
        with open(faapath, "r") as hin_genome:
            for record in SeqIO.parse(hin_genome, "fasta"):
                orthogroup = gene_to_orthogroup.get(record.id)
                if orthogroup is None:
                    continue
                fout_orthogroup = os.path.join(dout_orthogroups,
                    orthogroup + ".fasta")
                with open(fout_orthogroup, "a+") as hout_orthogroup:
                    SeqIO.write(record, hout_orthogroup, "fasta")

def read_domtbl(fin_domtbl):
    domtbl = pd.read_csv(fin_domtbl, delim_whitespace = True, comment = '#',
        usecols = [0, 3, 13, 15, 16],
        names = ['gene', 'profile', 'score', 'hmm_from', 'hmm_to'],
        dtype = {'gene': str, 'profile': str, 'score': float, 'hmm_from': int,
            'hmm_to': int},
        engine = 'python'
    )
    return(domtbl)

def read_genes(fin):
    genes = pd.read_csv(fin, sep = "\t",
        names = ["gene", "genome", "orthogroup"])
    return(genes)

def read_orthogroups(fin):
    orthogroups = pd.read_csv(fin, sep = "\t", names = ["orthogroup", "cutoff"])
    return(orthogroups)

def read_species(fin):
    genomes = pd.read_csv(fin, sep = "\t", names = ["genome", "species"])
    return(genomes)

def write_tsv(genes, fout):
    genes.to_csv(fout, sep = "\t", index = False, header = False)

def write_lines(lines, fout):
    lines = pd.DataFrame({"line": lines})
    lines.to_csv(fout, index = False, header = False)

def read_lines(fin):
    with open(fin) as hin:
        lines = [line.strip() for line in hin.readlines()]
    return(lines)

def read_orthogroups_orthofinder(fin_orthogroups):
    with open(fin_orthogroups, 'r') as hin_orthogroups:
        pangenome = {"gene": [], "genome": [], "orthogroup": []}
        header = hin_orthogroups.readline().strip()
        genomes = header.split("\t")[1:]
        for line in hin_orthogroups.readlines():
            line = line.strip().split("\t")
            orthogroup = line[0]
            for genome_ix, copies in enumerate(line[1:]):
                if copies != "":
                    for copy in copies.split(", "):
                        pangenome["gene"].append(copy)
                        pangenome["genome"].append(genomes[genome_ix])
                        pangenome["orthogroup"].append(orthogroup)
    pangenome = pd.DataFrame.from_dict(pangenome)
    return(pangenome)

def read_pangenome_orthofinder(din_orthofinder):
    din_pangenome = os.path.join(din_orthofinder,
        "OrthoFinder/Results_*/Orthogroups")
    din_pangenome = glob.glob(din_pangenome)[0]
    fin_genes_1 = os.path.join(din_pangenome, "Orthogroups.tsv")
    genes_assigned = read_orthogroups_orthofinder(fin_genes_1)
    fin_genes_2 = os.path.join(din_pangenome,
        "Orthogroups_UnassignedGenes.tsv")
    genes_unassigned = read_orthogroups_orthofinder(fin_genes_2)
    pangenome = genes_assigned.append(genes_unassigned)
    return(pangenome)

def extract_genes(fins_genomes):

    genes = pd.DataFrame()

    for fin_genome in fins_genomes:
        genome = filename_from_path(fin_genome)
        with open(fin_genome, mode = "rt") as hin:
            fasta = hin.read()
        genes_genome = re.compile(">([^ ]+)").findall(fasta)
        genes_genome = pd.DataFrame({"gene": genes_genome})
        try:
            genes_genome.loc[:, "genome"] = genome
            genes = genes.append(genes_genome)
        except ValueError:
            pass

    return(genes)

def collapse_pangenome(pangenomefin, faapathsfin, reprsfout, famprefix,
    tempdio):
    os.makedirs(tempdio)
    pangenome = read_genes(pangenomefin)
    faapaths = read_lines(faapathsfin)
    gather_orthogroup_sequences(pangenome, faapaths, tempdio)
    famfins = [os.path.join(tempdio, file) for file in os.listdir(tempdio)]
    reprs = []
    for famfin in famfins:
        fam = filename_from_path(famfin)
        with open(famfin, "r") as famhin:
            records = list(SeqIO.parse(famhin, "fasta"))
            repr = select_representative(records)
            repr.id = famprefix + "-" + fam
            reprs.append(repr)
    with open(reprsfout, "w") as reprshout:
        SeqIO.write(reprs, reprshout, "fasta")
    shutil.rmtree(tempdio)

"""
Function to construct a supermatrix given alignments of individual single-copy
core orthogroups. If there are two or more copies of an orthogroup in a genome,
all copies will be dropped.
Input:
  - coregenome: DataFrame with columns gene, genome, orthogroup
  - alifins: input files with alignments of single-copy core orthogroups
  - supermatrixfout: fasta file containing the supermatrix alignment
"""
def construct_supermatrix(coregenome, alifins, supermatrixfout):

    supermatrix = {}
    genomes = list(set(coregenome.genome))
    n_genomes = len(genomes)
    for genome in genomes:
        supermatrix[genome] = SeqRecord.SeqRecord(id = genome, seq = "",
            description = "")

    alifindict = {filename_from_path(alifin): alifin for alifin in alifins}

    n_fams_sc = 0

    for orthogroup, rows in coregenome.groupby("orthogroup"):
        alifin = alifindict[orthogroup]
        sequencedict = {}
        for record in SeqIO.parse(alifin, "fasta"):
            alilen = len(record.seq)
            sequencedict[record.id] = record.seq
        rows = rows.drop_duplicates("genome", keep = False)
        rows = pd.merge(pd.DataFrame({"genome": genomes}), rows, how = "left")
        for ix, row in rows.iterrows():
            sequence_to_add = sequencedict.get(row.gene, "-" * alilen)
            supermatrix[row.genome] = supermatrix[row.genome] + sequence_to_add

    with open(supermatrixfout, "a") as supermatrixhout:
        for genome in supermatrix:
            SeqIO.write(supermatrix[genome], supermatrixhout, "fasta")
