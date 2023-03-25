import glob
import gzip
import os
import re
import pandas as pd

from Bio import SeqIO, AlignIO
from concurrent.futures import ProcessPoolExecutor
from io import StringIO

from utils import *
        
def read_fastapaths(path):
    # when path is a file
    if os.path.isfile(path):
        # store fastapaths in list
        fastapaths = [fp.strip() for fp in open(path)]
    # when path is a folder
    elif os.path.isdir(path):
        # store fastapaths in list
        fastapaths = [os.path.join(path, file) for file in os.listdir(path)]
        fastapaths = [fp for fp in fastapaths if os.path.isfile(fp)]
        extensions = ("fasta", "fa", "faa", "ffn", "fasta.gz", "fa.gz", 
            "faa.gz", "ffn.gz")
        fastapaths = [fp for fp in fastapaths if fp.endswith(extensions)]
    return(fastapaths)

def read_mmseqs_table(fin):
    names = ["query", "target", "pident", "alnlen", "mismatch", "gapopen", 
        "qstart", "qend", "tstart", "tend", "evalue", "bits"]
    table = pd.read_csv(fin, sep = "\t", names = names)
    return(table)

def read_fasta(fin):
    with open_smart(fin) as hin:
        seqs = list(SeqIO.parse(hin, "fasta"))
    return(seqs)
    
def write_fasta(seqs, fout):
    open(fout, "w").close()
    with open(fout, "a") as hout:
        for seq in seqs:
            SeqIO.write(seq, hout, "fasta")

def read_domtbl(fin_domtbl):
    domtbl = pd.read_csv(fin_domtbl, delim_whitespace = True, comment = '#',
        usecols = [0, 3, 13, 15, 16],
        names = ['gene', 'profile', 'score', 'hmm_from', 'hmm_to'],
        dtype = {'gene': str, 'profile': str, 'score': float, 'hmm_from': int,
            'hmm_to': int},
        engine = 'c'
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

def extract_genetable(fin_genome):
    "Extracts a gene table from an ffn/faa file."
    
    genome = filename_from_path(fin_genome)
    with open_smart(fin_genome) as hin:
        fasta = hin.read()
    regex = re.compile("^>([^\n\t ]+)", re.MULTILINE)
    genes_genome = regex.findall(fasta)
    genes_genome = pd.DataFrame({"gene": genes_genome})
    genes_genome.loc[:, "genome"] = genome
    return(genes_genome)

def extract_genes(fins_genomes, threads = 1):
    "Extracts a gene table from multiple ffn/faa files"
    
    with ProcessPoolExecutor(max_workers = threads) as executor:
        genes = executor.map(extract_genetable, fins_genomes)
    genes = pd.concat(genes)
    return(genes)
  
def fastas2stockholm(fins_alis, fout_sto):
    with open(fout_sto, "a") as hout_sto:
        for fin_ali in fins_alis:
            orthogroup = os.path.splitext(os.path.basename(fin_ali))[0]
            with open(fin_ali) as hin_ali:
                ali = AlignIO.read(hin_ali, "fasta")
                with StringIO() as vf:
                    AlignIO.write(ali, vf, "stockholm")
                    sto = vf.getvalue().splitlines()
                sto.insert(1, "#=GF ID " + orthogroup)
                for line in sto:
                    hout_sto.write(line + "\n")
