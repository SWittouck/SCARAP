import concurrent.futures
import os
import pandas as pd

from Bio import SeqIO, SeqRecord
from random import shuffle
from statistics import mean

from utils import *
    
def calc_pgo(genomes_fam1, genomes_fam2):
    """Calculates the proportion of genome overlap (pgo).
    
    Calculates the observed proportion of genomes present in both subfamilies
    of a gene family.
    
    Args:
        genomes_fam1: A list of genomes that the genes of family 1 belong 
            to.
        genomes_fam2: A list of genomes that the genes of family 2 belong 
            to.
            
    Returns:
        The observed pgo
    """
    n_unique_fam1 = len(set(genomes_fam1))
    n_unique_fam2 = len(set(genomes_fam2))
    n_unique_tot = len(set(genomes_fam1 + genomes_fam2))
    pgo = (n_unique_fam1 + n_unique_fam2 - n_unique_tot) / n_unique_tot
    return(pgo)
    
def ncat_exp(n, probs):
    """Calculates the expected number of distinct categories.
    
    Calculates the expected number of distinct categories present in a 
    multinomial sample. See Emigh, 1983, Biometrics. 
    http://www.jstor.org/stable/2531019
    
    Args:
        n (int): The sample size
        probs: A list of probabilities that define the multinomial 
            distribution.
            
    Returns:
        The expected number of distinct categories.
    """
    ncat_exp = len(probs) - sum([(1 - prob) ** n for prob in probs])
    return(ncat_exp)
    
def pred_pgo(genomes_fam1, genomes_fam2):
    """Predicts the proportion of genome overlap (pgo).
    
    Predicts the percentage of unique genomes that are expected to be present 
    in both subfamilies of a gene family, give a model where genomes are 
    assigned randomly to the genes in each subfamily.
    
    Args:
        genomes_fam1: A list of genomes that the genes of family 1 belong 
            to.
        genomes_fam2: A list of genomes that the genes of family 2 belong 
            to.
            
    Returns:
        The predicted pgo.
    """
    genomes = genomes_fam1 + genomes_fam2
    freqs = pd.Series(genomes).value_counts().tolist()
    tot = sum(freqs)
    probs = [freq / tot for freq in freqs]
    ncat_exp_tot = ncat_exp(len(genomes), probs)
    ncat_exp_fam1 = ncat_exp(len(genomes_fam1), probs)
    ncat_exp_fam2 = ncat_exp(len(genomes_fam2), probs)
    pgo = (ncat_exp_fam1 + ncat_exp_fam2 - ncat_exp_tot) / len(set(genomes))
    return(pgo)

def decide_split(genomes_fam1, genomes_fam2):
    """Decides if a gene family should be split or not.
    
    Args:
        genomes_fam1: A list of genomes that the genes of family 1 belong 
            to.
        genomes_fam2: A list of genomes that the genes of family 2 belong 
            to.
        
    Returns:
        A boolean value reflecting if the family should be split or not.
    """
    pgo_obs = calc_pgo(genomes_fam1, genomes_fam2)
    pgo_exp = pred_pgo(genomes_fam1, genomes_fam2)
    return(pgo_obs > pgo_exp)
    
def train_cutoffs(hits):
    """
    The table "hits" should have the columns "gene", "profile", "score" and
    "positive".
    """
    profile_to_cutoff = {}
    for profile, profile_hits in hits.groupby("profile"):
        scores_positive = profile_hits["score"][profile_hits["positive"]]
        scores_negative = profile_hits["score"][~ profile_hits["positive"]]
        if scores_negative.empty:
            profile_to_cutoff[profile] = mean([0, mean(scores_positive)])
        elif scores_positive.empty:
            profile_to_cutoff[profile] = 0
        else:
            profile_to_cutoff[profile] = mean([mean(scores_positive),
                mean(scores_negative)])
    profiles = pd.DataFrame(profile_to_cutoff.items(), columns = ["profile",
        "cutoff"])
    return(profiles)

def train_cutoffs_pan(hits, pangenome):
    hits = pd.merge(hits, pangenome[["gene", "orthogroup"]], how = "left")
    hits["positive"] = hits["profile"] == hits["orthogroup"]
    profiles = train_cutoffs(hits[["gene", "profile", "score", "positive"]])
    return(profiles)

def train_cutoffs_core(hits, genes_genomes):
    hits = pd.merge(hits, genes_genomes[["gene", "genome"]], how = "left")
    ixs_positive = (hits[["gene", "genome", "score"]].
        groupby(["gene", "genome"]).idxmax().score)
    hits["positive"] = False
    hits.loc[ixs_positive, "positive"] = True
    profiles = train_cutoffs(hits[["gene", "profile", "score", "positive"]])
    return(profiles)

def process_scores(hits, orthogroups):
    hits = hits.rename(columns = {"profile": "orthogroup"})
    hits = pd.merge(hits, orthogroups, how = "left")
    hits = hits[hits.score >= hits.cutoff]
    ixs_to_keep = hits[["gene", "score"]].groupby("gene").idxmax().score
    genes = hits.loc[ixs_to_keep, :][["gene", "orthogroup"]]
    return(genes)

def checkgenomes(coregenome):
    n_groups = len(set(coregenome.orthogroup))
    genomes = {"genome": [], "completeness": [], "redundancy": []}
    for genome, genes in coregenome.groupby("genome"):
        genomes["genome"].append(genome)
        genomes["completeness"].append(len(set(genes.orthogroup)) / n_groups)
        genomes["redundancy"].append(sum(genes.orthogroup.value_counts() > 1) /
            n_groups)
    genomes = pd.DataFrame(genomes)
    return(genomes)

def checkgroups(coregenome):
    n_genomes = len(set(coregenome.genome))
    orthogroups = {"orthogroup": [], "coreness": []}
    for orthogroup, genes in coregenome.groupby("orthogroup"):
        orthogroups["orthogroup"].append(orthogroup)
        orthogroups["coreness"].append(sum(genes.genome.value_counts() == 1) /
            n_genomes)
    orthogroups = pd.DataFrame(orthogroups)
    return(orthogroups)

def filter_genomes(pangenome, genomes):
    pangenome = pangenome[pangenome.genome.isin(genomes)]
    return(pangenome)

def filter_groups(pangenome, orthogroups):
    pangenome = pangenome[pangenome.orthogroup.isin(orthogroups)]
    return(pangenome)

def select_representative(records):
    records = sorted(records, key = lambda x: len(x))
    return(records[len(records) // 2])

"""
Function to align a SeqRecord of nucleotide sequences, given a SeqRecord
of aligned amino acid sequences.
Input:
  - og_nucs: SeqRecord of nucleotide sequences to align
  - og_aas_aligned = SeqRecord of amino acid sequences that are aligned
    (reference alignment)
Output:
  - SeqRecord of aligned nucleotide sequences
Naming:
  - sr = SeqRecord (Biopython)
  - og = orthogroup
  - nucs = nucleotides
  - aas = amino acids
"""
def reverse_align(og_nucs, og_aas_aligned):
    og_aas_aligned = {seq.id: seq.seq for seq in og_aas_aligned}
    og_nucs_aligned = []
    for nucs_sr in og_nucs:
        nucs_id = nucs_sr.id
        nucs_seq = nucs_sr.seq
        aas_gapped_seq = og_aas_aligned[nucs_id]
        nucs_gapped_seq = ""
        nucs_pos = 0
        for aa in aas_gapped_seq:
            if aa == "-":
                nucs_gapped_seq += "---"
            else:
                nucs_gapped_seq += nucs_seq[nucs_pos:(nucs_pos + 3)]
                nucs_pos += 3
        nucs_gapped_sr = SeqRecord.SeqRecord(id = nucs_id,
            seq = nucs_gapped_seq, description = "")
        og_nucs_aligned.append(nucs_gapped_sr)
    return(og_nucs_aligned)

def gather_orthogroup_sequences(pangenome, faapaths, dout_orthogroups,
    min_genomes = 1):
    
    # make output folder and empty if already exists
    makedirs_smart(dout_orthogroups)

    # filter orthogroups based on number of genomes they occur in
    if min_genomes > 1:
        pangenome = pangenome.groupby("orthogroup").\
            filter(lambda x: len(set(x["genome"])) >= min_genomes)

    # make dictionary to store faapaths of genomes
    genome_to_faapath = {filename_from_path(faapath): faapath for faapath in
        faapaths}

    # gather sequences and store in file per orthogroup
    for genome, genes in pangenome.groupby("genome"):
        faapath = genome_to_faapath[genome]
        gene_to_orthogroup = dict(zip(genes.gene, genes.orthogroup))
        with open_smart(faapath) as hin_genome:
            for record in SeqIO.parse(hin_genome, "fasta"):
                orthogroup = gene_to_orthogroup.get(record.id)
                if orthogroup is None:
                    continue
                fout_orthogroup = os.path.join(dout_orthogroups,
                    orthogroup + ".fasta")
                with open(fout_orthogroup, "a+") as hout_orthogroup:
                    SeqIO.write(record, hout_orthogroup, "fasta")

def collapse_pangenome(pangenome, faapathsfin, reprsfout, famprefix,
    tempdio):
    os.makedirs(tempdio)
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

"""
Function to align a nucleotide fasta using an amino acid fasta as a reference
alignment.
Arguments:
  - fin_og_nucs: nucleotide fasta to align (.fasta)
  - fin_og_aas_aligned: amino acid fasta that is aligned (.aln)
  - fout_og_nucs_aligned: file to store the aligned nucleotide fasta (.aln)
"""
def reverse_align_helper(fin_og_nucs, fin_og_aas_aligned, fout_og_nucs_aligned):
    og_nucs = list(SeqIO.parse(fin_og_nucs, "fasta"))
    og_aas_aligned = list(SeqIO.parse(fin_og_aas_aligned, "fasta"))
    og_nucs_aligned = reverse_align(og_nucs, og_aas_aligned)
    with open(fout_og_nucs_aligned, "w") as hout_og_nucs_aligned:
        SeqIO.write(og_nucs_aligned, hout_og_nucs_aligned, "fasta")

"""
Function to align a set of nucleotide fastas using a set of amino acid fastas
as reference alignments.
"""
def reverse_align_parallel(fins_ogs_nucs, fins_ogs_aas_aligned,
    fouts_ogs_nucs_aligned):

    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(reverse_align_helper, fins_ogs_nucs, fins_ogs_aas_aligned,
            fouts_ogs_nucs_aligned)
