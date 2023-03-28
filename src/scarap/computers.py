import concurrent.futures
import logging
import numpy as np
import os
import pandas as pd

from Bio import SeqIO, SeqRecord
from random import shuffle
from statistics import mean

from readerswriters import *
from utils import *

def subset_idmat(matrix, rownames, rownames_sub):
    """Subsets an identity matrix. 
    
    Args:
        matrix: An identity matrix as a numpy array
        rownames: A list of rownames (= colnames)
        rownames_sub: A list of rownames to select
        
    Returns:
        A subset of the identity matrix
        A subset of the rownames
    """
    ixs = [ix for ix, gene in enumerate(rownames) if gene in rownames_sub]
    ixs = np.array(ixs)
    matrix_sub = matrix[ixs[:, None], ixs]
    rownames_sub = [rownames[ix] for ix in ixs]
    return(matrix_sub, rownames_sub)

def distmat_from_idmat(idmat):
    """Calculates a distance matrix.
    
    Args:
        idmat: An identity matrix as a numpy array
    
    Returns:
        A distance matrix in the format needed by the 
            cluster.hierarchy.linkage function of scipy.
    """
    n = np.size(idmat, 0)
    dm = []
    for r in range(n - 1):
        for c in range(r + 1, n):
            dm.append(1 - idmat[r, c])
    return(dm) 
    
def identity_matrix(seqrecords):
    '''Calculates a pairwise sequence identity matrix.
    
    Probably approximately as fast as the distmat tool from EMBOSS, but this 
    function can easily be parallelized on (batches of) columns.
    
    Args:
        seqrecords (list): A list of DNA sequences or SeqRecords.
        
    Returns:
        A matrix with sequence identity values. 
    '''

    n_seqs = len(seqrecords)
    n_cols = len(seqrecords[0].seq)

    # matrix with pairwise counts of alignment positions that are
    # identical
    identity_counts = np.zeros((n_seqs, n_seqs), np.uint32)
    # matrix with pairwise counts of alignment positions that are
    # incomparable
    comparable_counts = np.zeros((n_seqs, n_seqs), np.uint32)

    for col in range(n_cols):
        # print("col ", col, " of ", n_cols)
        # dictionary to store for each possible character the indices
        # of sequences with that character
        chars = {}
        for sr_ix, sr in enumerate(seqrecords):
            char = sr.seq[col]
            if char != "-":
                chars.setdefault(char, []).append(sr_ix)
        # for letters in the alphabet, add identity counts
        ixs_comp = np.array([], np.uint32)
        for char in chars.keys():
            ixs = np.array(chars[char], np.uint32)
            identity_counts[ixs[:, None], ixs] += 1
            ixs_comp = np.hstack((ixs_comp, ixs))
        # for all sequence combinations that both have a letter in the alphabet
        # add that they are comparable
        comparable_counts[ixs_comp[:, None], ixs_comp] += 1
        
    # avoid division by zero
    comparable_counts[comparable_counts == 0] += 1

    identity_matrix = identity_counts / comparable_counts

    return(identity_matrix)

def split_possible(genomes):
    """Determines whether splitting the family is an option. 
    
    Args:
        genomes (list): For each gene in the family, the genome that it belongs
            to.
            
    Returns:
        true/false
    """
    n_genes = len(genomes)
    n_genomes = len(set(genomes))
    if n_genes <= 3: return(False)
    if n_genomes == 1: return(False)
    if n_genes == n_genomes: return(False)
    return(True)

def calc_pgo(genomes_fam1, genomes_fam2):
    """Calculates the proportion of genome overlap (pgo).
    
    Calculates the observed proportion of genomes present in both subfamilies
    of a gene family.
    
    Args:
        genomes_fam1 (list): For each gene in family 1, the genome that it 
            belongs to.
        genomes_fam2 (list): For each gene in family 2, the genome that it 
            belongs to.
            
    Returns:
        The observed pgo.
    """
    if 0 in [len(genomes_fam1), len(genomes_fam2)]: return 0
    n_unique_fam1 = len(set(genomes_fam1))
    n_unique_fam2 = len(set(genomes_fam2))
    n_unique_tot = len(set(genomes_fam1 + genomes_fam2))
    pgo = (n_unique_fam1 + n_unique_fam2 - n_unique_tot) / \
        min([n_unique_fam1, n_unique_fam2])
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
    if 0 in [len(genomes_fam1), len(genomes_fam2)]: return 0
    genomes = genomes_fam1 + genomes_fam2
    freqs = pd.Series(genomes).value_counts().tolist()
    tot = sum(freqs)
    probs = [freq / tot for freq in freqs]
    ncat_exp_tot = ncat_exp(len(genomes), probs)
    ncat_exp_fam1 = ncat_exp(len(genomes_fam1), probs)
    ncat_exp_fam2 = ncat_exp(len(genomes_fam2), probs)
    pgo = (ncat_exp_fam1 + ncat_exp_fam2 - ncat_exp_tot) / \
        min([ncat_exp_fam1, ncat_exp_fam2])
    if pgo > 1: pgo = 1
    return(pgo)

def train_cutoffs(hits, pangenome):
    """Trains a cutoff per profile (orthogroup).
    
    Args:
        hits (DataFrame): A table with the columns gene, profile, score and 
            positive.
        pangenome (DataFrame): A table with the columns gene, genome and 
            orthogroup. 
            
    Returns:
        A table with the columns profile and cutoff.
    """
    hits = pd.merge(hits, pangenome[["gene", "orthogroup"]], how = "left")
    hits["positive"] = hits["profile"] == hits["orthogroup"]
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

def process_scores(hits, cutoffs, top_profiles = True):
    hits = pd.merge(hits, cutoffs, how = "left")
    hits = hits[hits.score >= hits.cutoff]
    if top_profiles:
        hits = hits.sort_values("score").\
            drop_duplicates(["gene"], keep = "last")
    hits = hits[["gene", "profile"]]
    hits = hits.rename(columns = {"profile": "orthogroup"})
    return(hits)

def determine_corefams(pan, core_filter, max_cores = 0):
    fams = checkgroups(pan)
    corefams = fams[fams.occurrence_singlecopy >= core_filter]
    if (max_cores > 0):
        corefams = corefams\
            .sort_values("occurrence_singlecopy", ascending = False)\
            .head(max_cores)
    return(corefams.orthogroup.tolist())

def checkgenomes(genes_core):
    """Computes genome statistics from a core genome table.
    
    Args:
        genes_core (Data Frame): A table with columns gene, genome and 
            orthogroup.
        
    Returns:
        A pandas Data Frame with the columns genome, completeness and 
            contamination.
    """
    n_groups = len(set(genes_core.orthogroup))
    genomes = genes_core\
        .groupby(["genome", "orthogroup"])\
        .agg(n_copies = ("gene", len))\
        .reset_index()\
        .groupby("genome")\
        .agg(completeness = ("n_copies", lambda x: len(x) / n_groups),
            contamination = ("n_copies", lambda x: sum(x > 1) / n_groups))\
        .reset_index()
    return(genomes)

def checkgroups(genes):
    """Computes orthogroup statistics from a pangenome table.
    
    Args:
        genes (Data Frame): A table with columns gene, genome and orthogroup.
        
    Returns:
        A pandas Data Frame with the columns orthogroup, occurrence and 
            occurrence_singlecopy.
    """
    n_genomes = len(set(genes.genome))
    orthogroups = genes\
        .groupby(["genome", "orthogroup"])\
        .agg(n_copies = ("gene", len))\
        .reset_index()\
        .groupby("orthogroup")\
        .agg(occurrence = ("n_copies", len), 
            occurrence_singlecopy = ("n_copies", lambda x: sum(x == 1)))\
        .apply(lambda x: x / n_genomes)\
        .reset_index()
    return(orthogroups)

def filter_genomes(pangenome, genomes):
    pangenome = pangenome[pangenome.genome.isin(genomes)]
    return(pangenome)

def filter_groups(pangenome, orthogroups):
    pangenome = pangenome[pangenome.orthogroup.isin(orthogroups)]
    return(pangenome)

def select_representative(records, longest = False):
    records = sorted(records, key = lambda x: len(x))
    if longest:
        return(records[-1])
    else: 
        return(records[len(records) // 2])


def reverse_align(og_nucs, og_aas_aligned):
    """Aligns a list of nucleotide sequences. 
    
    This function aligns a list of nucleotide sequences, given a list of 
        aligned amino acid sequences.
    
    Naming:
    
    * sr = SeqRecord (Biopython)
    * og = orthogroup
    * nucs = nucleotides
    * aas = amino acids
    
    Ags: 
        og_nucs (list): The SeqRecord list of nucleotide sequences to align.
        og_aas_aligned (list): The SeqRecord list of amino acid sequences that 
            are already aligned (reference alignment).
            
    Returns: 
        A list with seqrecords of aligned nucleotide sequences.
    """
    og_aas_aligned = {seq.id: seq.seq for seq in og_aas_aligned}
    og_nucs_aligned = []
    for nucs_sr in og_nucs:
        nucs_id = nucs_sr.id
        nucs_seq = nucs_sr.seq
        if not nucs_id in og_aas_aligned.keys():
            logging.warning(f"no translation found for gene {nucs_id}")
            continue
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
        pangenome = pangenome\
            .groupby("orthogroup")\
            .filter(lambda x: len(set(x["genome"])) >= min_genomes)\
            .reset_index()

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
    
def create_pseudogenome(pangenome, faapaths, tmpdio):
    "Returns a pseudogenome with one representative gene per orthogroup."
    os.makedirs(tmpdio)
    gather_orthogroup_sequences(pangenome, faapaths, tmpdio)
    genes = [None] * pangenome["orthogroup"].nunique()
    for ix, ogfin in enumerate(listpaths(tmpdio)):
        seqs = read_fasta(ogfin)
        repr = select_representative(seqs)
        genes[ix] = repr.id
    genetbl = pangenome[pangenome["gene"].isin(genes)]
    genetbl = genetbl[["gene", "genome"]] 
    shutil.rmtree(tmpdio)
    return(genetbl)

def construct_supermatrix(coregenome, alifins, supermatrixfout):
    """Constructs a contatenated alignment (= supermatrix).
    
    Function to construct a supermatrix given alignments of individual single-copy
    core orthogroups. If there are two or more copies of an orthogroup in a genome,
    all copies will be dropped.
    
    Args:
        coregenome (DataFrame): A table with the columns gene, genome and 
            orthogroup.
        alifins (list): A list of input files containing alignments of 
            single-copy core orthogroups.
        supermatrixfout (str): A path to a fasta file to write the supermatrix 
            to.
    """

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
