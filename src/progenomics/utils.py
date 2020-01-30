import os
import pandas as pd

from Bio import SeqRecord
from statistics import mean

def make_paths(filenames, folder, extension):
    paths = [os.path.join(folder, filename + extension) for filename in
        filenames]
    return(paths)

def filename_from_path(path):
    filename = os.path.basename(path)
    filename = os.path.splitext(filename)[0]
    return(filename)

"""
The table "hits" should have the columns "gene", "profile", "score" and
"positive".
"""
def train_cutoffs(hits):
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
