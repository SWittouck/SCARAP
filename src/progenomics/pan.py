import collections
import logging
import math
import numpy as np
import os
import shutil

from argparse import Namespace
from Bio import AlignIO, Align, SeqIO
from copy import copy
from ete3 import Tree
from concurrent.futures import ProcessPoolExecutor
from scipy import cluster

from utils import *
from readerswriters import *
from computers import *
from callers import *

## strategy-independent splitting helpers
        
def generate_thresholds(f = 2):
    thresholds = [0.01]
    while thresholds[-1] < 0.50:
        thresholds.append(thresholds[-1] * f)
    thresholds.reverse()
    thresholds[0] = 0.50
    thresholds = [1 - t for t in thresholds]
    return(thresholds)
    
def split_pan(pan, tree, outgr):
    """Split pan into pan1 and pan2 based on a tree and an outgroup node. 
    
    Args:
        pan (DataFrame): A gene table with at least the columns reprf and 
            orthogroup.
        tree: An ete3 tree (= the root node of a tree)
        outgr: A node in the tree
        
    Returns:
        [pan1, pan2]
    """
    reps_subfam1 = outgr.get_leaf_names()
    parent = outgr.up
    outgr.detach()
    reps_subfam2 = tree.get_leaf_names()
    parent.add_child(outgr) # re-add the outgroup
    pan1 = pan[pan["rep"].isin(reps_subfam1)].copy()
    pan2 = pan[pan["rep"].isin(reps_subfam2)].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    return([pan1, pan2])

def lowest_cn_roots(tree, pan):
    """Determine the set of lowest copy-number roots.
    
    Args:
        tree: ete3 tree object where the leaf names correspond to the values of
            the reprf column in pan.
        pan (DataFrame): Table with at least the columns reprf and genome.
        
    Returns:
        A list with references to the nodes that would be minimal copy number 
            roots.
    """

    # initialize list with reprfs and corresponding genomes 
    reprfs_genomes = {}
    for index, row in pan.iterrows():
        reprfs_genomes.setdefault(row["reprf"], []).append(row["genome"])
    reprfs_genomes = list(map(list, reprfs_genomes.items()))
    
    roots = []
    min_av_cn = 100000000
    
    # loop over all nodes except the root
    for node in tree.iter_descendants():
        # initialize empty genome lists for partition 1 and 2
        genomes1 = []
        genomes2 = []
        # loop over list with reprfs and corresponding genomes
        for element in reprfs_genomes:
            reprf = element[0]
            genomes = element[1]
            # append genomes of reprf to correct partition
            if reprf in node:
                genomes1.extend(genomes)
            else:
                genomes2.extend(genomes)
        cn1 = collections.Counter(genomes1).values()
        cn2 = collections.Counter(genomes2).values()
        av_cn = (sum(cn1) + sum(cn2)) / (len(cn1) + len(cn2))
        if av_cn < min_av_cn:
            roots = [node]
            min_av_cn = av_cn
        elif av_cn == min_av_cn:
            roots.append(node)
            
    return(roots)
    
def partition_genomes(reprfs_genomes, node):
    """Split genomes into two partitions using a tree. 
    
    Split genomes into two partitions based on their presence/absence in a 
        phylogenetic tree. 
    
    Args: 
        reprfs_genomes (list): List where each element is a list with two 
            elements: a representative gene (reprf) and a list with the
            genomes of all genes represented by reprf.
        node: An ete3 node (= tree) object.
        
    Returns: 
        A list with the genomes of the genes in the first partition.
        A list with the gneomes of the genes in the second partition. 
    """
    # initialize empty genome lists for partition 1 and 2
    genomes1 = []
    genomes2 = []
    # loop over list with reprfs and corresponding genomes
    for element in reprfs_genomes:
        reprf = element[0]
        genomes = element[1]
        # append genomes of reprf to correct partition
        if reprf in node:
            genomes1.extend(genomes)
        else:
            genomes2.extend(genomes)
    return(genomes1, genomes2)
    
def correct_root(root, tree, pan):
    """Correct the root of a phylogenetic tree.
    
    Correct the root of a tree by selecting the root that shows the lowest
    average copy number across both sides of the corresponding bipartion,
    while satisfying the restriction that all genomes present in both sides of
    the original bipartion (= genome overlap) are also in the genome overlap
    of the selected node.
    
    Args: 
        root: An ete3 node (= tree) that is present in tree.
        tree: An ete3 node (= tree).
        pan: A table with at least the columns reprf (corresponding to the tips
            of tree) and genome.
            
    Returns:
        An ete3 node representing an outgroup to the corrected root. 
    """

    # initialize list with reprfs and corresponding genomes 
    reprfs_genomes = {}
    for index, row in pan.iterrows():
        reprfs_genomes.setdefault(row["reprf"], []).append(row["genome"])
    reprfs_genomes = list(map(list, reprfs_genomes.items()))
    
    genomes1, genomes2 = partition_genomes(reprfs_genomes, root)
    midpoint_overlap = set(genomes1) & set(genomes2) # intersection
    
    roots = []
    min_av_cn = 100000000
    
    # loop over all nodes except the root
    for node in tree.iter_descendants():
        genomes1, genomes2 = partition_genomes(reprfs_genomes, node)
        overlap = set(genomes1) & set(genomes2) # intersection
        # if genomes that overlap in the midpoint bipartition do not all
        # overlap in this partition --> split
        if not midpoint_overlap.issubset(overlap):
            continue
        cn1 = collections.Counter(genomes1).values()
        cn2 = collections.Counter(genomes2).values()
        av_cn = (sum(cn1) + sum(cn2)) / (len(cn1) + len(cn2))
        if av_cn < min_av_cn:
            roots = [node]
            min_av_cn = av_cn
            # logging.info(f"new min av cn: {min_av_cn}")
        elif av_cn == min_av_cn:
            roots.append(node)  
            
    fam = pan.orthogroup[0]
  
    if root in roots:
        cnoutgr = root
        # logging.info(f"{fam}: cn root is midpoint root")
    else:
        cnoutgr = roots[0]
        # logging.info(f"{fam}: cn root is not midpoint root")
        
    return(cnoutgr)
    
def assess_split(pan1, pan2, family):
    """Assess whether a family should be split into two subfamilies.
    
    Remark: technically, this function only needs [genomes1, genomes2]. 
    However, extracting these from [pan1, pan2] is included here to avoid
    repeating this code in each split_family_recursive_STR function. The family
    argument is needed for logging.
    
    Args:
        pan1 (DataFrame): Gene table for the first subfamily, with the column
            genome.
        pan2 (DataFrame): Gene table for the second subfamily, with the column
            genome.
        family (str): Name of the gene family.
    
    Returns:
        True/False value indicating whether the split should happen.
    """
    genomes1 = pan1["genome"].tolist()
    genomes2 = pan2["genome"].tolist()
    pgo_obs = calc_pgo(genomes1, genomes2)
    pgo_exp = pred_pgo(genomes1, genomes2)
    split = pgo_obs >= pgo_exp and not pgo_exp == 0
    n1 = len(set(genomes1))
    n2 = len(set(genomes2))
    # if split:
    #     logging.info(f"{family}: {n1}/{n2}; "
    #         f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
    #         f"decision = split")
    # else:
    #     logging.info(f"{family}: {n1}/{n2}; "
    #         f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
    #         f"decision = no split")
    return(split)
    
## strategy-specific splitting helpers

def update_linreps(pan, thresholds, sequences, dio_tmp, threads):
    """Updates current and next linclust-based representatives. 
    
    Helper function for the L-based splitting strategies.
    
    Args:
        pan (DataFrame): A pangenome table with at least a "gene" column.
        thresholds (list): A list of percentage identity thresholds to try, in
            ascending order.
        sequences (list): A list of SeqRecords.
        dio_tmp (str): A folder to store temporary files.
        threads (int): Number of threads to use.
        
    Returns:
        A pan table with updated columns linrepc and linrepn.
        An updated list of thresholds where the thresholds that have already
            been tried are removed.
    """
    
    # prepare stuff
    thresholds = thresholds.copy()
    write_fasta(sequences, f"{dio_tmp}/seqs.fasta")        
    for dir in ["sequenceDB", "logs", "tmp"]:
        makedirs_smart(f"{dio_tmp}/{dir}")
    run_mmseqs(["createdb", f"{dio_tmp}/seqs.fasta", 
        f"{dio_tmp}/sequenceDB/db"], f"{dio_tmp}/logs/createdb.log", 
        threads = threads)
    
    # count number of linrepcs
    linrepcs = pan["linrepc"].unique().tolist()
    n_linrepcs = len(linrepcs) # linreps for current threshold
      
    # perform linclust until n_linrepcs > 1 or no thresholds left
    while n_linrepcs == 1 and thresholds:
        
        # previous linrepns become the linrepcs
        pan.linrepc = pan.linrepn
        n_linrepcs = len(linrepcs)
        
        # remove the previous linrepns
        pan = pan.drop("linrepn", axis = 1, errors = "ignore")
        
        # perform linclust on the next threshold
        id = str(thresholds[0])
        del thresholds[0]
        makedirs_smart(f"{dio_tmp}/clusterDB")
        run_mmseqs(["linclust", f"{dio_tmp}/sequenceDB/db", 
            f"{dio_tmp}/clusterDB/db", f"{dio_tmp}/tmp", "--min-seq-id", 
            id], f"{dio_tmp}/logs/linclust_{id}.log", threads = threads)
        run_mmseqs(["createtsv", f"{dio_tmp}/sequenceDB/db", 
            f"{dio_tmp}/sequenceDB/db", f"{dio_tmp}/clusterDB/db", 
            f"{dio_tmp}/clusters.tsv"], f"{dio_tmp}/logs/createtsv.log")
        genes_clusters = pd.read_csv(f"{dio_tmp}/clusters.tsv", sep = "\t",
            names = ["linrepn", "gene"])
                
        # add new linrepns to pan table
        pan = pan.reset_index().merge(genes_clusters, on = "gene")
        pan = pan.set_index("gene")
            
    return(pan, thresholds)     
    
def select_hclustreps(hclust, tip_names, n_reps, sequences):            
    hclust_groups = cluster.hierarchy.cut_tree(hclust, n_clusters = n_reps)
    hclust_groups = [el[0] for el in hclust_groups]
    tips = pd.DataFrame({"tip": tip_names, "hclust_group": hclust_groups})
    groups = []
    reps = []
    for hclust_group, tips_sub in tips.groupby("hclust_group"):
        groups.append(hclust_group)
        seqs = [seq for seq in sequences if seq.id in tips_sub.tip.tolist()]
        rep = select_representative(seqs)
        reps.append(rep.id)
    groups = pd.DataFrame({"hclust_group": groups, "hclustrep": reps})
    tips = tips.merge(groups, on = "hclust_group")
    tips = tips.drop(["hclust_group"], axis = 1)
    return(tips)

## family splitting functions

def split_family_H(pan, sequences, threads, dio_tmp): 
    """Splits a family in two subfamilies.
    
    See split_family_recursive_H.
    """
    write_fasta(sequences, f"{dio_tmp}/seqs.fasta")
    run_mafft(f"{dio_tmp}/seqs.fasta", f"{dio_tmp}/seqs.aln", threads)
    with open(f"{dio_tmp}/seqs.aln", "r") as fin:
        aln = AlignIO.read(fin, "fasta")
    n = len(aln)
    im = identity_matrix(aln)
    dm = []
    for r in range(n - 1):
        for c in range(r + 1, n):
            dm.append(1 - im[r, c])
    link = cluster.hierarchy.linkage(dm, method = "average")
    clusters = cluster.hierarchy.cut_tree(link, n_clusters = 2)
    genes_subfam1 = [seq.id for seq, cl in zip(aln, clusters) if cl[0] == 0]
    genes_subfam2 = [seq.id for seq, cl in zip(aln, clusters) if cl[0] == 1]
    pan1 = pan.loc[genes_subfam1].copy()
    pan2 = pan.loc[genes_subfam2].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    return([pan1, pan2])

def split_family_T(pan, sequences, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_T.
    """
    write_fasta(sequences, f"{dio_tmp}/seqs.fasta")
    run_mafft(f"{dio_tmp}/seqs.fasta", f"{dio_tmp}/seqs.aln", threads)
    run_iqtree(f"{dio_tmp}/seqs.aln", f"{dio_tmp}/tree", threads, 
        ["-m", "LG"])
    tree = Tree(f"{dio_tmp}/tree/tree.treefile")
    midoutgr = tree.get_midpoint_outgroup()
    genes_subfam1 = midoutgr.get_leaf_names()
    midoutgr.detach()
    genes_subfam2 = tree.get_leaf_names()
    pan1 = pan.loc[genes_subfam1].copy()
    pan2 = pan.loc[genes_subfam2].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    return([pan1, pan2])
    
def split_family_H_F(pan, hclust, genes):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_H_F.
    """
    hclust1 = hclust.get_left()
    hclust2 = hclust.get_right()
    genes1 = [genes[i] for i in hclust1.pre_order(lambda x: x.id)]
    genes2 = [genes[i] for i in hclust2.pre_order(lambda x: x.id)]
    pan1 = pan.loc[genes1].copy()
    pan2 = pan.loc[genes2].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    return([pan1, pan2, hclust1, hclust2])

def split_family_LHT_F(pan, strategy, sequences, thresholds, idmat, 
    idmat_names, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_LHT_F.
    """
    
    # some extra parameters
    max_seqs_hclust = 50
    max_seqs_tree = 100
    
    # determine initial step requests
    linclust_requested = 'L' in strategy
    hclust_requested = 'H' in strategy
    tree_requested = 'T' in strategy
    
    # turn off linclust and/or hclust if very few sequences
    n_seqs = len(pan.index)
    if hclust_requested and n_seqs <= max_seqs_hclust:
        linclust_requested = False
    if tree_requested and n_seqs <= max_seqs_tree:
        linclust_requested = False 
        hclust_requested = False
        
    # set default for linrepns_updated
    linrepns_updated = False
    
    # perform linclust steps if requested
    if linclust_requested:
      
        # initialize linreps to constants if not yet present
        if not "linrepc" in pan:
            pan["linrepc"] = "c" # linreps for current threshold
            pan["linrepn"] = "c" # linreps for next threshold
      
        linrepcs = pan["linrepc"].unique().tolist()
        n_linrepcs = len(linrepcs)
      
        # update linreps if necessary and possible
        if n_linrepcs == 1 and thresholds:
            
            # perform linclust until n_linrepcs > 1 or no thresholds left
            pan, thresholds = update_linreps(pan, thresholds, sequences, 
                dio_tmp, threads)
            
            # indicate that linrepns were updated
            linrepns_updated = True
      
            linrepcs = pan["linrepc"].unique().tolist()
            n_linrepcs = len(linrepcs)
            
        # turn off hclust if trees requested and very few linrepns
        if tree_requested and n_linrepns <= max_seqs_tree:
            hclust_requested = False
            
    # perform hclust steps if requested
    if hclust_requested:
    
        # check if idmat needs an update
        update_idmat = False
        if linrepns_updated:
            update_idmat = True
            if n_linrepcs > 1:
                linrepns = pan["linrepn"].unique().tolist()
                idmat_names = linrepns
                idmat_seqs = [seq for seq in sequences if seq.id in linrepns]
            else:
                idmat_names = [seq.id for seq in sequences]
                idmat_seqs = sequences
        # elif idmat is None:
        elif not linclust_requested:
            if idmat is None or len(idmat_names) < len(pan.index):
                update_idmat = True
                idmat_names = [seq.id for seq in sequences]
                idmat_seqs = sequences
            
        # update idmat if necessary
        if update_idmat:
            logging.info(f"aligning {len(idmat_names)} sequences")
            write_fasta(idmat_seqs, f"{dio_tmp}/idmat_seqs.fasta")
            run_mafft(f"{dio_tmp}/idmat_seqs.fasta", 
                f"{dio_tmp}/idmat_seqs.aln", threads)
            aln = read_fasta(f"{dio_tmp}/idmat_seqs.aln")
            idmat = identity_matrix(aln)
            
        # perform hclust
        distmat = distmat_from_idmat(idmat)
        hclust = cluster.hierarchy.linkage(distmat, method = "average")
        
        # select hclustreps if necessary 
        if tree_requested:
            tips = select_hclustreps(hclust, idmat_names, 
                max_seqs_tree, sequences)
            pan = pan.drop(["hclustrep"], axis = 1, errors = "ignore")
            if linclust_requested:
                tips = tips.rename(columns = {"tip": "linrepn"})
                pan = pan.merge(tips, on = "linrepn")
            else:
                tips = tips.rename(columns = {"tip": "gene"})
                pan = pan.reset_index().merge(tips, on = "gene")
                pan = pan.set_index("gene")
            
    # infer tree if requested 
    if tree_requested:
      
        if hclust_requested:
            seqs = [seq for seq in sequences if seq.id in 
                pan.hclustrep.unique().tolist()]
            pan.rep = pan.hclustrep
            
        elif linclust_requested:
            seqs = [seq for seq in sequences if seq.id in 
                pan.linrepn.unique().tolist()]
            pan.rep = pan.linrepn
        
        # remark: this can only happen in T-F, which should re-use the aln
        else:
            seqs = sequences
            pan.rep = pan.gene
      
        write_fasta(seqs, f"{dio_tmp}/seqs.fasta")
        run_mafft(f"{dio_tmp}/seqs.fasta", f"{dio_tmp}/seqs.aln", threads)
        run_iqtree(f"{dio_tmp}/seqs.aln", f"{dio_tmp}/tree", threads, 
            ["-m", "LG"])
        tree = Tree(f"{dio_tmp}/tree/tree.treefile")
        
    # split pan based on hclust
    if not tree_requested:
      
        clusters = cluster.hierarchy.cut_tree(hclust, n_clusters = 2)
        tips_subfam1 = [tip for tip, cl in zip(idmat_names, clusters) 
            if cl[0] == 0]
        tips_subfam2 = [tip for tip, cl in zip(idmat_names, clusters) 
            if cl[0] == 1]
            
        if linclust_requested:
            pan1 = pan[pan["linrepn"].isin(tips_subfam1)].copy()
            pan2 = pan[pan["linrepn"].isin(tips_subfam2)].copy()
        else:
            pan1 = pan.loc[tips_subfam1].copy()
            pan2 = pan.loc[tips_subfam2].copy()
        
        family = pan.orthogroup.tolist()[0]    
        pan1.loc[:, "orthogroup"] = family + "_1"
        pan2.loc[:, "orthogroup"] = family + "_2"
        
    # split pan based on tree
    else:
    
        # split pan based on midpoint root
        midoutgr = tree.get_midpoint_outgroup()
        pan1, pan2 = split_pan(pan, tree, midoutgr)
        
        # # split pan based on minimal copy number root
        # cnoutgr = correct_root(midoutgr, tree, pan)
        # pan_cn1, pan_cn2 = split_pan(pan, tree, cnoutgr)
    
    return([pan1, pan2, thresholds, idmat, idmat_names])

def split_family_P(pan, sequences, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_P.
    """
  
    report = False
        
    # align sequences
    write_fasta(sequences, f"{dio_tmp}/seqs.fasta")
    run_mafft(f"{dio_tmp}/seqs.fasta", f"{dio_tmp}/seqs.aln", threads)

    # create sequence database
    for dir in ["sequenceDB", "tmp", "resultDB", "logs"]:
        makedirs_smart(f"{dio_tmp}/{dir}")
    run_mmseqs(["createdb", f"{dio_tmp}/seqs.fasta", 
        f"{dio_tmp}/sequenceDB/db"], f"{dio_tmp}/logs/createdb.log")
    
    # initialize subfamilies
    with open(f"{dio_tmp}/seqs.aln", "r") as fin:
        aln = AlignIO.read(fin, "fasta")
    genes = pd.DataFrame({"gene": [seq.id for seq in aln],
        "profile": ["profile2"] * len(aln)})
    repr = select_representative(aln)
    if report: print(repr.id)
    genes.loc[genes.gene == repr.id, "profile"] = "profile1"
    genes_profile1 = genes.loc[genes["profile"] == "profile1", "gene"].tolist()
    
    # update profiles
    for i in range(5):
        if report: print(f"iteration {i}")
        for dir in ["msaDB", "profileDB", "searchDB", "tmp"]:
            makedirs_smart(f"{dio_tmp}/{dir}")
        open(f"{dio_tmp}/profiles.sto", "w").close()
        for profile, profilegenes in genes.groupby("profile"):
            records = [copy(rec) for rec in aln if rec.id in \
                profilegenes.gene.tolist()]
            # mmseqs takes id of first sequence in alignment as profile name
            records[0].id = profile
            aln_sub = Align.MultipleSeqAlignment(records)
            with open(f"{dio_tmp}/profiles.sto", "a") as hout:
                AlignIO.write(aln_sub, hout, "stockholm")
        run_mmseqs(["convertmsa", f"{dio_tmp}/profiles.sto", 
            f"{dio_tmp}/msaDB/db"], f"{dio_tmp}/logs/convertmsa.log")
        run_mmseqs(["msa2profile", f"{dio_tmp}/msaDB/db",
            f"{dio_tmp}/profileDB/db", "--match-mode", "1", 
            "--filter-msa", "1", "--diff", "1000", "--qsc", "-50", 
            "--match-ratio", "0.5"], f"{dio_tmp}/logs/msa2profile.log")
        run_mmseqs(["search", f"{dio_tmp}/profileDB/db", 
            f"{dio_tmp}/sequenceDB/db", f"{dio_tmp}/searchDB/db", 
            f"{dio_tmp}/tmp", "--max-seqs", "1000000", "-s", "7.5"], 
            f"{dio_tmp}/logs/search.log")
        run_mmseqs(["convertalis", f"{dio_tmp}/profileDB/db", 
            f"{dio_tmp}/sequenceDB/db", f"{dio_tmp}/searchDB/db", 
            f"{dio_tmp}/results.tsv"], f"{dio_tmp}/logs/convertalis.log")
        hits = read_mmseqs_table(f"{dio_tmp}/results.tsv")
        hits = hits.sort_values("pident", ascending = False)
        hits = hits.reset_index(drop = True)
        hits = hits.drop_duplicates(["target"])
        genes = hits.iloc[:, [0, 1]].copy()
        genes = genes.rename(columns = {"query": "profile", "target": "gene"})
        genes_profile1_new = \
            genes.loc[genes["profile"] == "profile1", "gene"].tolist()
        n_diff = len(set(genes_profile1) ^ set(genes_profile1_new))
        if report: print(f"difference in profile1: {n_diff}")
        genes_profile1 = genes_profile1_new
        if report: print(f"number of genes in profile1: {len(genes_profile1)}")
        if n_diff == 0: break
        if len(genes_profile1) in [len(genes.index), 0]: break

    genes_profile2 = [seq.id for seq in aln if not seq.id in genes_profile1]
    pan1 = pan.loc[genes_profile1].copy()
    pan2 = pan.loc[genes_profile2].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    return([pan1, pan2])
    
## recursive family splitting functions

def split_family_recursive_H(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the H strategy. 
    
    Args:
        pan (DataFrame): A gene table for a single gene family, with the
            columns gene, genome and orthogroup. 
        sequences (list): A list with one SeqRecord object per row in pan.
        threads (int): The number of threads to use. 
        dio_tmp (str): Folder to store temporary files. 
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    family = pan.orthogroup.tolist()[0]
    
    if not split_possible(pan.genome.tolist()):
        logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)
    
    pan1, pan2 = split_family_H(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_H(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_H(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_T(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the T strategy. 
    
    Args:
        pan (DataFrame): A gene table for a single gene family, with the
            columns gene, genome and orthogroup. 
        sequences (list): A list with one SeqRecord object per row in pan.
        threads (int): The number of threads to use. 
        dio_tmp (str): Folder to store temporary files. 
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    family = pan.orthogroup.tolist()[0]
    
    if not split_possible(pan.genome.tolist()):
        logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)
    
    pan1, pan2 = split_family_T(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_T(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_T(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_H_F(pan, hclust, genes):
    """Splits up a gene family using the H-F strategy. 
    
    See [https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.\
        hierarchy.to_tree.html#scipy.cluster.hierarchy.to_tree]. 
    
    Args:
        pan (DataFrame): A gene table for a single gene family, with the
            columns gene, genome and orthogroup. 
        hclust: An hclust object from scipy.hierarchy.cluster in the form of a 
            tree. 
        genes (list): Tip names of the hclust object; tip ids of hclust 
            correspond to indices in genes. Unlike the hclust object, the 
            genes aren't split in iterations of the recursion. Otherwise, the
            tip ids would match the wrong indices. 
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    family = pan.orthogroup.tolist()[0]
    
    if not split_possible(pan.genome.tolist()):
        logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)
    
    pan1, pan2, hclust1, hclust2 = split_family_H_F(pan, hclust, genes)
    split = assess_split(pan1, pan2, family)
    
    if split:

        pan1 = split_family_recursive_H_F(pan1, hclust1, genes)
        pan2 = split_family_recursive_H_F(pan2, hclust2, genes)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_LHT_F(pan, strategy, sequences, thresholds, idmat, 
    idmat_names, threads, dio_tmp):
    """Splits up a gene family using the LHT_F strategy, or subsets of it. 
    
    This splitting function implements three splitting strategy components:
        - linclust (L): the linear-time clustering algorithm of mmseqs2
        - hclust (H): average-linkage hieararchical clustering
        - tree inference (T): phylogeny inference by IQTREE
    
    Linclust can quickly cluster sequences giving a fixed sequence identity
    threshold. This can be used to select a representative subset of the
    sequences ("linreps"). It is super scalable. Hclust can be used to propose
    a bipartition of the sequences (by extracting the two top-level clusters)
    or to make a representative subset of a given number of sequences
    ("hclustreps"). It is relatively scalable. Tree inference can also be used
    to propose a bipartition. It is not very scalable.
    
    A splitting strategy needs the H and/or T component, but otherwise all
    combinations are theoretically possible: H, T, LH, LT, HT, LHT. The
    function makes sure that information from the parent gene family is
    re-used in the splitting process if available (linclusters, distmat;
    re-use of a tree is not yet implemented). Because of this, we add "-F"
    (for "fast") to the strategy names, as in "LHT-F".
    
    Args:
        pan (DataFrame): A table with the columns gene, genome and orthogroup, 
            indexed on the gene column.
        strategy (str): The name of the strategy to use. 
            ["H-F", "T-F", "LH-F", "LT-F", "HT-F", "LHT-F"]
        sequences (list): A list with one SeqRecord object per row in pan.
        thresholds (list): A list with linclust thresholds that have not yet 
            been used. 
        idmat (array): A matrix with percentage identity values for all repris.
        idmat_names (list): A list with the names of the genes of idmat, in 
            the same order as the rows/columns of idmat. 
        threads (int): The number of threads to use. 
        dio_tmp (str): Folder to store temporary files. 
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    max_seqs_tree = 100
    
    family = pan.orthogroup.tolist()[0]

    if not split_possible(pan.genome.tolist()):
        # logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)

    pan1, pan2, thresholds, idmat, idmat_names = split_family_LHT_F(pan, 
        strategy, sequences, thresholds, idmat, idmat_names, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        # idmat is not necessary if few sequences and tree requested
        if 'T' in strategy and len(pan.index) <= max_seqs_tree:
            idmat1 = None
            idmat2 = None
            idmat_names1 = None
            idmat_names2 = None
        else:
            idmat1, idmat_names1 = subset_idmat(idmat, idmat_names, genes1)
            idmat2, idmat_names2 = subset_idmat(idmat, idmat_names, genes2)
        pan1 = split_family_recursive_LHT_F(pan1, strategy, sequences1, 
            thresholds, idmat1, idmat_names1, threads, dio_tmp)
        pan2 = split_family_recursive_LHT_F(pan2, strategy, sequences2, 
            thresholds, idmat2, idmat_names2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_P(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the P strategy. 
    
    Args:
        pan (DataFrame): A gene table for a single gene family, with the
            columns gene, genome and orthogroup. 
        sequences (list): A list with one SeqRecord object per row in pan.
        threads (int): The number of threads to use. 
        dio_tmp (str): Folder to store temporary files. 
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    family = pan.orthogroup.tolist()[0]
    
    if not split_possible(pan.genome.tolist()):
        # logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)
    
    pan1, pan2 = split_family_P(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_P(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_P(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)
    
## top-level functions

def split_superfamily(pan, strategy, din_fastas, threads, dio_tmp):
    """Splits a gene superfamily in a set of gene families.
    
    Args:
        pan (Data Frame): Table with columns gene, genome and orthogroup 
            for the genes or a single gene family. 
        strategy (str): Splitting strategy. [H, T, H-F, LH-F, HT-F, LHT-F, P]
        din_fastas (str): Input folder with fasta file of gene family. 
        threads (int): Number of threads to use. 
        dio_tmp (str): Folder to store temporary files.
        
    Returns:
        A pandas Data Frame with the pangenome, where the orthogroup column has
            been updated to reflect the gene families. 
    """
    
    superfam = pan.orthogroup.tolist()[0]
    pan = pan.set_index("gene")
    dio_tmp = f"{dio_tmp}/{superfam}"
    makedirs_smart(dio_tmp)
    
    if strategy == "H":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_H(pan, sequences, threads, dio_tmp)
    elif strategy == "T":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_T(pan, sequences, threads, dio_tmp)
    elif strategy == "H-F":
        run_mafft(f"{din_fastas}/{superfam}.fasta", f"{dio_tmp}/seqs.aln",
            threads)
        aln = read_fasta(f"{dio_tmp}/seqs.aln")
        genes = [seq.id for seq in aln]
        n = len(aln)
        im = identity_matrix(aln)
        dm = distmat_from_idmat(im)
        hclust = cluster.hierarchy.linkage(dm, method = "average")
        hclust = cluster.hierarchy.to_tree(hclust)
        pan = split_family_recursive_H_F(pan, hclust, genes)
    elif strategy in ["LH-F", "HT-F", "LHT-F"]:
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        # thresholds = [0.50, 0.68, 0.84, 0.92, 0.96, 0.98, 0.99]
        thresholds = generate_thresholds()
        pan = split_family_recursive_LHT_F(pan, strategy, sequences, 
            thresholds = thresholds, idmat = None, idmat_names = None, 
            threads = threads, dio_tmp = dio_tmp)
        pan = pan[["genome", "orthogroup"]] # remark: "gene" is the index
    elif strategy == "P":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_P(pan, sequences, threads, dio_tmp)
    else:
        logging.error("pangenome strategy unknown")
            
    pan = pan.reset_index()
    shutil.rmtree(dio_tmp)
    return(pan)

def infer_superfamilies(faafins, dout, threads):
    """Infers superfamilies of a set of faa files. 
    
    Remarks: the pangenome with the superfamilies is written out to 
    f"{dout}/pangenome.tsv".
    
    Args:
        faafins (list): Paths to .faa fasta files, one per genome. 
        dout (str): Output folder for the pangenome with the superfamilies.
        threads (int): Number of threads to use.
    """
    
    threads = str(threads)
  
    # create output subfolders
    logging.info("creating subfolders for superfamily inference steps")
    for folder in ["logs", "sequenceDB", "preclusterDB", "profileDB", 
        "alignmentDB", "clusterDB", "tmp"]:
        os.makedirs(f"{dout}/{folder}", exist_ok = True)
        
    # create mmseqs sequence database
    logging.info("creating mmseqs sequence database")
    run_mmseqs(["createdb"] + faafins + [f"{dout}/sequenceDB/db"], 
        f"{dout}/logs/createdb.log")
    
    # create preclusters with cluster (could also be linclust)
    logging.info("creating preclusters")
    run_mmseqs(["cluster", f"{dout}/sequenceDB/db", f"{dout}/preclusterDB/db",
        f"{dout}/tmp", "--max-seqs", "1000000", "-c", "0.5", "--cov-mode", "0",
        "-e", "inf", "--min-seq-id", "0.1", "--cluster-mode", "0"], 
        f"{dout}/logs/cluster.log", 
        skip_if_exists = f"{dout}/preclusterDB/db.index", threads = threads)

    # cluster the preclusters into the final clusters
    logging.info("clustering the preclusters")
    run_mmseqs(["result2profile", f"{dout}/sequenceDB/db",
        f"{dout}/sequenceDB/db", f"{dout}/preclusterDB/db",
        f"{dout}/profileDB/db"], f"{dout}/logs/result2profile.log",
        skip_if_exists = f"{dout}/profileDB/db.index", threads = threads)
    run_mmseqs(["search", f"{dout}/profileDB/db",
        f"{dout}/profileDB/db_consensus", f"{dout}/alignmentDB/db",
        f"{dout}/tmp", "-c", "0.5", "--cov-mode", "1"], 
        f"{dout}/logs/search.log",
        skip_if_exists = f"{dout}/alignmentDB/db.index", threads = threads)
    run_mmseqs(["clust", f"{dout}/profileDB/db", f"{dout}/alignmentDB/db",
        f"{dout}/clusterDB/db", "--cluster-mode", "2"],
        f"{dout}/logs/clust.log",
        skip_if_exists = f"{dout}/clusterDB/db.index", threads = threads)

    # create the tsv files with the preclusters and clusters
    logging.info("compiling pangenome file with superfamilies")
    run_mmseqs(["createtsv", f"{dout}/sequenceDB/db", f"{dout}/sequenceDB/db",
        f"{dout}/preclusterDB/db", f"{dout}/preclusters.tsv"], 
        f"{dout}/logs/createtsv_preclusters.log")
    run_mmseqs(["createtsv", f"{dout}/sequenceDB/db", f"{dout}/sequenceDB/db",
        f"{dout}/clusterDB/db", f"{dout}/clusters.tsv"],
        f"{dout}/logs/createtsv_clusters.log")
    preclustertable = pd.read_csv(f"{dout}/preclusters.tsv", sep = "\t", 
        names = ["precluster", "gene"])
    clustertable = pd.read_csv(f"{dout}/clusters.tsv", sep = "\t",
        names = ["cluster", "precluster"])
    genes_ogs = pd.merge(preclustertable, clustertable, on = "precluster")
    genes_ogs = genes_ogs.rename(columns = {"cluster": "orthogroup"})
    genes_genomes = extract_genes(faafins)
    genes = pd.merge(genes_genomes, genes_ogs, on = "gene")
    genes = genes.drop(["precluster"], axis = 1)
    write_tsv(genes, f"{dout}/pangenome.tsv")

def infer_pangenome(faafins, splitstrategy, dout, threads):
    """Infers the pangenome of a set of genomes and writes it to disk. 
    
    Args:
        faafins (list): Paths to .faa fasta files, one per genome. 
        splitstrategy (str): Splitting strategy.
        dout (str): Output folder for the pangenome.
        threads (int): Number of threads to use.
    """
  
    # threads per process
    tpp = 1
  
    logging.info("STAGE 1: creation of superfamilies")
    os.makedirs(f"{dout}/superfamilies", exist_ok = True)
    infer_superfamilies(faafins, f"{dout}/superfamilies", threads)
      
    logging.info("STAGE 2: splitting of superfamilies")

    logging.info("gathering sequences of superfamilies")
    os.makedirs(f"{dout}/superfamilies/superfamilies", exist_ok = True)
    pangenome = read_genes(f"{dout}/superfamilies/pangenome.tsv")
    dio_fastas = f"{dout}/superfamilies/fastas"
    gather_orthogroup_sequences(pangenome, faafins, dio_fastas)

    logging.info("splitting superfamilies")
    os.makedirs(f"{dout}/tmp", exist_ok = True)
    pangenome = pangenome.groupby("orthogroup")
    orthogroups = pangenome.aggregate({"genome": split_possible})
    splitable = orthogroups[orthogroups.genome].index.tolist()
    pangenome_splitable = [pan for name, pan in pangenome if name in splitable]
    pangenome_splitable.sort(key = lambda pan: len(pan.index), reverse = True)
    n = len(pangenome_splitable)
    logging.info(f"splitting {n} splitable superfamilies")
    with ProcessPoolExecutor(max_workers = threads // tpp) as executor:
        pangenome_splitable = executor.map(split_superfamily, 
            pangenome_splitable, [splitstrategy] * n, [dio_fastas] * n, 
            [tpp] * n, [f"{dout}/tmp"] * n)
    pangenome = pd.concat(list(pangenome_splitable) +
        [pan for name, pan in pangenome if not name in splitable])
    write_tsv(pangenome, f"{dout}/pangenome.tsv")

def run_pan_nonhier(args):

    pangenomefout = os.path.join(args.outfolder, "pangenome.tsv")
    if os.path.isfile(pangenomefout):
        logging.info("existing pangenome detected - moving on")
        return()
    
    faafins = read_lines(args.faapaths)
    
    if args.method in ["O-B", "O-D"]:

        logging.info("creating orthofinder directory")
        dir_orthofinder = os.path.join(args.outfolder, "orthofinder")
        os.makedirs(dir_orthofinder, exist_ok = True)

        logging.info("running orthofinder")
        logfile = os.path.join(args.outfolder, "orthofinder.log")
        
        if args.method == "O-B":
            engine = "blast"
        else:
            engine = "diamond"
        run_orthofinder(faafins, dir_orthofinder, logfile, args.threads, 
            engine)
    
        logging.info("creating tidy pangenome file")
        pangenome = read_pangenome_orthofinder(dir_orthofinder)
        write_tsv(pangenome, pangenomefout)
        
    else:
      
        logging.info(f"constructing pangenome with {args.method} strategy")
        infer_pangenome(faafins, args.method, args.outfolder, args.threads)

def run_pan_hier(args):

    logging.info("processing species file")
    genomes_species = read_species(args.species)
    genomes_species.species = [species.replace(" ", "_") for species in
        genomes_species.species]
    speciesdict = {}
    for ix, row in genomes_species.iterrows():
        speciesdict.setdefault(row.species, []).append(row.genome)

    logging.info("processing faapaths")
    faafins = read_lines(args.faapaths)
    genomedict = {}
    for faafin in faafins:
        genome = filename_from_path(faafin)
        genomedict[genome] = faafin

    logging.info("started building pangenomes of species")
    speciespansdio = os.path.join(args.outfolder, "speciespangenomes")
    os.makedirs(speciespansdio, exist_ok = True)
    reprpaths = []
    speciespanfios = []
    for species, genomes in speciesdict.items():
        logging.info(f"inferring pangenome of {species}")
        faapaths = [genomedict[genome] for genome in genomes]
        dout = os.path.join(speciespansdio, species)
        os.makedirs(dout, exist_ok = True)
        faapathsfio = os.path.join(dout, "faapaths.txt")
        write_lines(faapaths, faapathsfio)
        run_pan_nonhier(Namespace(faapaths = faapathsfio, outfolder = dout,
            threads = args.threads, method = "LHT-F"))
        speciespanfio = os.path.join(dout, "pangenome.tsv")
        speciespanfios.append(speciespanfio)
        reprfio = os.path.join(dout, species + ".faa")
        reprpaths.append(reprfio)
        speciespan = read_genes(speciespanfio)
        collapse_pangenome(speciespan, faapathsfio, reprfio, species,
            os.path.join(dout, "temp"))
    reprpathsfio = os.path.join(args.outfolder, "reprpaths.txt")
    write_lines(reprpaths, reprpathsfio)

    logging.info("started building metapangenome using representatives")
    metapandio = os.path.join(args.outfolder, "metapangenome")
    os.makedirs(metapandio, exist_ok = True)
    run_pan_nonhier(Namespace(faapaths = reprpathsfio,
        outfolder = metapandio, threads = args.threads,
        method = "LHT-F"))

    logging.info("started inflating metapangenome with species pangenomes")
    speciespans = [read_genes(panfin) for panfin in speciespanfios]
    speciespans = pd.concat(speciespans)
    speciespans = speciespans.rename(columns = {"orthogroup": "speciesfam"})
    speciespans = pd.merge(speciespans, genomes_species)
    speciespans.speciesfam = [species + "-" + speciesfam for species,
        speciesfam in zip(speciespans.species, speciespans.speciesfam)]
    metapan = read_genes(os.path.join(metapandio, "pangenome.tsv"))
    metapan = metapan.rename(columns = {"gene": "speciesfam",
        "genome": "species"})
    pangenome = pd.merge(speciespans, metapan)
    pangenome = pangenome[["gene", "genome", "orthogroup"]]
    write_tsv(pangenome, os.path.join(args.outfolder, "pangenome.tsv"))
