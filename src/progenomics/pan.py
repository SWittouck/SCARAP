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
    reprfs_subfam1 = outgr.get_leaf_names()
    parent = outgr.up
    outgr.detach()
    reprfs_subfam2 = tree.get_leaf_names()
    parent.add_child(outgr) # re-add the outgroup
    pan1 = pan[pan["reprf"].isin(reprfs_subfam1)].copy()
    pan2 = pan[pan["reprf"].isin(reprfs_subfam2)].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    return([pan1, pan2])

def lowest_cn_roots(tree, pan):
    """Determine the set of lowest copy-number roots 
    
    Args:
        tree: ete3 tree object where the leaf names correspond to the values of
            the reprf column in pan
        pan (DataFrame): Table with at least the columns reprf and genome
        
    Returns:
        A list with references to the nodes that would be minimal copy number 
            roots
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
    if split:
        logging.info(f"{family}: {n1}/{n2}; "
            f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
            f"decision = split")
    else:
        logging.info(f"{family}: {n1}/{n2}; "
            f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
            f"decision = no split")
    return(split)
    
## strategy-specific splitting helpers

def update_linclusters(sequences, thresholds, dio_tmp, threads):
    """Updates linclusters until there are more than ten, if possible. 
    
    Helper function for the TRE-F(S) splitting strategy.
    
    Args:
        sequences (list): A list of SeqRecords.
        thresholds (list): A list of percentage identity thresholds to try, in
            ascending order.
        dio_tmp (str): A folder to store temporary files.
        threads (int): Number of threads to use.
        
    Returns:
        A pandas dataframe with the columns gene and lincluster.
        An updated list of thresholds where the thresholds that have already
            been tried are removed.
    """
    max_reprf = 10
    thresholds = thresholds.copy()
    write_fasta(sequences, f"{dio_tmp}/seqs.fasta")        
    for dir in ["sequenceDB", "logs", "tmp"]:
        makedirs_smart(f"{dio_tmp}/{dir}")
    run_mmseqs(["createdb", f"{dio_tmp}/seqs.fasta", 
        f"{dio_tmp}/sequenceDB/db"], f"{dio_tmp}/logs/createdb.log", 
        threads = threads)
    n_lin = 0
    while thresholds and n_lin <= max_reprf:
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
            names = ["lincluster", "gene"])                        
        n_lin = genes_clusters["lincluster"].nunique()
    return(genes_clusters, thresholds)
        
def add_reprfs(pan, idmat, repris, sequences):
    """Adds reprfs to a pan table.
    
    Selects final representative sequences (reprfs) by performing hierarchical 
    clustering on intermediate representative sequences (repris). 
    
    Helper function for the TRE-F(S) splitting strategy.
    
    Args:
        pan: A data frame with columns gene and repri.
        idmat: A percentage identity matrix of the repris.
        repris: The rownames (= colnames) of idmat.
        sequences: A list of SeqRecords.
        
    Returns:
        An updated pan table where a column "reprf" is added.
    """
    max_reprf = 10
    pan = pan.drop(["reprf"], axis = 1, errors = "ignore")
    distmat = distmat_from_idmat(idmat)
    hclust = cluster.hierarchy.linkage(distmat, method = "average")
    hclust_groups = cluster.hierarchy.cut_tree(hclust, n_clusters = max_reprf)
    hclust_groups = [el[0] for el in hclust_groups]
    repris_df = pd.DataFrame({"repri": repris, "hclust_group": hclust_groups})
    pan = pan.reset_index().merge(repris_df, on = "repri").set_index("gene")
    hclust_groups = []
    reprfs = []
    for hclust_group, pansub in pan.groupby("hclust_group"):
        hclust_groups.append(hclust_group)
        seqs = [seq for seq in sequences if seq.id in pansub.index.tolist()]
        reprf = select_representative(seqs)
        reprfs.append(reprf.id)
    hclust_groups_df = pd.DataFrame({"hclust_group": hclust_groups, 
        "reprf": reprfs})
    pan = pan.reset_index().merge(hclust_groups_df, on = "hclust_group")
    pan = pan.set_index("gene")
    pan = pan.drop(["hclust_group"], axis = 1)
    return(pan)

## family splitting functions

def split_family_TRE(pan, sequences, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_TRE.
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

def split_family_TRE_F(pan, sequences, idmat, repris, thresholds, threads, 
    dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_TRE_F.
    """
    
    max_reprf = 10
    
    # if ten or fewer genes --> all genes are reprfs
    if len(pan.index) <= max_reprf:
        
        pan["reprf"] = pan.index 
        
    # if more than max_reprf genes --> reprfs are a subset of repris
    else:
        
        # update repris and their idmat if necessary
        if idmat is None or len(repris) <= max_reprf:
            
            n_lin = 0
            
            # remove the previous repri column
            pan = pan.drop("repri", axis = 1, errors = "ignore")
            
            # 1) update linclusters if thresholds left
            if thresholds:
                genes_clusters, thresholds = update_linclusters(sequences, 
                    thresholds, dio_tmp, threads)    
                n_lin = genes_clusters["lincluster"].nunique()
                    
            # 2a) if enough linclusters: select repris from linclusters
            if n_lin > max_reprf:
                repris = genes_clusters["lincluster"].unique().tolist()
                pan = pan.reset_index().merge(genes_clusters, on = "gene")
                pan = pan.set_index("gene")
                pan = pan.rename(columns = {"lincluster": "repri"})
                
            # 2b) if not enough linclusters: repris are all genes
            else:
                repris = pan.index.tolist()
                pan["repri"] = pan.index
                
            # 3) update idmat
            seqs_repris = [seq for seq in sequences if seq.id in repris]
            write_fasta(seqs_repris, f"{dio_tmp}/repris.fasta")
            run_mafft(f"{dio_tmp}/repris.fasta", f"{dio_tmp}/repris.aln", 
                threads)
            aln = read_fasta(f"{dio_tmp}/repris.aln")
            idmat = identity_matrix(aln)

        # select (at max) max_reprf reprfs
        pan = add_reprfs(pan, idmat, repris, sequences)
        
    # build tree from reprfs
    seqs = [seq for seq in sequences if seq.id in pan.reprf.unique().tolist()]
    write_fasta(seqs, f"{dio_tmp}/reprfs.fasta")
    run_mafft(f"{dio_tmp}/reprfs.fasta", f"{dio_tmp}/reprfs.aln", threads)
    run_iqtree(f"{dio_tmp}/reprfs.aln", f"{dio_tmp}/tree", threads, 
        ["-m", "LG"])
    tree = Tree(f"{dio_tmp}/tree/tree.treefile")
    
    # split pan based on midpoint root
    midoutgr = tree.get_midpoint_outgroup()
    pan1, pan2 = split_pan(pan, tree, midoutgr)
    
    # split pan based on minimal copy number root
    roots = lowest_cn_roots(tree, pan)
    if midoutgr in roots:
        cnoutgr = midoutgr
        # logging.info("cn root is midpoint root")
    else:
        cnoutgr = roots[0]
        # logging.info("cn root is not midpoint root")
    pan_cn1, pan_cn2 = split_pan(pan, tree, cnoutgr)
    
    return([pan1, pan2, pan_cn1, pan_cn2, idmat, repris, thresholds])

def split_family_CLU(pan, sequences, threads, dio_tmp): 
    """Splits a family in two subfamilies.
    
    See split_family_recursive_CLU.
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
    
def split_family_CLU_F(pan, hclust, genes):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_CLU_F.
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

def split_family_PRO(pan, sequences, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_PRO.
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

def split_family_recursive_TRE(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the TRE strategy. 
    
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
    
    pan1, pan2 = split_family_TRE(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_TRE(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_TRE(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_TRE_F(pan, sequences, idmat, repris, thresholds, 
    threads, dio_tmp):
    """Splits up a gene family using the TRE-F strategy. 
    
    Two types of representative genes are used in this strategy:
        Intermediate representatives (repris): these are selected from the
            genes using linclust if possible; otherwise they are just all genes
            (no subset).
        Final representatives (reprfs): these are selected from the repris
            through hierarchical clustering. There are 10 or fewer of them.
            They will be used to construct a phylogeny that is then split in
            two. 
    
    If idmat is supplied, repris should also be supplied and pan should have 
        an additional column "repri".
    
    Args:
        pan (DataFrame): A table with the columns gene, genome and orthogroup, 
            indexed on the gene column.
        sequences (list): A list with one SeqRecord object per row in pan.
        idmat (array): A matrix with percentage identity values for all repris.
        repris (list): A list with repris in the same order as the rows and 
            columns of idmat.
        thresholds (list): A list with linclust thresholds that can be tried. 
        threads (int): The number of threads to use. 
        dio_tmp (str): Folder to store temporary files. 
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    max_reprf = 10
    
    family = pan.orthogroup.tolist()[0]

    if not split_possible(pan.genome.tolist()):
        logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)

    pan1, pan2, pan_cn1, pan_cn2, idmat, repris, thresholds = \
        split_family_TRE_F(pan, sequences, idmat, repris, thresholds, threads, 
        dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan_cn1.index.tolist()
        genes2 = pan_cn2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        # idmat and repris are only necessary when there are > max_reprf genes:
        if len(pan.index) > max_reprf:
            idmat1, repris1 = subset_idmat(idmat, repris, genes1)
            idmat2, repris2 = subset_idmat(idmat, repris, genes2)
        else:
            idmat1 = None
            idmat2 = None
            repris1 = None
            repris2 = None
        pan1 = split_family_recursive_TRE_F(pan_cn1, sequences1, idmat1, 
            repris1, thresholds, threads, dio_tmp)
        pan2 = split_family_recursive_TRE_F(pan_cn2, sequences2, idmat2, 
            repris2, thresholds, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_CLU(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the CLU strategy. 
    
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
    
    pan1, pan2 = split_family_CLU(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_CLU(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_CLU(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_CLU_F(pan, hclust, genes):
    """Splits up a gene family using the CLU-F strategy. 
    
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
    
    pan1, pan2, hclust1, hclust2 = split_family_CLU_F(pan, hclust, genes)
    split = assess_split(pan1, pan2, family)
    
    if split:

        pan1 = split_family_recursive_CLU_F(pan1, hclust1, genes)
        pan2 = split_family_recursive_CLU_F(pan2, hclust2, genes)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_PRO(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the PRO strategy. 
    
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
    
    pan1, pan2 = split_family_PRO(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_PRO(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_PRO(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)
    
## top-level functions

def split_superfamily(pan, strategy, din_fastas, threads, dio_tmp):
    """Splits a gene superfamily in a set of gene families.
    
    Args:
        pan (Data Frame): Table with columns gene, genome and orthogroup 
            for the genes or a single gene family. 
        strategy (str): Splitting strategy. [TRE, TRE-F, TRE-FS, CLU, CLU-F, 
            PRO]
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
    
    if strategy == "TRE":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_TRE(pan, sequences, threads, dio_tmp)
    elif strategy in ["TRE-F", "TRE-FS"]:
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        if strategy == "TRE-F":
            thresholds = []
        else:
            thresholds = [0.50, 0.68, 0.84, 0.92, 0.96, 0.98, 0.99]
        pan = split_family_recursive_TRE_F(pan, sequences, idmat = None, 
            repris = None, thresholds = thresholds, threads = threads, 
            dio_tmp = dio_tmp)
        pan = pan[["genome", "orthogroup"]] # remark: "gene" is the index
    elif strategy == "CLU":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_CLU(pan, sequences, threads, dio_tmp)
    elif strategy == "CLU-F":
        run_mafft(f"{din_fastas}/{superfam}.fasta", f"{dio_tmp}/seqs.aln", 
            threads)
        aln = read_fasta(f"{dio_tmp}/seqs.aln")
        genes = [seq.id for seq in aln]
        n = len(aln)
        im = identity_matrix(aln)
        dm = distmat_from_idmat(im)
        hclust = cluster.hierarchy.linkage(dm, method = "average")
        hclust = cluster.hierarchy.to_tree(hclust)
        pan = split_family_recursive_CLU_F(pan, hclust, genes)
    elif strategy == "PRO":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_PRO(pan, sequences, threads, dio_tmp)
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
        splitstrategy (str): Splitting strategy. ["trees", "profiles"]
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
    
    if args.method in ["ORT-B", "ORT-D"]:

        logging.info("creating orthofinder directory")
        dir_orthofinder = os.path.join(args.outfolder, "orthofinder")
        os.makedirs(dir_orthofinder, exist_ok = True)

        logging.info("running orthofinder")
        logfile = os.path.join(args.outfolder, "orthofinder.log")
        
        if args.method == "ORT-B":
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
            threads = args.threads, method = "TRE-F"))
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
        method = "TRE-F"))

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
