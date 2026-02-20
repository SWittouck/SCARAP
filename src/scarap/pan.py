import collections
import logging
import numpy as np
import os
import shutil
import sys

from Bio import AlignIO, Align
from copy import copy
from ete4 import Tree
from concurrent.futures import ProcessPoolExecutor
from scipy import cluster

from scarap.utils import *
from scarap.readerswriters import *
from scarap.computers import *
from scarap.callers import *
from scarap.helpers import archive_faafins, MAX_FAA_ARG_LEN

## helpers - ficlin module (F)
                
def update_seedmatrix(seedmatrix, sequences, dout_tmp, threads):
    """Updates the identity values in a seedmatrix. 
    
    A seed matrix is a matrix where the rows represent genes, while the
    columns represent a subset of these genes ("seeds"). The cells contain
    sequence identity values of the genes to the seeds. This function will
    identify seeds that are no longer present in the set of genes (= seeds
    that don't have an identity value of one to at least one of the genes),
    and replace them with new seeds.
    
    Args:
        seedmatrix (np.array): An array of identity values of genes to seeds. 
        sequences (list): A list of SeqRecords in the same order as the rows 
            of the seedmatrix. 
        dout_tmp (str): The path to a folder to store temporary files. 
        threads (int): The number of threads to use. 
        
    Returns:
        An updated seedmatrix. 
    """
  
    # construct target and prefilter dbs
    makedirs_smart(f"{dout_tmp}")
    for dir in ["sequenceDB", "prefDB", "logs", "tmp"]:
        makedirs_smart(f"{dout_tmp}/{dir}")
    write_fasta(sequences, f"{dout_tmp}/seqs.fasta")
    run_mmseqs(["createdb", f"{dout_tmp}/seqs.fasta", 
        f"{dout_tmp}/sequenceDB/db"], f"{dout_tmp}/logs/createdb.log", 
        threads = threads)
    create_prefdb("../sequenceDB/db", f"{dout_tmp}/prefDB/db")
    
    # identify seeds (= columns) to replace
    replace = np.apply_along_axis(lambda l: all(l != 1), 0, seedmatrix)
    seedstoreplace = [i for i, r in enumerate(replace) if r]
    # logging.info(f"adding {len(seedstoreplace)} seeds")
    
    # update seed matrix with identity values 
    genes = [seq.id for seq in sequences]
    lengths = [len(seq) for seq in sequences]
    for c in seedstoreplace:
      
        for dir in ["seedDB", "alignmentDB"]:
            makedirs_smart(f"{dout_tmp}/{dir}")
        
        # select longest next seed
        max_ids = np.amax(seedmatrix, 1) # max of each row = max along columns
        min_max_ids = np.amin(max_ids)
        cand_length = 0
        for i, s in enumerate(sequences):
            if max_ids[i] == min_max_ids and lengths[i] > cand_length:
                cand_length = lengths[i]
                seed = s
                seed_ix = i
        
        # create seed db
        write_fasta([seed], f"{dout_tmp}/seed.fasta")
        run_mmseqs(["createdb", f"{dout_tmp}/seed.fasta", 
            f"{dout_tmp}/seedDB/db"], f"{dout_tmp}/logs/createseeddb.log", 
            threads = threads)
            
        # run mmseqs align
        update_prefdb(f"{dout_tmp}/seedDB/db", f"{dout_tmp}/sequenceDB/db", 
            f"{dout_tmp}/prefDB/db")
        run_mmseqs(["align", f"{dout_tmp}/seedDB/db",
            f"{dout_tmp}/sequenceDB/db", f"{dout_tmp}/prefDB/db", 
            f"{dout_tmp}/alignmentDB/db", "--alignment-mode", "1"], 
            f"{dout_tmp}/logs/search.log", threads = threads)
        run_mmseqs(["createtsv", f"{dout_tmp}/seedDB/db", 
            f"{dout_tmp}/sequenceDB/db", f"{dout_tmp}/alignmentDB/db", 
            f"{dout_tmp}/hits.tsv", "--full-header"], 
            f"{dout_tmp}/logs/createtsv.log")
        hits_new = pd.read_csv(f"{dout_tmp}/hits.tsv", sep = "\t", 
            usecols = [1, 3], names = ["gene", "identity"])
        hits_new["gene"] = hits_new["gene"].apply(lambda x: x.split(" ")[0])
        
        # put the identity values in the seed matrix
        for index, row in hits_new.iterrows():
            gene_ix = genes.index(row["gene"])
            seedmatrix[gene_ix, c] = row["identity"]
            
        # set identity of seed to itself to one
        # (mmseqs align estimates the identity values from the scores)
        seedmatrix[seed_ix, c] = 1
    
    # give warning if some sequences don't align to their cluster seed
    ids_to_seed = np.amax(seedmatrix, 1)
    seeds_not_aligned = np.count_nonzero(ids_to_seed == 0)
    if seeds_not_aligned:
        logging.warning(f"ficlin: {seeds_not_aligned} sequences do not align to any "
            "seed")
    
    # remove temporary output folder
    shutil.rmtree(dout_tmp)
        
    return(seedmatrix)

def run_ficlin(sequences, n_clusters, dout_tmp, threads):
    """Partitions sequences in a fixed number of clusters. 
    
    Clusters a set of sequences into a fixed number of clusters in linear time
    and memory.
    
    Remark: this function hasn't been tested yet with multiple threads; I'm
    not sure if it can use them efficiently.
    
    Args:
        sequences (list): A list of SeqRecords to cluster.
        n_clusters (int): The number of clusters requested.
        dout_tmp (str): The path to a folder to store temporary files. 
        threads (int): The number of threads to use. 
        
    Returns:
        A list containing the cluster number for each sequence in sequences.
    """
    
    # initialize a seed matrix 
    seedmatrix = np.zeros([len(sequences), n_clusters])
    
    # update the seedmatrix
    seedmatrix = update_seedmatrix(seedmatrix, sequences, dout_tmp, threads)
    
    # determine the cluster of each gene
    clusters = np.argmax(seedmatrix, 1)
    
    return(clusters)
    
## helpers - hierarchical clustering module (H)
    
def hclusts(distmat, n_clusters):
  
    hclust = cluster.hierarchy.linkage(distmat, method = "average")
    clusters = cluster.hierarchy.cut_tree(hclust, n_clusters = n_clusters)
    clusters = [el[0] for el in clusters]
    
    return(clusters)

## helpers - tree inference module (T)
    
def split_pan(pan, tree):
    """Split pan into pan1 and pan2 based on a tree and an outgroup node. 
    
    Args:
        pan (DataFrame): A gene table with at least the columns reprf and 
            orthogroup.
        tree: An ete4 tree (= the root node of a tree)
        
    Returns:
        [pan1, pan2, tree1, tree2]
    """
    
    # midpoint root the tree
    midoutgr = tree.get_midpoint_outgroup()
    if tree.root.dist:
        tree.root.dist = None
    
    if midoutgr != tree:
        tree.set_outgroup(midoutgr) 
    # split tree at root
    tree1 = tree.children[0].copy()
    tree2 = tree.children[1].copy()
    
    # split pan
    reps_subfam1 = tree1.leaf_names()
    reps_subfam2 = tree2.leaf_names()
    pan1 = pan[pan["rep"].isin(reps_subfam1)].copy()
    pan2 = pan[pan["rep"].isin(reps_subfam2)].copy()
    
    # set subfamily names
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"

    return([pan1, pan2, tree1, tree2])

def lowest_cn_roots(tree, pan):
    """Determine the set of lowest copy-number roots.
    
    Args:
        tree: ete4 tree object where the leaf names correspond to the values of
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
    for node in tree.descendants():
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
    for node in tree.descendants():
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
    
## general helpers
    
def select_reps(genes, clusters, sequences):
    """Select representative sequences for a set of clusters.
    
    For each cluster, selects a sequence with median length as its
    representative.
    
    Args: 
        genes (list): A list of gene names.
        clusters (list): A list of cluster ids in the same order as genes.
        sequences (list): A list of SeqRecords containing at least the
            sequences of all genes; the order doesn't matter.
            
    Returns:
        A DataFrame with the columns gene and rep. 
    """
    genes = pd.DataFrame({"gene": genes, "cluster": clusters})
    n_genes = len(genes.index)
    # add sequences to gene table without changing order
    genes_seqs = pd.DataFrame({"gene": [s.id for s in sequences], 
        "sequence": sequences})
    genes = genes.merge(genes_seqs, on = "gene", how = "left")
    if (genes.sequence.isnull().values.any()):
        logging.error("select_reps: not all sequences were given")
    # the following doesn't change the order of the rows
    genes["rep"] = genes.groupby("cluster")["sequence"].\
        transform(lambda x: select_representative(x.tolist()).id)
    return(genes.rep.tolist())
    
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
    genomes = genomes1 + genomes2
    cns = pd.Series(genomes).value_counts().tolist()
    if all([cn > 10 for cn in cns]) and not split:
        logging.warning(f"{family}: all copy-numbers > 10 but no split")
        # logging.info(f"subfam1: {pan1.index.tolist()}")
        # logging.info(f"subfam2: {pan2.index.tolist()}")
    # n = len(set(genomes))
    # logging.info(f"{family}: {n} genomes; copy-number {min(cns)} - "
    #     f"{max(cns)}; pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
    #     f"split = {split}")
    return(split)

## family splitting functions

def split_family_H_nl(pan, sequences, threads, dio_tmp):
    """Splits a family in two subfamilies.

    See split_family_recursive_H_nl.
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

def split_family_T_nl(pan, sequences, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_T_nl.
    """
    write_fasta(sequences, f"{dio_tmp}/seqs.fasta")
    run_mafft(f"{dio_tmp}/seqs.fasta", f"{dio_tmp}/seqs.aln", threads)
    run_iqtree(f"{dio_tmp}/seqs.aln", f"{dio_tmp}/tree", threads, 
        ["-m", "LG+F+G4"])
    with open(f"{dio_tmp}/tree/tree.treefile", "r") as ftree:
        tree = Tree(ftree)
    midoutgr = tree.get_midpoint_outgroup()
    genes_subfam1 = midoutgr.leaf_names()
    midoutgr.detach()
    genes_subfam2 = tree.leaf_names()
    pan1 = pan.loc[genes_subfam1].copy()
    pan2 = pan.loc[genes_subfam2].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    return([pan1, pan2])
    
def split_family_FH(pan, sequences, hclust, ficlin, min_reps, max_reps, 
    max_align, seedmatrix, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_FH.
    """

    # determine necessary steps
    n_seqs = len(pan.index)
    n_reps = 0 if hclust is None else hclust.get_count()
    if ficlin and (n_reps >= min_reps or n_reps == n_seqs):
        update_reps = False
        cluster_reps = False
    elif ficlin and n_seqs > max_align:
        update_reps = True
        cluster_reps = True
    else:
        pan["rep"] = pan.index # all seqs become reps 
        update_reps = False
        cluster_reps = True

    # STEP 1: UPDATE REPRESENTATIVES

    if update_reps:

        seedmatrix = update_seedmatrix(seedmatrix, sequences, 
            f"{dio_tmp}/ficlin", threads)
        linclusters = np.argmax(seedmatrix, 1)
        # return emptiness if all sequences map to the same seed
        if (len(set(linclusters))) == 1: 
            return([None] * 8)
        pan["rep"] = select_reps(pan.index.tolist(), linclusters, sequences)
     
    # STEP 2: CLUSTER REPRESENTATIVES
    
    if cluster_reps:

        repseqs = [s for s in sequences if s.id in pan["rep"].unique()]
        reps = [s.id for s in repseqs]
        # logging.info(f"aligning {len(reps)} sequences")
        write_fasta(repseqs, f"{dio_tmp}/repseqs.fasta")
        run_mafft(f"{dio_tmp}/repseqs.fasta", f"{dio_tmp}/repseqs.aln", 
            threads, ["--amino", "--anysymbol"])
        aln = read_fasta(f"{dio_tmp}/repseqs.aln")
        idmat = identity_matrix(aln)
        distmat = distmat_from_idmat(idmat)
        hclust = cluster.hierarchy.linkage(distmat, method = "average")
        hclust = cluster.hierarchy.to_tree(hclust)
        for leaf in hclust.pre_order(lambda x: x):
            leaf.id = reps[leaf.id]
    
    # STEP 3: SPLIT DATA IN SUBFAMILIES

    # split hclust
    hclust1 = hclust.get_left()
    hclust2 = hclust.get_right()
    
    # split pan
    reps1 = hclust1.pre_order(lambda x: x.id)
    reps2 = hclust2.pre_order(lambda x: x.id)
    pan1 = pan[pan["rep"].isin(reps1)].copy()
    pan2 = pan[pan["rep"].isin(reps2)].copy()
    family = pan.orthogroup.tolist()[0]
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"

    # split sequences
    sequences1 = [s for s in sequences if s.id in pan1.index]
    sequences2 = [s for s in sequences if s.id in pan2.index]

    # split seedmatrix
    if ficlin:
        seedmatrix1 = seedmatrix[pan["rep"].isin(reps1).tolist(), :]
        seedmatrix2 = seedmatrix[pan["rep"].isin(reps2).tolist(), :]
    else: 
        seedmatrix1 = None
        seedmatrix2 = None
    
    return([pan1, pan2, sequences1, sequences2, hclust1, hclust2, seedmatrix1, 
        seedmatrix2])
    
def split_family_FT(pan, sequences, tree, ficlin, min_reps, max_reps, 
    seedmatrix, threads, dio_tmp):
    """Splits a family in two subfamilies.
    
    See split_family_recursive_FT.
    """
    
    update_tree = False
    # if reps not necessary, for the first time: ...
    if not ficlin or len(pan.index) <= max_reps:
        if tree is None or len(tree) != len(pan.index):
            pan["rep"] = pan.index
            update_tree = True
    # if parent reps are too few to re-use: ...
    else:
        n_reps = len(pan["rep"].unique())
        if n_reps < min_reps:
            seedmatrix = update_seedmatrix(seedmatrix, sequences, 
                f"{dio_tmp}/ficlin", threads)
            linclusters = np.argmax(seedmatrix, 1)
            # return emptiness if all sequences map to the same seed
            if (len(set(linclusters))) == 1: 
                return([None] * 8)
            pan["rep"] = select_reps(pan.index.tolist(), linclusters, sequences)
            update_tree = True
            
    # update tree if requested, otherwise use parent tree
    if update_tree:
        repseqs = [s for s in sequences if s.id in pan["rep"].unique()]
        reps = [s.id for s in repseqs]
        # logging.info(f"aligning {len(reps)} sequences")
        write_fasta(repseqs, f"{dio_tmp}/repseqs.fasta")
        run_mafft(f"{dio_tmp}/repseqs.fasta", f"{dio_tmp}/repseqs.aln", 
            threads, ["--amino"])
        run_iqtree(f"{dio_tmp}/repseqs.aln", f"{dio_tmp}/tree", threads, 
            ["-m", "LG"])
        with open(f"{dio_tmp}/tree/tree.treefile", "r") as ftree:
            tree = Tree(ftree)
    # split pan based on midpoint root
    pan1, pan2, tree1, tree2 = split_pan(pan, tree)
    sequences1 = [s for s in sequences if s.id in pan1.index]
    sequences2 = [s for s in sequences if s.id in pan2.index]
    if ficlin:
        seedmatrix1 = seedmatrix[pan.index.isin(pan1.index).tolist(), :]
        seedmatrix2 = seedmatrix[pan.index.isin(pan2.index).tolist(), :]
    else: 
        seedmatrix1 = None
        seedmatrix2 = None
    
    # # split pan based on minimal copy number root
    # cnoutgr = correct_root(midoutgr, tree, pan)
    # pan_cn1, pan_cn2 = split_pan(pan, tree, cnoutgr)
    
    # give the subfamilies names 
    family = pan.orthogroup.tolist()[0]    
    pan1.loc[:, "orthogroup"] = family + "_1"
    pan2.loc[:, "orthogroup"] = family + "_2"
    
    return([pan1, pan2, sequences1, sequences2, tree1, tree2, seedmatrix1, 
        seedmatrix2])

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

def split_family_recursive_H_nl(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the H_nl strategy.

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

    pan1, pan2 = split_family_H_nl(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)

    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_H_nl(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_H_nl(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])

    return(pan)

def split_family_recursive_T_nl(pan, sequences, threads, dio_tmp):
    """Splits up a gene family using the T-nl strategy. 
    
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
    
    pan1, pan2 = split_family_T_nl(pan, sequences, threads, dio_tmp)
    split = assess_split(pan1, pan2, family)
    
    if split:

        genes1 = pan1.index.tolist()
        genes2 = pan2.index.tolist()
        sequences1 = [seq for seq in sequences if seq.id in genes1]
        sequences2 = [seq for seq in sequences if seq.id in genes2]
        pan1 = split_family_recursive_T_nl(pan1, sequences1, threads, dio_tmp)
        pan2 = split_family_recursive_T_nl(pan2, sequences2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_FH(pan, sequences, hclust, ficlin, min_reps, 
    max_reps, max_align, seedmatrix, threads, dio_tmp):
    """Splits up a gene family using the H or FH strategy. 
    
    See [https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.\
        hierarchy.to_tree.html#scipy.cluster.hierarchy.to_tree]. 
    
    Args:
        pan (DataFrame): A gene table for a single gene family, with the
            columns gene, genome and orthogroup. 
        sequences (list): A list with one SeqRecord object per row in pan, in 
            the same order. 
        hclust: An hclust object from scipy.hierarchy.cluster in the form of a 
            tree; the tip ids should be gene names. 
        finclin (bool): Should ficlin be used to pick representatives?
        min_reps (int): The minimum number of representatives to use. 
        max_reps (int): The maximum number of representatives to use. 
        seedmatrix (np.array): An array of identity values of genes to seeds.
        threads (int): The number of threads to use.
        dio_tmp (str): Folder to store temporary files.
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    family = pan.orthogroup.tolist()[0]
    
    if not split_possible(pan.genome.tolist()):
        # logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)
    
    pan1, pan2, sequences1, sequences2, hclust1, hclust2, seedmatrix1, \
        seedmatrix2 = split_family_FH(pan, sequences, hclust, ficlin, min_reps,
        max_reps, max_align, seedmatrix, threads, dio_tmp)
    if pan1 is None:
        split = False # don't split if all sequences are identical
    else:
        split = assess_split(pan1, pan2, family)
    
    if split:

        pan1 = split_family_recursive_FH(pan1, sequences1, hclust1, ficlin, 
            min_reps, max_reps, max_align, seedmatrix1, threads, dio_tmp)
        pan2 = split_family_recursive_FH(pan2, sequences2, hclust2, ficlin, 
            min_reps, max_reps, max_align, seedmatrix2, threads, dio_tmp)
        pan = pd.concat([pan1, pan2])
        
    return(pan)

def split_family_recursive_FT(pan, sequences, tree, ficlin, min_reps, 
    max_reps, seedmatrix, threads, dio_tmp):
    """Splits up a gene family using the T or FT strategy. 
    
    Args:
        pan (DataFrame): A gene table for a single gene family, with the
            columns gene, genome and orthogroup. 
        sequences (list): A list with one SeqRecord object per row in pan, in 
            the same order. 
        tree: An ete4 tree object. 
        finclin (bool): Should ficlin be used to pick representatives?
        min_reps (int): The minimum number of representatives to use. 
        max_reps (int): The maximum number of representatives to use. 
        seedmatrix (np.array): An array of identity values of genes to seeds.
        threads (int): The number of threads to use.
        dio_tmp (str): Folder to store temporary files.
        
    Returns: 
        An pan object where the orthogroup column has been updated. 
    """
    
    family = pan.orthogroup.tolist()[0]
    
    if not split_possible(pan.genome.tolist()):
        # logging.info(f"{family}: {len(pan.index)} genes - split not an option")
        return(pan)
    
    pan1, pan2, sequences1, sequences2, tree1, tree2, seedmatrix1, \
        seedmatrix2 = split_family_FT(pan, sequences, tree, ficlin, min_reps, 
        max_reps, seedmatrix, threads, dio_tmp)
    if pan1 is None:
        split = False # don't split if all sequences are identical
    else:
        split = assess_split(pan1, pan2, family)
    
    if split:

        pan1 = split_family_recursive_FT(pan1, sequences1, tree1, ficlin, 
            min_reps, max_reps, seedmatrix1, threads, dio_tmp)
        pan2 = split_family_recursive_FT(pan2, sequences2, tree2, ficlin, 
            min_reps, max_reps, seedmatrix2, threads, dio_tmp)
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

def split_superfamily(pan, strategy, din_fastas, min_reps, max_reps, max_align, 
    threads, dio_tmp):
    """Splits a gene superfamily in a set of gene families.
    
    Args:
        pan (Data Frame): Table with columns gene, genome and orthogroup 
            for the genes or a single gene family. 
        strategy (str): Splitting strategy. [H-nl, T-nl, H, LH, HT, LHT, P]
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
    
    if strategy == "H-nl":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_H_nl(pan, sequences, threads, dio_tmp)
    elif strategy == "T-nl":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_T_nl(pan, sequences, threads, dio_tmp)
    elif strategy in ["H", "FH", "T", "FT"]:
        ficlin = "F" in strategy
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        genes = pan.index.tolist()
        sequences.sort(key = lambda s: genes.index(s.id))
        pan["rep"] = "c" 
        if "F" in strategy:
            seedmatrix = np.zeros([len(sequences), max_reps])
        else:
            seedmatrix = None
        if "H" in strategy:
            pan = split_family_recursive_FH(pan, sequences, None, ficlin,
                min_reps, max_reps, max_align, seedmatrix, threads, dio_tmp)
        else:
            pan = split_family_recursive_FT(pan, sequences, None, ficlin,
                min_reps, max_reps, seedmatrix, threads, dio_tmp)
        pan = pan[["genome", "orthogroup"]] # remark: "gene" is the index
    elif strategy == "P":
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        pan = split_family_recursive_P(pan, sequences, threads, dio_tmp)
    else:
        logging.error("pangenome strategy unknown")
            
    pan = pan.reset_index()
    shutil.rmtree(dio_tmp)
    
    print("|", end = "", flush = True)
    return(pan)

def infer_superfamilies(faafins, dout, speciesmode, threads):
    """Infers superfamilies of a set of faa files. 
    
    Remarks: the pangenome with the superfamilies is written out to 
    f"{dout}/pangenome.tsv".
    
    Args:
        faafins (list): Paths to .faa fasta files, one per genome. 
        dout (str): Output folder for the pangenome with the superfamilies.
        speciesmode (bool): Whether to run in species mode.
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
    sequence_db_path =  f"{dout}/sequenceDB"
    if len(faafins) > MAX_FAA_ARG_LEN:
        faafins = archive_faafins(faafins, sequence_db_path, f"{dout}/logs")
    run_mmseqs(["createdb"] + faafins + [f"{sequence_db_path}/db"], 
        f"{dout}/logs/createdb.log")
    
    # create preclusters with mmseqs cluster module (includes mmseqs linclust)
    logging.info("creating preclusters")
    min_seq_id = "0.7" if speciesmode else "0.2"
    run_mmseqs(["cluster", f"{dout}/sequenceDB/db", f"{dout}/preclusterDB/db",
        f"{dout}/tmp", "--max-seqs", "1000000", "-c", "0.5", "--cov-mode", "0",
        "-e", "inf", "--min-seq-id", min_seq_id, "--cluster-mode", "0"],
        f"{dout}/logs/cluster.log", 
        skip_if_exists = f"{dout}/preclusterDB/db.index", threads = threads)
    # why --full-header option? --> to avoid MMseqs2 extracting the
    # UniqueIdentifier part of sequences in UniProtKB format 
    # (see https://www.uniprot.org/help/fasta-headers)
    run_mmseqs(["createtsv", f"{dout}/sequenceDB/db", f"{dout}/sequenceDB/db",
        f"{dout}/preclusterDB/db", f"{dout}/preclusters.tsv", "--full-header"], 
        f"{dout}/logs/createtsv_preclusters.log")

    # cluster the preclusters into the final clusters
    if not speciesmode:
        logging.info("clustering the preclusters using profiles")
        run_mmseqs(["result2profile", f"{dout}/sequenceDB/db",
            f"{dout}/sequenceDB/db", f"{dout}/preclusterDB/db",
            f"{dout}/profileDB/db"], f"{dout}/logs/result2profile.log",
            skip_if_exists = f"{dout}/profileDB/db.index", threads = threads)
        run_mmseqs(["profile2consensus", f"{dout}/profileDB/db",
            f"{dout}/profileDB/db_consensus"],
            f"{dout}/logs/profile2consensus.log")
        run_mmseqs(["search", f"{dout}/profileDB/db",
            f"{dout}/profileDB/db_consensus", f"{dout}/alignmentDB/db",
            f"{dout}/tmp", "-c", "0.5", "--cov-mode", "1"],
            f"{dout}/logs/search.log",
            skip_if_exists = f"{dout}/alignmentDB/db.index", threads = threads)
        run_mmseqs(["clust", f"{dout}/profileDB/db", f"{dout}/alignmentDB/db",
            f"{dout}/clusterDB/db", "--cluster-mode", "2"],
            f"{dout}/logs/clust.log",
            skip_if_exists = f"{dout}/clusterDB/db.index", threads = threads)
        run_mmseqs(["createtsv", f"{dout}/sequenceDB/db",
            f"{dout}/sequenceDB/db", f"{dout}/clusterDB/db",
            f"{dout}/clusters.tsv", "--full-header"],
            f"{dout}/logs/createtsv_clusters.log")

    # create the pangenome file
    logging.info("compiling pangenome file with superfamilies")
    preclustertable = read_mmseqs_clustertable(f"{dout}/preclusters.tsv")
    preclustertable = preclustertable.rename(columns = {"cluster": "precluster"})
    if speciesmode:
        genes = preclustertable.rename(columns = {"precluster": "orthogroup"})
    else:
        clustertable = read_mmseqs_clustertable(f"{dout}/clusters.tsv")
        clustertable = clustertable.rename(columns = {"gene": "precluster"})
        genes = pd.merge(preclustertable, clustertable, on = "precluster")
        genes = genes.rename(columns = {"cluster": "orthogroup"})
        genes = genes.drop(["precluster"], axis = 1)
    
    # rename the superfamilies 
    famnames_old = genes["orthogroup"].unique()
    famnames_new = [f"F{c}" for c in padded_counts(len(famnames_old))]
    namedict = dict(zip(famnames_old, famnames_new))
    genes["orthogroup"] = [namedict[f] for f in genes["orthogroup"]]
    
    # write pangenome file
    write_tsv(genes, f"{dout}/genes.tsv")

def infer_pangenome(faafins, splitstrategy, min_reps, max_reps, max_align,
    speciesmode, dout, threads):
    """Infers the pangenome of a set of genomes and writes it to disk. 
    
    Args:
        faafins (list): Paths to .faa fasta files, one per genome. 
        splitstrategy (str): Splitting strategy.
        speciesmode (bool): Whether to run in species mode.
        dout (str): Output folder for the pangenome.
        threads (int): Number of threads to use.
    """
  
    # threads per process
    tpp = 1
    
    logging.info(f"{len(faafins)} genomes were supplied")
    
    logging.info("constructing gene table")
    genes = extract_genes(faafins)
    
    logging.info("checking if gene names are unique")
    if not genes.gene.is_unique:
        duplicated_genes = genes.gene.duplicated()
        genomes_with_duplication = genes.genome[duplicated_genes].unique()
        logging.error("gene names are not unique\n the following genomes " +
            f"have duplicated genes: {genomes_with_duplication}") 
        sys.exit(1)
    
    if (len(faafins)) == 1:
      
        logging.info("only one genome supplied - each gene will be its own "
            "orthogroup")

        logging.info("assiging names to gene families")
        pangenome = genes
        pangenome["orthogroup"] = [f"F{c}" for c in \
            padded_counts(len(pangenome.index))]
    
        logging.info("writing pangenome file")
        write_tsv(pangenome, f"{dout}/pangenome.tsv")
        
        return()
  
    logging.info("STAGE 1: creation of superfamilies")
    os.makedirs(f"{dout}/superfamilies", exist_ok = True)
    infer_superfamilies(faafins, f"{dout}/superfamilies", speciesmode, threads)
    genes_superfams = pd.read_csv(f"{dout}/superfamilies/genes.tsv", 
        sep = "\t", names = ["gene", "orthogroup"])
    pangenome = pd.merge(genes, genes_superfams, on = "gene")
    
    if splitstrategy == "S":
      
        logging.info("writing pangenome file")
        write_tsv(pangenome, f"{dout}/pangenome.tsv")
        logging.info("removing temporary folders")
        shutil.rmtree(f"{dout}/superfamilies")
        return()
      
    logging.info("STAGE 2: splitting of superfamilies")

    logging.info("gathering sequences of superfamilies")
    os.makedirs(f"{dout}/superfamilies/superfamilies", exist_ok = True)
    dio_fastas = f"{dout}/superfamilies/fastas"
    gather_orthogroup_sequences(pangenome, faafins, dio_fastas)

    logging.info("counting splitable superfamilies")
    os.makedirs(f"{dout}/tmp", exist_ok = True)
    pangenome = pangenome.groupby("orthogroup")
    orthogroups = pangenome.aggregate({"genome": split_possible})
    splitable = orthogroups[orthogroups.genome].index.tolist()
    pangenome_splitable = [pan for name, pan in pangenome if name in splitable]
    pangenome_splitable.sort(key = lambda pan: len(pan.index), reverse = True)
    n = len(pangenome_splitable)
    logging.info(f"found {n} splitable superfamilies")
    
    if n == 0:
  
        logging.info("writing pangenome file")
        pangenome = pangenome.obj # to "ungroup"
        write_tsv(pangenome, f"{dout}/pangenome.tsv")
        logging.info("removing temporary folders")
        shutil.rmtree(f"{dout}/superfamilies")
        shutil.rmtree(f"{dout}/tmp")
        return()
    
    logging.info("splitting superfamilies")
    with ProcessPoolExecutor(max_workers = threads // tpp) as executor:
        pangenome_splitable = executor.map(split_superfamily, 
            pangenome_splitable, [splitstrategy] * n, [dio_fastas] * n, 
            [min_reps] * n, [max_reps] * n, [max_align] * n, [tpp] * n, 
            [f"{dout}/tmp"] * n)
    print("")
    pangenome = pd.concat(list(pangenome_splitable) +
        [pan for name, pan in pangenome if not name in splitable])
    n_orthogroups = pangenome["orthogroup"].nunique()
    logging.info(f"Identified {n_orthogroups} orthogroups")

    logging.info("Determining single-copy core orthogroups")
    n_genomes = pangenome["genome"].nunique()
    copynumbers = (
        pangenome
        .groupby(["genome", "orthogroup"])
        .size()
        .rename("copies")
        .reset_index()
    )
    core = (
        copynumbers.loc[copynumbers["copies"] == 1, "orthogroup"]
        .value_counts()
        .div(n_genomes)
        .loc[lambda p: p >= 0.95]
        .index
        .tolist()
    )
    logging.info(f"Identified {len(core)} 95% single-copy core orthogroups")

    logging.info("assigning names to the gene families")
    nametable = pd.DataFrame({"old": pangenome["orthogroup"].unique()})
    nametable["sf"] = [f.split("_")[0] for f in nametable["old"]]
    nametable["new"] = nametable.groupby("sf")["sf"].transform(lambda sfs: \
        [f"{sfs.tolist()[0]}_{c}" for c in padded_counts(len(sfs))])
    namedict = dict(zip(nametable["old"], nametable["new"]))
    pangenome["orthogroup"] = [namedict[f] for f in pangenome["orthogroup"]]
    
    logging.info("writing pangenome file")
    write_tsv(pangenome, f"{dout}/pangenome.tsv")
    
    logging.info("removing temporary folders")
    shutil.rmtree(f"{dout}/superfamilies")
    shutil.rmtree(f"{dout}/tmp")
