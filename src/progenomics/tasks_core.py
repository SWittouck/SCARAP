import logging
import os
import shutil

from argparse import Namespace
from Bio import SeqIO, AlignIO, Align, SeqRecord
from copy import copy
from ete3 import Tree
from concurrent.futures import ProcessPoolExecutor
from scipy import cluster

from utils import *
from readerswriters import *
from computers import *
from callers import *

def infer_superfamilies(faafins, dout, threads):
    """Infers superfamilies of a set of faa files. 
    
    Args:
        faafins (list): Paths to .faa fasta files, one per genome. 
        dout (chr): Output folder for the pangenome with the superfamilies.
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
        f"{dout}/tmp", "-s", "7.5", "--max-seqs", "1000000", "-c", "0.5", 
        "-e", "inf"], f"{dout}/logs/cluster.log", 
        skip_if_exists = f"{dout}/preclusterDB/db.index", threads = threads)

    # cluster the preclusters into the final clusters
    logging.info("clustering the preclusters")
    run_mmseqs(["result2profile", f"{dout}/sequenceDB/db",
        f"{dout}/sequenceDB/db", f"{dout}/preclusterDB/db",
        f"{dout}/profileDB/db"], f"{dout}/logs/result2profile.log",
        skip_if_exists = f"{dout}/profileDB/db.index", threads = threads)
    run_mmseqs(["search", f"{dout}/profileDB/db",
        f"{dout}/profileDB/db_consensus", f"{dout}/alignmentDB/db",
        f"{dout}/tmp", "-s", "7.5"], f"{dout}/logs/search.log",
        skip_if_exists = f"{dout}/alignmentDB/db.index", threads = threads)
    run_mmseqs(["clust", f"{dout}/profileDB/db", f"{dout}/alignmentDB/db",
        f"{dout}/clusterDB/db"], f"{dout}/logs/clust.log",
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
    
def generate_subfams_treestrat(fin_seqs, dio_tmp, threads):
    """Generates two subfamilies for a gene family.
    
    Generates subfamilies using a tree-based strategy. 
    
    Args:
        fin_seqs (chr): Path to fasta file with sequences of genes.
        dio_tmp (chr): Path to folder to store temporary files in.
        threads (int): Number of threads to use.
        
    Returns:
        A list of two elements: the gene ids for subfamily one and the gene ids
            for subfamily two. 
    """
    run_mafft(fin_seqs, f"{dio_tmp}/seqs.aln", threads)
    # infer tree with iqtree
    run_iqtree(f"{dio_tmp}/seqs.aln", f"{dio_tmp}/tree", threads, 
        ["-m", "LG"])
    tree = Tree(f"{dio_tmp}/tree/tree.treefile")
    # midpoint root the tree and split in two
    midoutgr = tree.get_midpoint_outgroup()
    # if midoutgr is None:
    #     logging.info(f"{family}: failed to midpoint root - moving on")
    #     return(pangenome)
    genes_subfam1 = midoutgr.get_leaf_names()
    midoutgr.detach()
    genes_subfam2 = tree.get_leaf_names()
    return([genes_subfam1, genes_subfam2])
    
def generate_subfams_cluststrat(fin_seqs, dio_tmp, threads):
    """Generates two subfamilies for a gene family.
    
    Generates subfamilies using a hclust-based strategy. 
    
    Args:
        fin_seqs (chr): Path to fasta file with sequences of genes.
        
    Returns:
        A list of two elements: the gene ids for subfamily one and the gene ids
            for subfamily two. 
    """
    run_mafft(fin_seqs, f"{dio_tmp}/seqs.aln", threads = threads)
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
    # tree = to_tree(link)
    return([genes_subfam1, genes_subfam2])
        
def generate_subfams_profilestrat(fin_seqs, dio_tmp, threads):
    """Generates two subfamilies for a gene family.
    
    Generates subfamilies using a profile-based strategy. 
    
    Args:
        fin_seqs (chr): Path to fasta file with sequences of genes.
        dio_tmp (chr): Path to folder to store temporary files in.
        threads (int): Number of threads to use.
        
    Returns:
        A list of two elements: the gene ids for subfamily one and the gene ids
            for subfamily two. 
    """
  
    report = False
        
    # align sequences
    run_mafft(fin_seqs, f"{dio_tmp}/seqs.aln", threads)

    # create sequence database
    for dir in ["sequenceDB", "tmp", "resultDB", "logs"]:
        makedirs_smart(f"{dio_tmp}/{dir}")
    run_mmseqs(["createdb", fin_seqs, f"{dio_tmp}/sequenceDB/db"], 
        f"{dio_tmp}/logs/createdb.log")
    
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
    # print(f"{len(genes_profile1)}, {len(genes_profile2)}")
    return([genes_profile1, genes_profile2])
    
def split_family_recursive(pangenome, sequences, family, strategy, threads, 
    dio_tmp):
    """Implements the main functionality of split_family through recursion.
    """
    
    # return pangenome if split is not possible
    genomes = pangenome[pangenome.orthogroup == family].genome.tolist()
    if not split_possible(genomes):
        logging.info(f"{family}: {len(genomes)} genes - split not an option")
        return(pangenome)
    
    # write sequences of family to file
    genes = pangenome[pangenome.orthogroup == family].index.tolist()
    seqs = [sequences[gene] for gene in genes]
    write_fasta(seqs, f"{dio_tmp}/seqs.fasta")

    # suggest subfamilies
    if strategy == "trees":
        genes_subfam1, genes_subfam2 = generate_subfams_treestrat(
            f"{dio_tmp}/seqs.fasta", f"{dio_tmp}", threads)
    elif strategy == "profiles":
        genes_subfam1, genes_subfam2 = generate_subfams_profilestrat(
            f"{dio_tmp}/seqs.fasta", f"{dio_tmp}", threads)
    elif strategy == "clusters":
        genes_subfam1, genes_subfam2 = generate_subfams_cluststrat(
            f"{dio_tmp}/seqs.fasta", f"{dio_tmp}", threads)
    else:
        loggig.error("unknown splitting strategy")
        
    # apply splitting criterium
    genomes_subfam1 = pangenome.loc[genes_subfam1, "genome"].tolist()
    genomes_subfam2 = pangenome.loc[genes_subfam2, "genome"].tolist()
    pgo_obs = calc_pgo(genomes_subfam1, genomes_subfam2)
    pgo_exp = pred_pgo(genomes_subfam1, genomes_subfam2)
    split = pgo_obs >= pgo_exp and not pgo_exp == 0
    n1 = len(set(genomes_subfam1))
    n2 = len(set(genomes_subfam2))
    
    # split if necessary
    if split:
        logging.info(f"{family}: {n1}/{n2}; "
            f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
            f"decision = split")
        # replace family ids with newly generated subfamily ids
        subfam_1 = family + "_1"
        subfam_2 = family + "_2"
        pangenome.loc[genes_subfam1, "orthogroup"] = subfam_1
        pangenome.loc[genes_subfam2, "orthogroup"] = subfam_2
        # attempt to split the two subfamilies
        pangenome = split_family_recursive(pangenome, sequences, subfam_1, 
            strategy, threads, dio_tmp)
        pangenome = split_family_recursive(pangenome, sequences, subfam_2, 
            strategy, threads, dio_tmp)
    else:
        logging.info(f"{family}: {n1}/{n2}; "
            f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
            f"decision = no split")
            
    return(pangenome)
    
def split_family_recursive_clust(pangenome, family, hclust, genes):
    
    # return pangenome if split is not possible
    genomes = pangenome[pangenome.orthogroup == family].genome.tolist()
    if not split_possible(genomes):
        logging.info(f"{family}: {len(genomes)} genes - split not an option")
        return(pangenome)
  
    hclust1 = hclust.get_left()
    hclust2 = hclust.get_right()
    
    genes_subfam1 = [genes[i] for i in hclust1.pre_order(lambda x: x.id)]
    genes_subfam2 = [genes[i] for i in hclust2.pre_order(lambda x: x.id)]
        
    # apply splitting criterium
    genomes_subfam1 = pangenome.loc[genes_subfam1, "genome"].tolist()
    genomes_subfam2 = pangenome.loc[genes_subfam2, "genome"].tolist()
    pgo_obs = calc_pgo(genomes_subfam1, genomes_subfam2)
    pgo_exp = pred_pgo(genomes_subfam1, genomes_subfam2)
    split = pgo_obs >= pgo_exp and not pgo_exp == 0
    n1 = len(set(genomes_subfam1))
    n2 = len(set(genomes_subfam2))
    
    # split if necessary
    if split:
        logging.info(f"{family}: {n1}/{n2}; "
            f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
            f"decision = split")
        # replace family ids with newly generated subfamily ids
        subfam1 = family + "_1"
        subfam2 = family + "_2"
        pangenome.loc[genes_subfam1, "orthogroup"] = subfam1
        pangenome.loc[genes_subfam2, "orthogroup"] = subfam2
        # attempt to split the two subfamilies
        pangenome = split_family_recursive_clust(pangenome, subfam1, hclust1,
            genes)
        pangenome = split_family_recursive_clust(pangenome, subfam2, hclust2,
            genes)
    else:
        logging.info(f"{family}: {n1}/{n2}; "
            f"pgo_obs = {pgo_obs:.3f}; pgo_exp = {pgo_exp:.3f}; "
            f"decision = no split")
    
    return(pangenome)
    
def split_family(pangenome, threads, din_fastas, strategy, dio_tmp):
    """Splits a gene family in one, two or more subfamilies.
    
    Args:
        pangenome (Data Frame): Table with columns gene, genome and orthogroup 
            for the genes or a single gene family. 
        threads (int): Number of threads to use. 
        din_fastas (chr): Input folder with fasta file of gene family. 
        strategy (chr): Splitting strategy. ["trees", "profiles", "clusters", 
            "clusters_fast"]
        dio_tmp (chr): Folder to store temporary files.
        
    Returns:
        A pandas Data Frame with the pangenome, where the orthogroups have been
            split. 
    """
    
    superfam = pangenome.orthogroup.tolist()[0]
    pangenome = pangenome.set_index("gene")
    dio_tmp = f"{dio_tmp}/{superfam}"
    makedirs_smart(dio_tmp)
    
    if strategy == "clusters_fast": 
        run_mafft(f"{din_fastas}/{superfam}.fasta", f"{dio_tmp}/seqs.aln", 
            threads)
        with open(f"{dio_tmp}/seqs.aln", "r") as fin:
            aln = AlignIO.read(fin, "fasta")
        genes = [seq.id for seq in aln]
        n = len(aln)
        im = identity_matrix(aln)
        dm = []
        for r in range(n - 1):
            for c in range(r + 1, n):
                dm.append(1 - im[r, c])
        hclust = cluster.hierarchy.linkage(dm, method = "average")
        hclust = cluster.hierarchy.to_tree(hclust)
        pangenome = split_family_recursive_clust(pangenome, superfam, hclust, 
            genes)
        
    else:
        sequences = read_fasta(f"{din_fastas}/{superfam}.fasta")
        sequences = {record.id: record for record in sequences}
        pangenome = split_family_recursive(pangenome, sequences, superfam, 
            strategy, threads = threads, dio_tmp = dio_tmp)
            
    pangenome = pangenome.reset_index()
    shutil.rmtree(dio_tmp)
    return(pangenome)

def infer_pangenome(faafins, dout, splitstrategy, threads):
    """Infers the pangenome of a set of genomes and writes it to disk. 
    
    Args:
        faafins (list): Paths to .faa fasta files, one per genome. 
        dout (chr): Output folder for the pangenome.
        splitstrategy (chr): Splitting strategy. ["trees", "profiles"]
        threads (int): Number of threads to use.
    """
  
    # threads per process
    tpp = 2
  
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
        pangenome_splitable = executor.map(split_family, pangenome_splitable,
            [tpp] * n, [dio_fastas] * n, [splitstrategy] * n,
            [f"{dout}/tmp"] * n)
    pangenome = pd.concat(list(pangenome_splitable) +
        [pan for name, pan in pangenome if not name in splitable])
    write_tsv(pangenome, f"{dout}/pangenome.tsv")

def run_pan_nonhier(args):

    pangenomefout = os.path.join(args.outfolder, "pangenome.tsv")
    if os.path.isfile(pangenomefout):
        logging.info("existing pangenome detected - moving on")
        return()
    
    faafins = read_lines(args.faapaths)
    
    if args.method in ["of_blast", "of_diamond"]:

        logging.info("creating orthofinder directory")
        dir_orthofinder = os.path.join(args.outfolder, "orthofinder")
        os.makedirs(dir_orthofinder, exist_ok = True)

        logging.info("running orthofinder")
        logfile = os.path.join(args.outfolder, "orthofinder.log")
        engine = args.method.split("_")[1]
        run_orthofinder(faafins, dir_orthofinder, logfile, args.threads, 
            engine)
    
        logging.info("creating tidy pangenome file")
        pangenome = read_pangenome_orthofinder(dir_orthofinder)
        write_tsv(pangenome, pangenomefout)
        
    elif args.method == "builtin_profiles":
        
        logging.info("constructing pangenome with profile strategy")
        infer_pangenome(faafins, args.outfolder, "profiles", args.threads)
        
    elif args.method == "builtin_trees":
        
        logging.info("constructing pangenome with tree strategy")
        infer_pangenome(faafins, args.outfolder, "trees", args.threads)
        
    elif args.method == "builtin_clusters":
        
        logging.info("constructing pangenome with cluster strategy")
        infer_pangenome(faafins, args.outfolder, "clusters", args.threads)
        
    elif args.method == "builtin_clusters_fast":
        
        logging.info("constructing pangenome with cluster strategy")
        infer_pangenome(faafins, args.outfolder, "clusters_fast", args.threads)

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
            threads = args.threads, method = "of_diamond"))
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
        method = "of_blast"))

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

def run_pan(args):
    if "species" in args and not args.species is None:
        run_pan_hier(args)
    else:
        run_pan_nonhier(args)

def run_build(args):

    if os.path.isfile(os.path.join(args.outfolder, "hmm_db.h3f")):
        logging.info("existing database detected - moving on")
        return()

    logging.info("creating output subfolders")
    orthogroupsdio = os.path.join(args.outfolder, "orthogroups")
    alignmentsdio = os.path.join(args.outfolder, "alignments")
    profilesdio = os.path.join(args.outfolder, "profiles")
    os.makedirs(orthogroupsdio, exist_ok = True)
    os.makedirs(alignmentsdio, exist_ok = True)
    os.makedirs(profilesdio, exist_ok = True)

    logging.info("gathering sequences of orthogroups")
    pangenome = read_genes(args.pangenome)
    faafins = read_lines(args.faapaths)
    gather_orthogroup_sequences(pangenome, faafins, orthogroupsdio,
        args.min_genomes)
    logging.info(f"gathered sequences for {len(os.listdir(orthogroupsdio))} "
        f"orthogroups occurring in at least {args.min_genomes} genome(s)")

    logging.info("aligning orthogroups")
    orthogroups = [os.path.splitext(file)[0] for file in
        os.listdir(orthogroupsdio)]
    orthogroupfouts = make_paths(orthogroups, orthogroupsdio, ".fasta")
    alifouts = make_paths(orthogroups, alignmentsdio, ".aln")
    run_mafft_parallel(orthogroupfouts, alifouts)

    logging.info("building profile hmms")
    profilefouts = make_paths(orthogroups, profilesdio, ".hmm")
    run_hmmbuild_parallel(alifouts, profilefouts)

    logging.info("pressing profile hmm database")
    run_hmmpress(profilefouts, args.outfolder)

def run_search(args):

    genesfout = os.path.join(args.outfolder, "genes.tsv")
    if os.path.isfile(genesfout):
        logging.info("existing search results detected - moving on")
        return()

    cutoffsfio = os.path.join(args.db, "orthogroups.tsv")

    logging.info("performing hmmsearch")
    queryfins = read_lines(args.qpaths)
    domtblfio = os.path.join(args.outfolder, "hmmer_domtbl.tmp")
    if os.path.isfile(domtblfio):
        logging.info("existing hmmer domtbl detected - skipping hmmsearch")
    else:
        run_hmmsearch(args.db, queryfins, domtblfio)
    hits = read_domtbl(domtblfio)
    # write_tsv(hits, os.path.join(args.outfolder, "hits.tsv"))
    genes_genomes = extract_genes(queryfins)

    if args.trainstrategy == "pan":

        logging.info("traninig profile-specific hmmer cutoffs with the pan "
            "strategy")
        pangenome = read_genes(args.pangenome)
        cutoffs = train_cutoffs_pan(hits, pangenome)
        write_tsv(cutoffs, cutoffsfio)

    elif args.trainstrategy == "core":

        logging.info("traninig profile-specific hmmer cutoffs with the core "
            "strategy")
        cutoffs = train_cutoffs_core(hits, genes_genomes)
        write_tsv(cutoffs, cutoffsfio)

    logging.info("applying hmmer score cutoffs")
    orthogroups = read_orthogroups(cutoffsfio)
    genes = process_scores(hits, orthogroups)
    genes = pd.merge(genes, genes_genomes, how = "left")
    genes = genes[["gene", "genome", "orthogroup"]]
    write_tsv(genes, genesfout)

    logging.info("removing temporary files")
    os.remove(domtblfio)

def run_checkgenomes(args):

    logging.info("checking genomes")
    coregenome = read_genes(args.coregenome)
    genomes = checkgenomes(coregenome)
    write_tsv(genomes, os.path.join(args.outfolder, "genomes.tsv"))

def run_checkgroups(args):

    logging.info("checking core orthogroups")
    coregenome = read_genes(args.coregenome)
    orthogroups = checkgroups(coregenome)
    write_tsv(orthogroups, os.path.join(args.outfolder, "orthogroups.tsv"))

def run_filter(args):

    logging.info("reading pangenome")
    pangenome = read_genes(args.pangenome)
    if not args.genomes is None:
        logging.info("filtering genomes")
        genomes = read_lines(args.genomes)
        pangenome = filter_genomes(pangenome, genomes)
    if not args.orthogroups is None:
        logging.info("filtering orthogroups")
        orthogroups = read_lines(args.orthogroups)
        pangenome = filter_groups(pangenome, orthogroups)
    write_tsv(pangenome, os.path.join(args.outfolder, "pangenome.tsv"))

def run_supermatrix(args):

    sm_aas_fout = os.path.join(args.outfolder, "supermatrix_aas.fasta")
    sm_nucs_fout = os.path.join(args.outfolder, "supermatrix_nucs.fasta")
    seqs_aas_dio = os.path.join(args.outfolder, "seqs_aas")
    seqs_nucs_dio = os.path.join(args.outfolder, "seqs_nucs")
    alis_aas_dio = os.path.join(args.outfolder, "alis_aas")
    alis_nucs_dio = os.path.join(args.outfolder, "alis_nucs")

    logging.info("creating output subfolders")
    os.makedirs(seqs_aas_dio, exist_ok = True)
    os.makedirs(alis_aas_dio, exist_ok = True)

    if os.path.isfile(sm_aas_fout):

        logging.info("existing amino acid supermatrix detected - moving on")

    else:

        logging.info("gathering amino acid sequences of orthogroups")
        coregenome = read_genes(args.coregenome)
        faa_fins = read_lines(args.faapaths)
        gather_orthogroup_sequences(coregenome, faa_fins, seqs_aas_dio)
        logging.info(f"gathered sequences for {len(os.listdir(seqs_aas_dio))} "
            f"core orthogroups")

        logging.info("aligning orthogroups on the amino acid level")
        orthogroups = [os.path.splitext(file)[0] for file in
            os.listdir(seqs_aas_dio)]
        seqs_aas_fios = make_paths(orthogroups, seqs_aas_dio, ".fasta")
        alis_aas_fios = make_paths(orthogroups, alis_aas_dio, ".aln")
        run_mafft_parallel(seqs_aas_fios, alis_aas_fios)

        logging.info("concatenating amino acid alignments")
        construct_supermatrix(coregenome, alis_aas_fios, sm_aas_fout)

    if not args.ffnpaths is None and os.path.isfile(sm_nucs_fout):

        logging.info("existing nucleotide supermatrix detected - moving on")

    if not args.ffnpaths is None and not os.path.isfile(sm_nucs_fout):

        logging.info("creating output subfolders")
        os.makedirs(seqs_nucs_dio, exist_ok = True)
        os.makedirs(alis_nucs_dio, exist_ok = True)

        logging.info("gathering nucleotide sequences of orthogroups")
        coregenome = read_genes(args.coregenome)
        ffn_fins = read_lines(args.ffnpaths)
        gather_orthogroup_sequences(coregenome, ffn_fins, seqs_nucs_dio)

        logging.info("aligning orthogroups on the nucleotide level")
        orthogroups = [os.path.splitext(file)[0] for file in
            os.listdir(seqs_aas_dio)]
        seqs_nucs_fios = make_paths(orthogroups, seqs_nucs_dio, ".fasta")
        alis_aas_fios = make_paths(orthogroups, alis_aas_dio, ".aln")
        alis_nucs_fios = make_paths(orthogroups, alis_nucs_dio, ".aln")
        reverse_align_parallel(seqs_nucs_fios, alis_aas_fios, alis_nucs_fios)

        logging.info("concatenating nucleotide alignments")
        construct_supermatrix(coregenome, alis_nucs_fios, sm_nucs_fout)
