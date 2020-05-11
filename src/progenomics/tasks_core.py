import logging
import os
import shutil

from argparse import Namespace
from ete3 import Tree
from concurrent.futures import ProcessPoolExecutor

from utils import *
from readerswriters import *
from computers import *
from callers import *

def infer_superfamilies(faafins, dout, threads):
    
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
        f"{dout}/tmp", "--min-seq-id", "0.1", "-c", "0.5", "--threads", 
        threads], f"{dout}/logs/cluster.log", 
        skip_if_exists = f"{dout}/preclusterDB/db.index")

    # cluster the preclusters into the final clusters
    logging.info("clustering the preclusters")
    run_mmseqs(["result2profile", f"{dout}/sequenceDB/db", 
        f"{dout}/sequenceDB/db", f"{dout}/preclusterDB/db", 
        f"{dout}/profileDB/db", "--threads", threads], 
        f"{dout}/logs/result2profile.log", 
        skip_if_exists = f"{dout}/profileDB/db.index")
    run_mmseqs(["search", f"{dout}/profileDB/db", 
        f"{dout}/profileDB/db_consensus", f"{dout}/alignmentDB/db", 
        f"{dout}/tmp", "--threads", threads, "-s", "7.5"], 
        f"{dout}/logs/search.log", 
        skip_if_exists = f"{dout}/alignmentDB/db.index")
    run_mmseqs(["clust", f"{dout}/profileDB/db", f"{dout}/alignmentDB/db",
        f"{dout}/clusterDB/db", "--threads", threads], 
        f"{dout}/logs/clust.log", 
        skip_if_exists = f"{dout}/clusterDB/db.index")

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
    
def split_family(pangenome, sequences, family, threads, dio_tmp):
    # if family is fully single-copy: just return pangenome
    genomes = pangenome[pangenome.orthogroup == family].genome.tolist()
    if len(genomes) <= 3 or len(set(genomes)) == 1:
        logging.info(f"{family} has not enough copies/genomes - moving on")
        return(pangenome)
    if len(genomes) == len(set(genomes)):
        logging.info(f"{family} is single copy - moving on")
        return(pangenome)
    # write sequences of family to file
    genes = pangenome[pangenome.orthogroup == family].index.tolist()
    seqs = [sequences[gene] for gene in genes]
    write_fasta(seqs, f"{dio_tmp}/seqs.fasta")
    # align sequences of family
    run_mafft(f"{dio_tmp}/seqs.fasta", f"{dio_tmp}/seqs.aln", threads)
    # infer tree of family
    run_iqtree(f"{dio_tmp}/seqs.aln", f"{dio_tmp}/tree", threads, 
        ["-m", "LG"])
    # read tree of family 
    tree = Tree(f"{dio_tmp}/tree/tree.treefile")
    # midpoint root the tree and split in two
    midoutgr = tree.get_midpoint_outgroup()
    if midoutgr is None:
        logging.info(f"{family} failed to midpoint root - moving on")
        return(pangenome)
    genes_subfam1 = midoutgr.get_leaf_names()
    midoutgr.detach()
    genes_subfam2 = tree.get_leaf_names()
    # check if splitting is necessary
    genomes_subfam1 = pangenome.loc[genes_subfam1, "genome"].tolist()
    genomes_subfam2 = pangenome.loc[genes_subfam2, "genome"].tolist()
    split = decide_split(genomes_subfam1, genomes_subfam2)
    if split:
        logging.info(f"{family} will be split")
        # replace family ids with newly generated subfamily ids
        subfam_1 = family + "_1"
        subfam_2 = family + "_2"
        pangenome.loc[genes_subfam1, "orthogroup"] = subfam_1
        pangenome.loc[genes_subfam2, "orthogroup"] = subfam_2
        # attempt to split the two subfamilies
        pangenome = split_family(pangenome, sequences, subfam_1, threads, 
            dio_tmp)
        pangenome = split_family(pangenome, sequences, subfam_2, threads, 
            dio_tmp)
    else:
        logging.info(f"{family} will be not be split")
    return(pangenome)
    
def split_superfamily(pangenome, threads, dio_fastas, dio_tmp):
    superfam = pangenome.orthogroup.tolist()[0]
    pangenome = pangenome.set_index("gene")
    dio_tmp = f"{dio_tmp}/{superfam}"
    makedirs_smart(dio_tmp)
    sequences = read_fasta(f"{dio_fastas}/{superfam}.fasta")
    sequences = {record.id: record for record in sequences}
    pangenome = split_family(pangenome, sequences, superfam, 
        threads = threads, dio_tmp = dio_tmp)
    pangenome = pangenome.reset_index()
    shutil.rmtree(dio_tmp)
    return(pangenome)

def infer_pangenome(faafins, dout, threads):
  
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
    pangenome = [pan for name, pan in pangenome]
    n = len(pangenome)
    with ProcessPoolExecutor(max_workers = threads) as executor:
        pangenome = executor.map(split_superfamily, pangenome, [1] * n, 
        [dio_fastas] * n, [f"{dout}/tmp"] * n)
    pangenome = pd.concat(pangenome)
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
        run_orthofinder(faafins, dir_orthofinder, logfile, args.threads, engine)
    
        logging.info("creating tidy pangenome file")
        pangenome = read_pangenome_orthofinder(dir_orthofinder)
        write_tsv(pangenome, pangenomefout)
        
    elif args.method == "builtin":
        
        logging.info("constructing pangenome")
        infer_pangenome(faafins, args.outfolder, args.threads)

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
