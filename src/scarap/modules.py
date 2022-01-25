import logging
import numpy as np
import os
import pandas as pd
import shutil

from argparse import Namespace
from concurrent.futures import ProcessPoolExecutor
from statistics import median, mean
from random import sample

from callers import *
from computers import *
from helpers import * 
from pan import *
from readerswriters import *
from utils import *

def run_pan(args):
    if "species" in args and not args.species is None:
        run_pan_hier(args)
    else:
        run_pan_nonhier(args)

def run_pan_nonhier(args):

    pangenomefout = os.path.join(args.outfolder, "pangenome.tsv")
    if os.path.isfile(pangenomefout):
        logging.info("existing pangenome detected - moving on")
        return()
    
    faafins = read_fastapaths(args.faapaths)
    
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
      
        logging.info(f"pangenome will be constructed with the {args.method} "
            "strategy")
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
    faafins = read_fastapaths(args.faapaths)
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
        logging.info(f"started inferring pangenome of {species}")
        faapaths = [genomedict[genome] for genome in genomes]
        dout = os.path.join(speciespansdio, species)
        os.makedirs(dout, exist_ok = True)
        faapathsfio = os.path.join(dout, "faapaths.txt")
        write_lines(faapaths, faapathsfio)
        run_pan_nonhier(Namespace(faapaths = faapathsfio, outfolder = dout,
            threads = args.threads, method = args.method))
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
        method = args.method))

    logging.info("inflating metapangenome with species pangenomes")
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
    
def run_build(args): 
    
    fin_faapaths = args.faapaths
    fin_pangenome = args.pangenome
    dout = args.outfolder
    core_prefilter = args.core_prefilter
    core_filter = args.core_filter
    max_cores = args.max_cores
    threads = args.threads
    
    # define output paths/folders
    dout_ogseqs = os.path.join(dout, "orthogroups")
    dout_alis = os.path.join(dout, "alignments")
    dout_tmp = os.path.join(dout, "tmp")
    fout_cutoffs = os.path.join(dout, "cutoffs.csv")
    fout_hits = os.path.join(dout, "hits.tsv")
    fout_genes = os.path.join(dout, "genes.tsv")

    if os.path.isfile(fout_cutoffs):
        logging.info("existing database detected - moving on")
        return()

    logging.info("creating output subfolders")
    for dir in [dout_ogseqs, dout_alis]:
        os.makedirs(dir, exist_ok = True)
        
    logging.info("reading pangenome")
    pangenome = read_genes(fin_pangenome)
    # read separate table with genomes of all genes, in case the pangenome
    # table is not complete or core genes get selected
    fins_faas = read_fastapaths(fin_faapaths)
    genes_genomes = extract_genes(fins_faas)
        
    if core_prefilter != 0:
        logging.info(f"applying core prefilter of {core_prefilter}")
        corefams = determine_corefams(pangenome, core_filter)
        pangenome = filter_groups(pangenome, corefams)

    logging.info("gathering sequences of orthogroups")
    gather_orthogroup_sequences(pangenome, fins_faas, dout_ogseqs, 1)
    logging.info(f"gathered sequences for {len(os.listdir(dout_ogseqs))} "
        f"orthogroups")

    logging.info("aligning orthogroups")
    orthogroups = [os.path.splitext(f)[0] for f in os.listdir(dout_ogseqs)]
    fouts_ogseqs = make_paths(orthogroups, dout_ogseqs, ".fasta")
    fouts_alis = make_paths(orthogroups, dout_alis, ".aln")
    run_mafft_parallel(fouts_ogseqs, fouts_alis)
    
    # run profile search (function does its own logging)
    run_profilesearch(fins_faas, fouts_alis, fout_hits, dout_tmp, threads)
    
    logging.info("training score cutoffs for profiles") 
    colnames = ["gene", "profile", "score"]
    hits = pd.read_csv(fout_hits, sep = "\t", names = colnames, 
        usecols = [0, 1, 2])
    hits[["gene", "profile"]] = hits[["gene", "profile"]].\
        applymap(lambda x: x.split(" ")[0])
    cutoffs = train_cutoffs(hits, pangenome)
    
    if core_filter != 0 or max_cores != 0:
        logging.info(f"applying core filter of {core_filter} and maximum "
            f"number of core genes of {max_cores}")
        genes = process_scores(hits, cutoffs, top_profiles = False)
        genes = pd.merge(genes_genomes, genes, how = "right")
        corefams = determine_corefams(genes, core_filter = core_filter,
            max_cores = max_cores)
        hits = hits[hits["profile"].isin(corefams)]
        cutoffs = cutoffs[cutoffs["profile"].isin(corefams)]
        for file in os.listdir(dout_alis):
            if not os.path.splitext(file)[0] in corefams:
                os.remove(os.path.join(dout_alis, file))
    
    logging.info("applying score cutoffs to pangenome")
    genes = process_scores(hits, cutoffs)
    genes = pd.merge(genes_genomes, genes, how = "right")
    
    logging.info("writing output files")
    write_tsv(genes, fout_genes)
    write_tsv(cutoffs, fout_cutoffs)
    
    logging.info("removing temporary files and folders")
    shutil.rmtree(dout_ogseqs)
    os.remove(fout_hits)
    
def run_search(args): 
    
    fin_qpaths = args.qpaths
    din_db = args.db
    dout = args.outfolder
    threads = args.threads
    
    # define output paths/folders
    dout_tmp = os.path.join(dout, "tmp")
    din_alis = os.path.join(din_db, "alignments")
    fin_cutoffs = os.path.join(din_db, "cutoffs.csv")
    fout_hits = os.path.join(dout, "hits.tsv")
    fout_genes = os.path.join(dout, "genes.tsv")

    if os.path.isfile(fout_genes):
        logging.info("existing search results detected - moving on")
        return()
    
    # run profile search (function does its own logging)
    fins_alis = [os.path.join(din_alis, f) for f in os.listdir(din_alis)]
    fins_queries = read_fastapaths(fin_qpaths)
    run_profilesearch(fins_queries, fins_alis, fout_hits, dout_tmp, threads)
    
    logging.info("reading hits and score cutoffs")
    colnames = ["gene", "profile", "score"]
    hits = pd.read_csv(fout_hits, sep = "\t", names = colnames, 
        usecols = [0, 1, 2])
    hits[["gene", "profile"]] = hits[["gene", "profile"]].\
        applymap(lambda x: x.split(" ")[0])
    colnames = ["profile", "cutoff"]
    cutoffs = pd.read_csv(fin_cutoffs, sep = "\t", names = colnames)
    
    logging.info("applying score cutoffs to hits")
    genes = process_scores(hits, cutoffs)
    genes_genomes = extract_genes(fins_queries)
    genes = pd.merge(genes_genomes, genes, how = "right")
    
    logging.info("writing output files")
    write_tsv(genes, fout_genes)
    
    logging.info("removing temporary files and folders")
    os.remove(fout_hits)

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
  
    fin_faapaths = args.faapaths
    fin_coregenome = args.coregenome
    dout = args.outfolder
    core_filter = args.core_filter
    max_cores = args.max_cores
    fin_ffnpaths = args.ffnpaths

    sm_aas_fout = os.path.join(dout, "supermatrix_aas.fasta")
    sm_nucs_fout = os.path.join(dout, "supermatrix_nucs.fasta")
    seqs_aas_dio = os.path.join(dout, "seqs_aas")
    seqs_nucs_dio = os.path.join(dout, "seqs_nucs")
    alis_aas_dio = os.path.join(dout, "alis_aas")
    alis_nucs_dio = os.path.join(dout, "alis_nucs")
    
    # exit if requested supermatrices already exist
    if os.path.isfile(sm_aas_fout) and (fin_ffnpaths is None or \
        os.path.isfile(sm_nucs_fout)):
        logging.info("requested supermatrices already exist - moving on")
        return()
    
    logging.info("reading core genome") 
    coregenome = read_genes(fin_coregenome)
    orthogroups = coregenome["orthogroup"].unique()
    genomes = coregenome["genome"].unique()
    logging.info(f"detected {len(orthogroups)} orthogroups in "
        f"{len(genomes)} genomes")
    
    if core_filter != 0 or max_cores != 0:
        logging.info(f"applying core filter of {core_filter} and maximum "
            f"number of core genes of {max_cores}")
        corefams = determine_corefams(coregenome, core_filter, max_cores)
        coregenome = filter_groups(coregenome, corefams)
        orthogroups = coregenome["orthogroup"].unique()

    logging.info("removing same-genome copies of core genes")
    coregenome = coregenome.drop_duplicates(["genome", "orthogroup"], 
        keep = False)
    
    seqs_aas_fios = make_paths(orthogroups, seqs_aas_dio, ".fasta")
    alis_aas_fios = make_paths(orthogroups, alis_aas_dio, ".aln")
        
    # move on if dir already exists and is not empty
    if os.path.isdir(seqs_aas_dio) and os.listdir(seqs_aas_dio):
      
        logging.info("existing amino acid sequences detected - moving on")
        
    else:
    
        logging.info("gathering amino acid sequences of orthogroups")
        faa_fins = read_fastapaths(fin_faapaths)
        os.makedirs(seqs_aas_dio, exist_ok = True)
        gather_orthogroup_sequences(coregenome, faa_fins, seqs_aas_dio)
        
    if os.path.isdir(alis_aas_dio) and os.listdir(alis_aas_dio):
      
        logging.info("existing amino acid alignments detected - moving on")
        
    else: 

        logging.info("aligning orthogroups on the amino acid level")
        os.makedirs(alis_aas_dio, exist_ok = True)
        run_mafft_parallel(seqs_aas_fios, alis_aas_fios)

    logging.info("concatenating amino acid alignments")
    construct_supermatrix(coregenome, alis_aas_fios, sm_aas_fout)

    if not fin_ffnpaths is None:
      
        seqs_nucs_fios = make_paths(orthogroups, seqs_nucs_dio, ".fasta")
        alis_aas_fios = make_paths(orthogroups, alis_aas_dio, ".aln")
        alis_nucs_fios = make_paths(orthogroups, alis_nucs_dio, ".aln")
        
        if os.path.isdir(seqs_nucs_dio) and os.listdir(seqs_nucs_dio):
          
            logging.info("existing nucleotide sequences detected - moving on")
            
        else:
    
            logging.info("gathering nucleotide sequences of orthogroups")
            ffn_fins = read_fastapaths(fin_ffnpaths)
            os.makedirs(seqs_nucs_dio, exist_ok = True)
            gather_orthogroup_sequences(coregenome, ffn_fins, seqs_nucs_dio)
            
        if os.path.isdir(alis_nucs_dio) and os.listdir(alis_nucs_dio):
          
            logging.info("existing nucleotide alignments detected - moving on")
            
        else:

            logging.info("aligning orthogroups on the nucleotide level")
            os.makedirs(alis_nucs_dio, exist_ok = True)
            reverse_align_parallel(seqs_nucs_fios, alis_aas_fios, 
                alis_nucs_fios)

        logging.info("concatenating nucleotide alignments")
        construct_supermatrix(coregenome, alis_nucs_fios, sm_nucs_fout)
        
    logging.info("removing temporary folders")
    shutil.rmtree(seqs_aas_dio)
    shutil.rmtree(alis_aas_dio)
    try:
        shutil.rmtree(seqs_nucs_dio)
        shutil.rmtree(alis_nucs_dio)
    except FileNotFoundError:
        pass

def run_sample(args):
  
    if "exact" in args and args.exact:
        ali_mode = 3
        logging.info("identity calculation set to exact")
    else:
        ali_mode = 1
        logging.info("identity calculation set to approximate")

    logging.info("creating output subfolders")
    dio_seqs = os.path.join(args.outfolder, "tmp_seqs")
    dio_alis = os.path.join(args.outfolder, "tmp_alis")
    os.makedirs(dio_seqs, exist_ok = True)
    os.makedirs(dio_alis, exist_ok = True)
    
    fout_clusters = os.path.join(args.outfolder, "clusters.tsv")
    fout_seeds = os.path.join(args.outfolder, "seeds.txt")
    fout_identities = os.path.join(args.outfolder, "identities.tsv")
    
    if os.path.isfile(fout_clusters):
        logging.info("existing results detected - moving on")
        return()

    logging.info("reading core genome") 
    core = read_genes(args.coregenome)
    fams = core["orthogroup"].unique()
    genomes = core["genome"].unique()
    logging.info(f"detected {len(fams)} orthogroups in "
        f"{len(genomes)} genomes")
    
    if args.max_genomes == 0:
        args.max_genomes = len(genomes)
        logging.info(f"max_genomes set to {args.max_genomes}")
    
    logging.info("removing same-genome copies of core genes")
    core = core.drop_duplicates(["genome", "orthogroup"], keep = False)
    fams = core["orthogroup"].unique() # in case some fams were fully removed
    
    if os.path.isdir(dio_seqs) and os.listdir(dio_seqs):
      
        logging.info("existing sequences detected - moving on")
        
    else:
    
        logging.info("gathering sequences of orthogroups")
        fins_faas = read_fastapaths(args.fastapaths)
        gather_orthogroup_sequences(core, fins_faas, dio_seqs)
    
    logging.info("creating database for alignments")
    for dir in ["sequenceDB", "logs"]:
        makedirs_smart(f"{dio_alis}/{dir}")
    fio_seqs = f"{dio_alis}/seqs.fasta"
    open(fio_seqs, "w").close()
    for fam in fams:
        with open(f"{dio_seqs}/{fam}.fasta", "r") as hin_fam:
            with open(fio_seqs, "a+") as hout_seqs:
                hout_seqs.write(hin_fam.read())
    run_mmseqs(["createdb", fio_seqs, f"{dio_alis}/sequenceDB/db"], 
        f"{dio_alis}/logs/createdb.log", threads = args.threads)
    os.remove(fio_seqs)
    
    logging.info("initializing empty identity matrix")
    id_m = np.zeros([len(genomes), 0])
    
    logging.info("adding genomes until max_genomes or identity threshold is "
        "reached")
    core_grouped = core.copy().groupby("orthogroup")
    core_grouped = [core_fam for fam, core_fam in core_grouped]
    seeds = []
    while True:
      
        print("|", end = "", flush = True)
      
        # add column of zeros to identity matrix
        id_m = np.hstack((id_m, np.zeros([id_m.shape[0], 1])))
      
        # select seed
        max_ids = np.amax(id_m, 1) # max of each row = max along columns
        for seed in seeds: max_ids[seed] = 1.1 # avoid seeds being reused
        s = np.argmin(max_ids)
        seeds.append(s)
        min_max_ids = max_ids[s]
        seed_genome = genomes[s]
        
        # if sampled genomes close enough to each other: break the loop
        if min_max_ids > args.identity:
            print("")
            logging.info("identity threshold reached")
            # remove last column from identity matrix
            id_m = np.delete(id_m, -1, 1)
            seeds = seeds[:-1]
            break
        
        # prepare prefilter database
        for dir in ["prefDB", "alignmentDB"]:
            makedirs_smart(f"{dio_alis}/{dir}")
        pref = core_grouped.copy()
        for p, fam in enumerate(pref):
            if not seed_genome in fam.genome.tolist():
                pref[p] = None
                continue
            fam["query"] = fam[fam.genome == seed_genome].iloc[0]["gene"]
        pref = pd.concat(pref)
        pref = pref.rename(columns = {"gene": "target"})
        lookup = pd.read_csv(f"{dio_alis}/sequenceDB/db.lookup", sep = "\t", 
            usecols = [0, 1], names = ["id", "sequence"])
        lookup = dict(zip(lookup["sequence"].tolist(), lookup["id"].tolist()))
        pref["query_id"] = [lookup[q] for q in pref["query"].tolist()]
        pref["target_id"] = [lookup[t] for t in pref["target"].tolist()]
        pref = pref[["query_id", "target_id"]]
        write_tsv(pref, f"{dio_alis}/pref.tsv")
        run_mmseqs(["tsv2db", f"{dio_alis}/pref.tsv", f"{dio_alis}/prefDB/db", 
            "--output-dbtype", "7"], f"{dio_alis}/logs/tsv2db.log")
        
        # perform alignments 
        run_mmseqs(["align", f"{dio_alis}/sequenceDB/db",
            f"{dio_alis}/sequenceDB/db", f"{dio_alis}/prefDB/db", 
            f"{dio_alis}/alignmentDB/db", "--alignment-mode", str(ali_mode)], 
            f"{dio_alis}/logs/align.log", threads = args.threads)
        run_mmseqs(["createtsv", f"{dio_alis}/sequenceDB/db", 
            f"{dio_alis}/sequenceDB/db", f"{dio_alis}/alignmentDB/db", 
            f"{dio_alis}/hits.tsv"], f"{dio_alis}/logs/createtsv.log")
        
        # parse alignment results
        hits = pd.read_csv(f"{dio_alis}/hits.tsv", sep = "\t", 
            usecols = [1, 3], names = ["gene", "identity"])
        hits = core.merge(hits, on = "gene", how = "right")
        hits = hits.drop(["gene"], axis = 1)

        # calculate mean/median identity with seed per genome
        if args.method == "median":
            ids = hits.groupby("genome").aggregate(median)
        elif args.method[:4] == "mean":
            p = args.method[4:]
            if (p == ""): p = "100"
            p = int(p) / 100
            def meanp(l, p):
                start = int(round(((1 - p) / 2) * len(l)))
                stop = int(round((1 - (1 - p) / 2) * len(l)))
                if start == stop: return(median(l))
                return(mean(sorted(l)[start:stop]))
            ids = hits.groupby("genome").aggregate(lambda l: meanp(l, p))
                
        # add mean/median identity to identity matrix
        for g, genome in enumerate(genomes):
            id_m[g, -1] = ids.loc[genome, "identity"]
            
        # if enough sampled genomes: break the loop
        if np.size(id_m, 1) == args.max_genomes:
            print("")
            logging.info("max number of genomes sampled")
            break
            
    logging.info("writing seeds.txt")
    with open(fout_seeds, "a") as hout_seeds:
        for s in seeds: hout_seeds.write(genomes[s] + "\n")
            
    logging.info("writing identities.tsv")
    id_df = pd.DataFrame(id_m)
    id_df.columns = [genomes[s] for s in seeds]
    id_df["genome"] = genomes 
    cols = id_df.columns.tolist()
    cols = [cols[-1]] + cols[:-1]
    id_df = id_df[cols]
    id_df.to_csv(fout_identities, sep = "\t", index = False, header = True)
        
    logging.info("clustering the genomes")
    clusters = np.argmax(id_m, 1)
    genomes_clusters = pd.DataFrame({"genome": genomes, "cluster": clusters})
    logging.info(f"{len(set(clusters))} clusters found")
    
    logging.info("writing clusters.tsv")
    write_tsv(genomes_clusters, fout_clusters)
        
    logging.info("removing temporary folders")
    shutil.rmtree(dio_seqs)
    shutil.rmtree(dio_alis)
    
def run_fetch(args):

    logging.info("creating subfolder for fastas")
    dout_seqs = os.path.join(args.outfolder, "fastas")
    os.makedirs(dout_seqs, exist_ok = True)

    logging.info("reading genes") 
    genes = read_genes(args.genes)
    fams = genes["orthogroup"].unique()
    genomes = genes["genome"].unique()
    logging.info(f"detected {len(fams)} orthogroups in "
        f"{len(genomes)} genomes")
        
    logging.info("gathering sequences of orthogroups")
    fins_fastas = read_fastapaths(args.fastapaths)
    gather_orthogroup_sequences(genes, fins_fastas, dout_seqs)
    
def run_core(args):
  
    fin_faapaths = args.faapaths
    dout = args.outfolder
    method = args.method
    seeds = args.seeds
    core_prefilter = args.core_prefilter
    core_filter = args.core_filter
    max_cores = args.max_cores
    threads = args.threads
    
    # define paths
    dout_seedpan = os.path.join(dout, "seedpan")
    dout_seedcore = os.path.join(dout, "seedcore")
    fout_seedpaths = os.path.join(dout, "seedpaths.txt")
    fout_nonseedpaths = os.path.join(dout, "nonseedpaths.txt")
    fout_seedpan = os.path.join(dout_seedpan, "pangenome.tsv")
    fout_genes_core = os.path.join(dout_seedcore, "genes.tsv")
    fout_genes = os.path.join(dout, "genes.tsv")
    
    # make output subfolders 
    for dir in [dout_seedpan, dout_seedcore]:
        os.makedirs(dir, exist_ok = True)

    logging.info("selecting random seed genomes")
    fins_faas = read_fastapaths(fin_faapaths)
    if not os.path.isfile(fout_seedpaths):
        fins_seeds = sample(fins_faas, seeds)
        fins_nonseeds = [fin for fin in fins_faas if not fin in fins_seeds]
        write_tsv(pd.DataFrame({"path": fins_seeds}), fout_seedpaths)
        write_tsv(pd.DataFrame({"path": fins_nonseeds}), fout_nonseedpaths)

    logging.info("STEP 1 - inferring pangenome of seed genomes")
    args_pan = Namespace(faapaths = fout_seedpaths, outfolder = dout_seedpan,
        method = method, threads = threads)
    run_pan(args_pan)
    
    logging.info("STEP 2 - building database of seed core genes and searching "
        "in seed faas")
    args_build = Namespace(faapaths = fout_seedpaths, pangenome = fout_seedpan,
        outfolder = dout_seedcore, core_prefilter = core_prefilter,
        core_filter = core_filter, max_cores = max_cores, threads = threads)
    run_build(args_build)
    
    if os.stat(fout_nonseedpaths).st_size == 0:
        logging.info("all faa files are seeds - skipping search in non-seeds")
        shutil.copy(fout_genes_core, fout_genes)
        return()

    logging.info("STEP 3 - identifying core orthogroups in non-seed genomes")
    args_search = Namespace(qpaths = fout_nonseedpaths, db = dout_seedcore,
        outfolder = dout, threads = threads)
    run_search(args_search)
    
    logging.info("adding search results in seeds to search results of "
        "non-seeds")
    with open(fout_genes, "a") as hout_genes:
        with open(fout_genes_core) as hout_genes_core:
            for line in hout_genes_core:
                hout_genes.write(line)
    
