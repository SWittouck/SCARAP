import logging
import numpy as np
import os
import pandas as pd
import shutil

from argparse import Namespace
from concurrent.futures import ProcessPoolExecutor
from statistics import median, mean
from random import sample

from scarap.callers import *
from scarap.computers import *
from scarap.helpers import * 
from scarap.pan import *
from scarap.readerswriters import *
from scarap.utils import *

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
    
    faafins = read_fastapaths(args.faa_files)
    
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
  
    logging.info("creating genome table with species and faa paths")
    genometbl1 = read_species(args.species)
    genometbl1.species = [sp.replace(" ", "_") for sp in genometbl1.species]
    faapaths = read_fastapaths(args.faa_files)
    genomes = [filename_from_path(path) for path in faapaths]
    genometbl2 = pd.DataFrame({"faapath": faapaths, "genome": genomes})
    genometbl = pd.merge(genometbl1, genometbl2)
    
    logging.info("creating the necessary subfolders")
    speciespansdio = os.path.join(args.outfolder, "speciespangenomes")
    pseudogenomesdio = os.path.join(args.outfolder, "pseudogenomes")
    pseudogenomesfio = os.path.join(args.outfolder, "pseudogenomes.tsv")
    pseudopandio = os.path.join(args.outfolder, "pseudopangenome")
    pseudopanfio = os.path.join(args.outfolder, "pseudopangenome.tsv")
    pangenomefout = os.path.join(args.outfolder, "pangenome.tsv")
    os.makedirs(speciespansdio, exist_ok = True)
    os.makedirs(pseudogenomesdio, exist_ok = True)
    os.makedirs(pseudopandio, exist_ok = True)
    
    logging.info("PHASE 1: inferring species-level pangenomes")
    for species, genomesubtbl in genometbl.groupby("species"):
        speciespanfio = os.path.join(speciespansdio, species + ".tsv")
        if os.path.exists(speciespanfio):
            logging.info(f"existing pangenome found for {species}")
            continue
        logging.info(f"started inferring pangenome of {species}")
        dout = os.path.join(speciespansdio, species)
        os.makedirs(dout, exist_ok = True)
        faapaths_sub = genomesubtbl["faapath"].values.tolist()
        faapathsfio = os.path.join(dout, "faapaths.txt")
        write_lines(faapaths_sub, faapathsfio)
        run_pan_nonhier(Namespace(faa_files = faapathsfio, outfolder = dout,
            threads = args.threads, method = args.method))
        shutil.move(os.path.join(dout, "pangenome.tsv"), speciespanfio)
        shutil.rmtree(dout)
    
    logging.info("PHASE 2: constructing pseudogenome per species")
    n_species = genometbl["species"].nunique()
    pseudogenomes = [None] * n_species
    for ix, speciespanfio in enumerate(listpaths(speciespansdio)):
        species = filename_from_path(speciespanfio)
        logging.info(f"constructing pseudogenome for {species}")
        speciespan = read_genes(speciespanfio)
        tmpdio = os.path.join(args.outfolder, "temp")
        pseudogenomes[ix] = create_pseudogenome(speciespan, faapaths, tmpdio)
        pseudogenomes[ix]["species"] = species
    pseudogenomes = pd.concat(pseudogenomes)
    write_tsv(pseudogenomes, pseudogenomesfio)
    
    logging.info("PHASE 3: inferring pseudopangenome")
    pseudogenomes = pseudogenomes.rename(columns = {"species": "orthogroup"})
    gather_orthogroup_sequences(pseudogenomes, faapaths, pseudogenomesdio)
    run_pan_nonhier(Namespace(faa_files = pseudogenomesdio, 
        outfolder = pseudopandio, threads = args.threads, method = args.method))
    shutil.move(os.path.join(pseudopandio, "pangenome.tsv"), pseudopanfio)
    shutil.rmtree(pseudogenomesdio)
    shutil.rmtree(pseudopandio)
    
    logging.info("PHASE 4: inflating pseudopangenome with species pangenomes")
    pangenomes = [None] * len(os.listdir(speciespansdio))
    for ix, speciespanfio in enumerate(listpaths(speciespansdio)):
        pangenomes[ix] = read_genes(speciespanfio)
        pangenomes[ix]["species"] = filename_from_path(speciespanfio)
    pangenomes = pd.concat(pangenomes)
    pangenomes = pangenomes.rename(columns = {"orthogroup": "speciesfam"})
    pseudopan = read_genes(pseudopanfio)
    pseudopan = pseudopan.rename(columns = {"genome": "species"})
    speciesfamtbl = pd.merge(pseudopan, pangenomes, on = ["gene", "species"], 
        how = "left")
    speciesfamtbl = speciesfamtbl[["speciesfam", "species", "orthogroup"]]
    pangenome = pd.merge(pangenomes, speciesfamtbl, 
        on = ["speciesfam", "species"])
    pangenome = pangenome[["gene", "genome", "orthogroup"]]
    write_tsv(pangenome, pangenomefout)

def run_build(args): 
    
    fin_faapaths = args.faa_files
    fin_pangenome = args.pangenome
    dout = args.outfolder
    core_prefilter = args.core_prefilter
    core_filter = args.core_filter
    max_cores = args.max_core_genes
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
        corefams = determine_corefams(pangenome, core_prefilter)
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
        logging.info(f"{len(corefams)} core genes were identified")
    
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
    
    fin_qpaths = args.faa_files
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
    if os.path.isfile(fout_hits):
        logging.info("existing mmseqs2 hits detected - moving on")
    else:
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
    logging.info("extracting genome names from faa files")
    genes_genomes = extract_genes(fins_queries, threads = threads)
    logging.info("merging search results and genome names")
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

def run_concat(args):
  
    fin_faapaths = args.faa_files
    fin_pangenome = args.pangenome
    dout = args.outfolder
    core_filter = args.core_filter
    max_cores = args.max_core_genes
    fin_ffnpaths = args.ffn_files

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
    coregenome = read_genes(fin_pangenome)
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
    core = read_genes(args.pangenome)
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
        fins_faas = read_fastapaths(args.fasta_files)
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
        hits = hits.drop(["gene", "orthogroup"], axis = 1)

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
    genes = read_genes(args.pangenome)
    fams = genes["orthogroup"].unique()
    genomes = genes["genome"].unique()
    logging.info(f"detected {len(fams)} orthogroups in "
        f"{len(genomes)} genomes")
        
    logging.info("gathering sequences of orthogroups")
    fins_fastas = read_fastapaths(args.fasta_files)
    gather_orthogroup_sequences(genes, fins_fastas, dout_seqs)
    
def run_core(args):
  
    fin_faapaths = args.faa_files
    dout = args.outfolder
    method = args.method
    seeds = args.seeds
    core_prefilter = args.core_prefilter
    core_filter = args.core_filter
    max_cores = args.max_core_genes
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
    args_pan = Namespace(faa_files = fout_seedpaths, outfolder = dout_seedpan,
        method = method, threads = threads)
    run_pan(args_pan)
    
    logging.info("STEP 2 - building database of seed core genes and searching "
        "in seed faas")
    args_build = Namespace(faa_files = fout_seedpaths, pangenome = fout_seedpan,
        outfolder = dout_seedcore, core_prefilter = core_prefilter,
        core_filter = core_filter, max_core_genes = max_cores, 
        threads = threads)
    run_build(args_build)
    
    if os.stat(fout_nonseedpaths).st_size == 0:
        logging.info("all faa files are seeds - skipping search in non-seeds")
        shutil.copy(fout_genes_core, fout_genes)
        return()

    logging.info("STEP 3 - identifying core orthogroups in non-seed genomes")
    args_search = Namespace(faa_files = fout_nonseedpaths, db = dout_seedcore,
        outfolder = dout, threads = threads)
    run_search(args_search)
    
    logging.info("adding search results in seeds to search results of "
        "non-seeds")
    with open(fout_genes, "a") as hout_genes:
        with open(fout_genes_core) as hout_genes_core:
            for line in hout_genes_core:
                hout_genes.write(line)
    
