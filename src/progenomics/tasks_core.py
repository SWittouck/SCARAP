import logging
import numpy as np
import os
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from statistics import median

from argparse import Namespace
from utils import *
from readerswriters import *
from computers import *
from callers import *
from pan import *

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
        run_hmmsearch(args.db, queryfins, domtblfio, args.threads)
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

        logging.info("reading core genome") 
        coregenome = read_genes(args.coregenome)
        orthogroups = coregenome["orthogroup"].unique()
        genomes = coregenome["genome"].unique()
        logging.info(f"detected {len(orthogroups)} orthogroups in "
            f"{len(genomes)} genomes")
    
        logging.info("removing same-genome copies of core genes")
        coregenome = coregenome.drop_duplicates(["genome", "orthogroup"], 
            keep = False)

        logging.info("gathering amino acid sequences of orthogroups")
        faa_fins = read_lines(args.faapaths)
        gather_orthogroup_sequences(coregenome, faa_fins, seqs_aas_dio)
        # logging.info(f"gathered sequences for {len(os.listdir(seqs_aas_dio))} "
        #     f"core orthogroups")

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
        
    logging.info("removing temporary folders")
    shutil.rmtree(seqs_aas_dio)
    shutil.rmtree(alis_aas_dio)
    try:
        shutil.rmtree(seqs_nucs_dio)
        shutil.rmtree(alis_nucs_dio)
    except FileNotFoundError:
        pass

def run_clust(args):
  
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
    if os.path.isfile(fout_clusters):
        logging.info("existing cluster file detected - moving on")
        return()

    logging.info("reading core genome") 
    core = read_genes(args.coregenome)
    fams = core["orthogroup"].unique()
    genomes = core["genome"].unique()
    logging.info(f"detected {len(fams)} orthogroups in "
        f"{len(genomes)} genomes")
    
    if args.max_clusters == 0:
        args.max_clusters = len(genomes)
        logging.info(f"set max_clusters to {args.max_clusters}")
    
    logging.info("removing same-genome copies of core genes")
    core = core.drop_duplicates(["genome", "orthogroup"], keep = False)
    fams = core["orthogroup"].unique() # in case some fams were fully removed
    
    logging.info("gathering sequences of orthogroups")
    fins_faas = read_lines(args.fastapaths)
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
    
    logging.info("adding cluster seeds until max_clusters or identity "
        "threshold is reached")
    core_grouped = core.copy().groupby("orthogroup")
    core_grouped = [core_fam for fam, core_fam in core_grouped]
    while True:
      
        print("|", end = "", flush = True)
      
        # add column of zeros to identity matrix
        id_m = np.hstack((id_m, np.zeros([id_m.shape[0], 1])))
      
        # select seed
        max_ids = np.amax(id_m, 1) # max of each row = max along columns
        s = np.argmin(max_ids)
        min_max_ids = max_ids[s]
        seed_genome = genomes[s]
        
        # if clusters small enough: break the loop
        if min_max_ids > args.identity:
            print("")
            logging.info("identity threshold reached")
            id_m = np.delete(id_m, -1, 1)
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
        
        # calculate median identity per genome and add to identity matrix
        ids = hits.groupby("genome").aggregate(lambda l: median(l))
        for g, genome in enumerate(genomes):
            id_m[g, -1] = ids.loc[genome, "identity"]
            
        # if enough clusters: break the loop
        if np.size(id_m, 1) == args.max_clusters:
            print("")
            logging.info("max number of clusters reached")
            break
        
    logging.info("assigning cluster numbers")
    clusters = np.argmax(id_m, 1)
    genomes_clusters = pd.DataFrame({"genome": genomes, "cluster": clusters})
    logging.info(f"{len(set(clusters))} clusters found")
    
    logging.info("writing clusters.tsv")
    write_tsv(genomes_clusters, fout_clusters)
        
    logging.info("removing temporary folders")
    shutil.rmtree(dio_seqs)
    shutil.rmtree(dio_alis)
