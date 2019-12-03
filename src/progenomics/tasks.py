# This script contains a function for each task, called run_<task>. These
# functions check input/output files/folders and dependencies. For most of them,
# the main functionality is outsourced to functions of the type run_<task>_nb
# (nb stands for "no bullshit").

import logging
import os
import shutil
import subprocess
import sys

from argparse import Namespace
from checkers import *
from random import sample
from readerswriters import *
from runners import *
from utils import *

def run_pan_nonhier(args):

    pangenomefout = os.path.join(args.outfolder, "pangenome.tsv")
    if os.path.isfile(pangenomefout):
        logging.info("existing pangenome detected - moving on")
        return()

    logging.info("creating orthofinder directory")
    dir_orthofinder = os.path.join(args.outfolder, "orthofinder")
    os.makedirs(dir_orthofinder)

    logging.info("running orthofinder")
    faafins = read_lines(args.faapaths)
    logfile = os.path.join(args.outfolder, "orthofinder.log")
    engine = args.method.split("_")[1]
    run_orthofinder(faafins, dir_orthofinder, logfile, args.threads, engine)

    logging.info("creating tidy pangenome file")
    pangenome = read_pangenome_orthofinder(dir_orthofinder)
    write_tsv(pangenome, pangenomefout)

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
        collapse_pangenome(speciespanfio, faapathsfio, reprfio, species,
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

def run_pan_nb(args):
    if "species" in args and not args.species is None:
        run_pan_hier(args)
    else:
        run_pan_nonhier(args)

def run_pan(args):

    logging.info("welcome to the pan task")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_fastas(args.faapaths)
    if not args.species is None:
        check_infile(args.species)

    logging.info("checking dependencies")
    check_tool("orthofinder")

    run_pan_nb(args)

def run_build_nb(args):

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

def run_build(args):

    logging.info("welcome to the build task")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_infile(args.pangenome)

    logging.info("checking dependencies")
    check_tool("hmmbuild", ["-h"])
    check_tool("mafft", ["--help"]) # --help avoids mafft interactive mode

    run_build_nb(args)

def run_search_nb(args):

    genesfout = os.path.join(args.outfolder, "genes.tsv")
    if os.path.isfile(genesfout):
        logging.info("existing search results detected - moving on")
        return()

    cutoffsfio = os.path.join(args.db, "orthogroups.tsv")

    logging.info("performing hmmsearch")
    queryfins = read_lines(args.qpaths)
    domtblfio = os.path.join(args.outfolder, "hmmer_domtbl.tmp")
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

def run_search(args):

    logging.info("welcome to the search task")

    logging.info("checking arguments other than output folder")
    check_infile(args.qpaths)
    check_fastas(args.qpaths)
    check_db(args.db)
    if args.trainstrategy is None:
        check_infile(os.path.join(args.db, "orthogroups.tsv"))
    else:
        check_outfile(os.path.join(args.db, "orthogroups.tsv"))
    if args.trainstrategy == "pan":
        if args.pangenome is None:
            logging.error("pangenome is required for cutoff training with the "
            " pan strategy")
            sys.exit(1)
        check_infile(args.pangenome)
    if not args.orthogroups is None:
        check_infile(args.orthogroups)
        logging.warning("the option to subset orthogroups is not yet "
            "implemented")

    logging.info("checking dependencies")
    check_tool("hmmsearch", ["-h"])

    run_search_nb(args)

def run_checkgenomes_nb(args):

    logging.info("checking genomes")
    coregenome = read_genes(args.coregenome)
    genomes = checkgenomes(coregenome)
    write_tsv(genomes, os.path.join(args.outfolder, "genomes.tsv"))

def run_checkgenomes(args):

    logging.info("welcome to the checkgenomes task")

    logging.info("checking arguments other than output folder")
    check_infile(args.coregenome)

    run_checkgenomes_nb(args)

def run_checkgroups_nb(args):

    logging.info("checking core orthogroups")
    coregenome = read_genes(args.coregenome)
    orthogroups = checkgroups(coregenome)
    write_tsv(orthogroups, os.path.join(args.outfolder, "orthogroups.tsv"))

def run_checkgroups(args):

    logging.info("welcome to the checkgroups task")

    logging.info("checking arguments other than output folder")
    check_infile(args.coregenome)

    run_checkgroups_nb(args)

def run_filter_nb(args):

    logging.info("reading pangenome")
    pangenome = read_genes(args.pangenome)
    if not args.genomes is None:
        logging.info("filtering genomes")
        genomes = read_lines(args.genomes)
        pangenome = filter_genomes(pangenome, genomes)
    if not args.orthogroups is None:
        logging.info("filtering orthogroups")
        orthogroups = read_lines(args.orthogroups, orthogroups)
        pangenome = filter_groups(pangenome, orthogroups)
    write_tsv(pangenome, os.path.join(args.outfolder, "pangenome.tsv"))

def run_filter(args):

    logging.info("welcome to the filter task")

    logging.info("checking arguments other than output folder")
    check_infile(args.pangenome)
    if not args.genomes is None:
        check_infile(args.genomes)
    if not args.orthogroups is None:
        check_infile(args.orthogroups)

    run_filter_nb(args)

def run_supermatrix(args):

    logging.info("welcome to the supermatrix task")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_infile(args.coregenome)
    if not args.ffnpaths is None:
        check_infile(args.ffnpaths)

    logging.info("checking dependencies")
    check_tool("mafft", ["--help"]) # --help avoids mafft interactive mode

    supermatrixfout = os.path.join(args.outfolder, "supermatrix.fasta")
    if os.path.isfile(supermatrixfout):
        logging.info("existing supermatrix detected - moving on")
        return()

    logging.info("creating output subfolders")
    orthogroupsdio = os.path.join(args.outfolder, "orthogroups")
    alignmentsdio = os.path.join(args.outfolder, "alignments")
    os.makedirs(orthogroupsdio, exist_ok = True)
    os.makedirs(alignmentsdio, exist_ok = True)

    logging.info("gathering sequences of orthogroups")
    coregenome = read_genes(args.coregenome)
    faafins = read_lines(args.faapaths)
    gather_orthogroup_sequences(coregenome, faafins, orthogroupsdio)
    logging.info(f"gathered sequences for {len(os.listdir(orthogroupsdio))} "
        f"core orthogroups")

    logging.info("aligning orthogroups")
    orthogroups = [os.path.splitext(file)[0] for file in
        os.listdir(orthogroupsdio)]
    orthogroupfouts = make_paths(orthogroups, orthogroupsdio, ".fasta")
    alifouts = make_paths(orthogroups, alignmentsdio, ".aln")
    run_mafft_parallel(orthogroupfouts, alifouts)

    logging.info("concatenating alignments")
    construct_supermatrix(coregenome, alifouts, supermatrixfout)

def run_pan_pipeline(args):

    logging.info("welcome to the pan pipeline")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_fastas(args.faapaths)

    logging.info("checking dependencies")
    check_tool("orthofinder")
    check_tool("hmmbuild", ["-h"])
    check_tool("mafft", ["--help"])
    check_tool("hmmsearch", ["-h"])

    logging.info("STEP 1 - inferring pangenome")
    run_pan_nb(args)

    logging.info("STEP 2 - building profile hmm database of pangenome")
    args.pangenome = os.path.join(args.outfolder, "pangenome.tsv")
    args.min_genomes = 1 # default not yet set in all situations
    if not args.species is None:
        args.faapaths = os.path.join(args.outfolder, "reprpaths.txt")
        args.pangenome = os.path.join(args.outfolder, "metapangenome",
            "pangenome.tsv")
    run_build_nb(args)

    logging.info("STEP 3 - training score cutoffs for profiles")
    args_search = Namespace(qpaths = args.faapaths, db = args.outfolder,
        outfolder = args.outfolder, trainstrategy = "pan",
        pangenome = args.pangenome)
    run_search_nb(args_search)

def run_core_pipeline(args):

    logging.info("welcome to the core pipeline")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_fastas(args.faapaths)

    logging.info("checking dependencies")
    check_tool("orthofinder")
    check_tool("hmmbuild", ["-h"])
    check_tool("mafft", ["--help"])
    check_tool("hmmsearch", ["-h"])

    logging.info("selecting random seed genomes")
    seedsdio = os.path.join(args.outfolder, "seeds")
    os.makedirs(seedsdio, exist_ok = True)
    faapaths = read_lines(args.faapaths)
    seedpathsfio = os.path.join(seedsdio, "seedpaths.txt")
    if not os.path.isfile(os.path.join(seedpathsfio)):
        seedpaths = sample(faapaths, args.seeds)
        write_tsv(pd.DataFrame({"path": seedpaths}), seedpathsfio)
    else:
        seedpaths = read_lines(seedpathsfio)

    logging.info("STEP 1 - inferring pangenome of seed genomes")
    args_pan = Namespace(faapaths = seedpathsfio, outfolder = seedsdio,
        method = args.method, threads = args.threads)
    run_pan_nb(args_pan)

    logging.info("STEP 2 - building database of candidate core genes")
    candsdio = os.path.join(args.outfolder, "cands")
    os.makedirs(candsdio, exist_ok = True)
    seedspanfio = os.path.join(seedsdio, "pangenome.tsv")
    args_build = Namespace(faapaths = seedpathsfio, pangenome = seedspanfio,
        outfolder = candsdio, min_genomes = args.seedfilter)
    run_build_nb(args_build)

    logging.info("STEP 3 - searching candidate core orthogroups in all genomes")
    args_search = Namespace(qpaths = args.faapaths, db = candsdio,
        outfolder = candsdio, trainstrategy = "core")
    run_search_nb(args_search)

    logging.info("STEP 4 - selecting core orthogroups from candidates")
    candcoregenome = read_genes(os.path.join(candsdio, "genes.tsv"))
    candcorefams = checkgroups(candcoregenome)
    corefams = candcorefams[candcorefams.coreness >= args.allfilter].orthogroup
    coregenome = filter_groups(candcoregenome, corefams)
    write_tsv(corefams, os.path.join(args.outfolder, "core_orthogroups.txt"))
    write_tsv(coregenome, os.path.join(args.outfolder, "coregenome.tsv"))
