import logging
import os

from argparse import Namespace
from random import sample

from utils import *
from tasks_core import *

def run_pan_pipeline(args):

    logging.info("STEP 1 - inferring pangenome")
    run_pan(args)

    logging.info("STEP 2 - building profile hmm database of pangenome")
    args.pangenome = os.path.join(args.outfolder, "pangenome.tsv")
    args.min_genomes = 1 # default not yet set in all situations
    if not args.species is None:
        args.faapaths = os.path.join(args.outfolder, "reprpaths.txt")
        args.pangenome = os.path.join(args.outfolder, "metapangenome",
            "pangenome.tsv")
    run_build(args)

    logging.info("STEP 3 - training score cutoffs for profiles")
    args_search = Namespace(qpaths = args.faapaths, db = args.outfolder,
        outfolder = args.outfolder, trainstrategy = "pan", 
        threads = args.threads, pangenome = args.pangenome)
    run_search(args_search)

def run_core_pipeline(args):

    logging.info("selecting random seed genomes")
    seedsdio = os.path.join(args.outfolder, "seeds")
    os.makedirs(seedsdio, exist_ok = True)
    faapaths = read_fastapaths(args.faapaths)
    seedpathsfio = os.path.join(seedsdio, "seedpaths.txt")
    if not os.path.isfile(os.path.join(seedpathsfio)):
        seedpaths = sample(faapaths, args.seeds)
        write_tsv(pd.DataFrame({"path": seedpaths}), seedpathsfio)
    else:
        seedpaths = read_lines(seedpathsfio)

    logging.info("STEP 1 - inferring pangenome of seed genomes")
    args_pan = Namespace(faapaths = seedpathsfio, outfolder = seedsdio,
        method = args.method, threads = args.threads)
    run_pan(args_pan)
    
    logging.info("STEP 2a - selecting candidate core genes")
    seedspanfio = os.path.join(seedsdio, "pangenome.tsv")
    seedspan = read_genes(seedspanfio)
    seedspanfams = checkgroups(seedspan)
    tokeep = seedspanfams.occurrence_singlecopy >= args.seedfilter / args.seeds
    seedscorefams = seedspanfams[tokeep]
    logging.info(f"{len(seedscorefams.index)} candidate core genes found")
    if (args.max_genes > 0):
        seedscorefams = seedscorefams\
            .sort_values("occurrence_singlecopy", ascending = False)\
            .head(args.max_genes)
    seedscore = filter_groups(seedspan, seedscorefams.orthogroup)
    logging.info(f"{len(seedscorefams.index)} candidate core genes selected")
    seedscorefio = os.path.join(seedsdio, "coregenome.tsv")
    write_tsv(seedscore, seedscorefio)

    logging.info("STEP 2b - building database of candidate core genes")
    candsdio = os.path.join(args.outfolder, "cands")
    os.makedirs(candsdio, exist_ok = True)
    args_build = Namespace(faapaths = seedpathsfio, pangenome = seedscorefio,
        outfolder = candsdio, min_genomes = 0)
    run_build(args_build)

    logging.info("STEP 3 - searching candidate core orthogroups in all genomes")
    args_search = Namespace(qpaths = args.faapaths, db = candsdio,
        outfolder = candsdio, trainstrategy = "core", threads = args.threads)
    run_search(args_search)

    logging.info("STEP 4 - selecting core orthogroups from candidates")
    candcoregenome = read_genes(os.path.join(candsdio, "genes.tsv"))
    candcorefams = checkgroups(candcoregenome)
    tokeep = candcorefams.occurrence_singlecopy >= args.allfilter
    corefams = candcorefams[tokeep].orthogroup
    coregenome = filter_groups(candcoregenome, corefams)
    write_tsv(corefams, os.path.join(args.outfolder, "core_orthogroups.txt"))
    write_tsv(coregenome, os.path.join(args.outfolder, "coregenome.tsv"))
