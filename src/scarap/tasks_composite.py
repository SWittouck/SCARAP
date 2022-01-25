import logging
import os
import shutil

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
