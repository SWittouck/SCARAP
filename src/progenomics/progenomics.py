#! /usr/bin/env python3

# This is the main script of progenomics; it only contains the commandline
# interface.

__author__ = "Stijn Wittouck"
__version__ = "0.2.0"

import argparse
import logging
import sys

from utils import *
from taskwrappers import *

def print_help():

    message = '''\
VERSION
    {0}
AUTHORS
    Stijn Wittouck (development)
    Sarah Lebeer (supervision)
USAGE
    progenomics [-h] <task> <task-specific arguments>
TASKS
    pan           --> infer a pangenome from a set of faa files
    build         --> build a profile HMM database for a core/pangenome
    search        --> search query genes in a core/pangenome database
    checkgenomes  --> assess the quality of genomes in a core genome
    checkgroups   --> assess the quality of orthogroups in a core genome
    filter        --> filter the genomes/orthogroups in a pangenome
    supermatrix   --> construct a concatenated core orthogroup alignment from a
                      core genome
PIPELINES
    pan-pipeline  --> infer a pangenome, build a profile HMM database and train
                      score cutoffs from a set of faa files
    core-pipeline --> infer a core genome, build a profile HMM database and
                      train score cutoffs from a set of faa files
DOCUMENTATION
    https://github.com/swittouck/progenomics\
'''

    print(message.format(__version__))

def print_intro():

    message = '''\

This is progenomics version {0}
'''

    print(message.format(__version__))

def parse_arguments():

    parser = argparse.ArgumentParser(
        add_help = False
    )

    parser.add_argument("-h", "--help", action = 'store_true')

    subparsers = parser.add_subparsers()

    # parent parser for the pan and pan-pipeline tasks
    parser_pan_parent = argparse.ArgumentParser(add_help = False)
    parser_pan_parent.add_argument("faapaths",
        help = "input file with paths to faa files of genomes")
    parser_pan_parent.add_argument("-d", "--method", default = "TRE-F",
        choices = ["TRE", "TRE-F", "TRE-FS", "CLU", "CLU-F", "PRO", "ORT-B", 
        "ORT-D"],
        help = "pangenome inference method [default: TRE-F]")
    parser_pan_parent.add_argument("outfolder",
        help = "output folder for pangenome file")
    parser_pan_parent.add_argument("-s", "--species",
        help = "input file with species of genomes; if given, a "
            "hierarchical pangenome strategy will be used")
    parser_pan_parent.add_argument("-t", "--threads", default = 8, type = int,
        help = "number of threads to use [default 8]")
    parser_pan_parent.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_pan = subparsers.add_parser('pan', parents = [parser_pan_parent])
    parser_pan.set_defaults(func = run_pan_withchecks)

    parser_build = subparsers.add_parser('build')
    parser_build.set_defaults(func = run_build_withchecks)
    parser_build.add_argument("faapaths",
        help = "input file with paths to faa files of genomes")
    parser_build.add_argument("pangenome",
        help = "input file with pangenome")
    parser_build.add_argument("outfolder",
        help = "output folder for hmm database")
    parser_build.add_argument("-m", "--min_genomes", default = 1, type = int,
        help = "minimum number of genomes an orthogroup should be present in "
            "to be included in the db [default 1]")
    parser_build.add_argument("-t", "--threads", default = 8, type = int,
        help = "number of threads to use [default 8]")
    parser_build.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_search = subparsers.add_parser('search')
    parser_search.set_defaults(func = run_search_withchecks)
    parser_search.add_argument("qpaths",
        help = "input file with paths to files with amino acid sequences of "
            "query genes")
    parser_search.add_argument("db",
        help = "input folder with hmm database")
    parser_search.add_argument("outfolder",
        help = "output folder for search hits")
    parser_search.add_argument("-y", "--trainstrategy",
        choices = ['pan', 'core'],
        help = "training strategy to use; if given, orthogroup score cutoff "
            "values will be determined [options: core, pan]")
    parser_search.add_argument("-p", "--pangenome",
        help = "input file with pangenome; required for cutoff training with "
            "the 'pan' strategy")
    parser_search.add_argument("-o", "--orthogroups",
        help = "input file with subset of orthogroups in the database to use")
    parser_search.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_checkgenomes = subparsers.add_parser("checkgenomes")
    parser_checkgenomes.set_defaults(func = run_checkgenomes_withchecks)
    parser_checkgenomes.add_argument("coregenome",
        help = "input file with a core genome")
    parser_checkgenomes.add_argument("outfolder",
        help = "output folder for genome statistics")
    parser_checkgenomes.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_checkgroups = subparsers.add_parser("checkgroups")
    parser_checkgroups.set_defaults(func = run_checkgroups_withchecks)
    parser_checkgroups.add_argument("coregenome",
        help = "input file with a core genome")
    parser_checkgroups.add_argument("outfolder",
        help = "output folder for orthogroup statistics")
    parser_checkgroups.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_filter = subparsers.add_parser('filter')
    parser_filter.set_defaults(func = run_filter_withchecks)
    parser_filter.add_argument("pangenome",
        help = "input file with pangenome")
    parser_filter.add_argument("outfolder",
        help = "output file for filtered pangenome")
    parser_filter.add_argument("-g", "--genomes",
        help = "input file with genomes to extract from pangenome")
    parser_filter.add_argument("-o", "--orthogroups",
        help = "input file with orthogroups to extract from pangenome")
    parser_filter.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_supermatrix = subparsers.add_parser('supermatrix')
    parser_supermatrix.set_defaults(func = run_supermatrix_withchecks)
    parser_supermatrix.add_argument("faapaths",
        help = "input file with paths to faa files of genomes")
    parser_supermatrix.add_argument("coregenome",
        help = "input file with core genome")
    parser_supermatrix.add_argument("outfolder",
        help = "output file for supermatrix fasta file based on core genome")
    parser_supermatrix.add_argument("-n", "--ffnpaths",
        help = "input file with paths to ffn files of genomes; if given, "
            "a nucleotide supermatrix will be constructed in addition to the "
            "amino acid suprmatrix")
    parser_supermatrix.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_pan_pipeline = subparsers.add_parser('pan-pipeline',
        parents = [parser_pan_parent])
    parser_pan_pipeline.set_defaults(func = run_pan_pipeline_withchecks)

    parser_core_pipeline = subparsers.add_parser('core-pipeline')
    parser_core_pipeline.set_defaults(func = run_core_pipeline_withchecks)
    parser_core_pipeline.add_argument("faapaths",
        help = "input file with paths to faa files of genomes")
    parser_core_pipeline.add_argument("outfolder",
        help = "output folder")
    parser_core_pipeline.add_argument("-d", "--method", default = "of_blast",
        choices = ["of_blast", "of_diamond"],
        help = "pangenome inference method [default: of_blast]")
    parser_core_pipeline.add_argument("-e", "--seeds", default = 50, type = int,
        help = "number of seed genomes to use [default 50]")
    parser_core_pipeline.add_argument("-f", "--seedfilter", default = 40,
        type = int,
        help = "minimum number of seed genomes where a candidate core gene "
            "should be present [default 40]")
    parser_core_pipeline.add_argument("-i", "--allfilter", default = 0.95,
        type = float,
        help = "minimum percentage of all genomes where a core gene should be "
            "present [default 0.95]")
    parser_core_pipeline.add_argument("-t", "--threads", default = 8,
        type = int, help = "number of threads to use [default 8]")
    parser_core_pipeline.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    args = parser.parse_args()

    return(args)

if __name__ == "__main__":

    args = parse_arguments()

    if not "func" in args:
        print_help()
        sys.exit()

    print_intro()

    logging.basicConfig(
        level = logging.INFO,
        format = '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt = '%d/%m %H:%M:%S'
    )
    logging.info("welcome to progenomics")

    if "cont" in args and args.cont and os.path.exists(args.outfolder):
        logging.info("continuing in existing output folder")
        check_outfile(os.path.join(args.outfolder, "progenomics.log"))
    else:
        logging.info("creating output folder and log file")
        check_outdir(args.outfolder)
        os.makedirs(args.outfolder, exist_ok = True)
        logging.info(f"output folder '{args.outfolder}' created")

    handler = logging.FileHandler(
        filename = os.path.join(args.outfolder, "progenomics.log"),
        mode = 'w'
    )
    handler.setFormatter(logging.getLogger().handlers[0].formatter)
    logging.getLogger().addHandler(handler)
    logging.info("log file created!")

    args.func(args)

    logging.info("progenomics out")
