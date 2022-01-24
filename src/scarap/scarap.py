#! /usr/bin/env python3

# This is the main script of SCARAP; it only contains the commandline
# interface.

__author__ = "Stijn Wittouck"
__version__ = "0.3.2"

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
    scarap [-h] <task> <task-specific arguments>
TASKS
    pan           --> infer a pangenome from a set of faa files
    build         --> build a profile HMM database for a core/pangenome
    search        --> search query genes in a core/pangenome database
    checkgenomes  --> assess the quality of genomes in a core genome
    checkgroups   --> assess the quality of orthogroups in a core genome
    filter        --> filter the genomes/orthogroups in a pangenome
    supermatrix   --> construct a concatenated core orthogroup alignment from a
                      core genome
    sample        --> sample a subset of representative genomes
    fetch         --> fetch sequences and store in fasta per orthogroup
PIPELINES
    core          --> infer a core genome from a set of faa files
    pan-pipeline  --> infer a pangenome, build a profile HMM database and train
                      score cutoffs from a set of faa files
DOCUMENTATION
    https://github.com/swittouck/scarap\
'''

    print(message.format(__version__))

def print_intro():

    message = '''\

This is SCARAP version {0}
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
        help = "file with paths to or folder with faa files of genomes")
    parser_pan_parent.add_argument("-d", "--method", default = "FH",
        choices = ["H", "FH", "FT", "H-nl", "T-nl", "P", "O-B", "O-D", "S"],
        help = "pangenome inference method [default: FH]")
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
        help = "file with paths to or folder with faa files of genomes")
    parser_build.add_argument("pangenome",
        help = "input file with pangenome")
    parser_build.add_argument("outfolder",
        help = "output folder for hmm database")
    parser_build.add_argument("-p", "--core_prefilter", default = 0,
        type = float,
        help = "minimum relative frequency of single-copy presence to be used "
            "for initial database construction [default 0]")
    parser_build.add_argument("-f", "--core_filter", default = 0,
        type = float,
        help = "minimum relative frequency of single-copy presence to be used "
            "for final database selection [default 0]")
    parser_build.add_argument("-m", "--max_cores", default = 0,
        type = int,
        help = "maximum number of core genes to retrieve (0 = no maximum) " 
            "[default 0]")
    parser_build.add_argument("-t", "--threads", default = 8, type = int,
        help = "number of threads to use [default 8]")
    parser_build.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_search = subparsers.add_parser('search')
    parser_search.set_defaults(func = run_search_withchecks)
    parser_search.add_argument("qpaths",
        help = "file with paths to or folder with files containing amino acid "
            "sequences of query genes")
    parser_search.add_argument("db",
        help = "input folder with hmm database")
    parser_search.add_argument("outfolder",
        help = "output folder for search hits")
    parser_search.add_argument("-t", "--threads", default = 8, type = int,
        help = "number of threads to use [default 8]")
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
        help = "file with paths to or folder with faa files of genomes")
    parser_supermatrix.add_argument("coregenome",
        help = "file with core genome")
    parser_supermatrix.add_argument("outfolder",
        help = "output file for supermatrix fasta file based on core genome")
    parser_supermatrix.add_argument("-n", "--ffnpaths",
        help = "file with paths to or folder with ffn files of genomes; if "
            "given, a nucleotide supermatrix will be constructed in addition "
            "to the amino acid suprmatrix")
    parser_supermatrix.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_sample = subparsers.add_parser('sample')
    parser_sample.set_defaults(func = run_sample_withchecks)
    parser_sample.add_argument("fastapaths",
        help = "file with paths to or folder with ffn/faa files of genomes")
    parser_sample.add_argument("coregenome",
        help = "file with core genome")
    parser_sample.add_argument("outfolder",
        help = "output folder")
    parser_sample.add_argument("-m", "--max_genomes", default = 0, type = int,
        help = "maximum number of genomes to sample (0 = no maximum) "
            "[default 0]")
    parser_sample.add_argument("-i", "--identity", default = 1, 
        type = float,
        help = "maximum sequence identity between sampled genomes [default 1]")
    parser_sample.add_argument("-x", "--exact", action = "store_true",
        help = "perform full alignments [default False]")
    parser_sample.add_argument("-t", "--threads", default = 8, type = int,
        help = "number of threads to use [default 8]")
    parser_sample.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")
    parser_sample.add_argument("-d", "--method", default = "mean",
        help = "genome-genome comparison method [default: mean]")

    parser_fetch = subparsers.add_parser('fetch')
    parser_fetch.set_defaults(func = run_fetch_withchecks)
    parser_fetch.add_argument("fastapaths",
        help = "file with paths to or folder with ffn/faa files of genomes")
    parser_fetch.add_argument("genes",
        help = "file with genes (columns gene, genome, orthogroup)")
    parser_fetch.add_argument("outfolder",
        help = "output folder for fasta files")
    parser_fetch.add_argument("-c", "--cont", action = "store_true",
        help = "continue in existing output folder [default False]")

    parser_pan_pipeline = subparsers.add_parser('pan-pipeline',
        parents = [parser_pan_parent])
    parser_pan_pipeline.set_defaults(func = run_pan_pipeline_withchecks)

    parser_core = subparsers.add_parser('core')
    parser_core.set_defaults(func = run_core_withchecks)
    parser_core.add_argument("faapaths",
        help = "file with paths to or folder with faa files of genomes")
    parser_core.add_argument("outfolder",
        help = "output folder")
    parser_core.add_argument("-d", "--method", default = "FH",
        choices = ["H", "FH", "H-nl", "T-nl", "P", "O-B", "O-D", "S"],
        help = "pangenome inference method [default: FH]")
    parser_core.add_argument("-e", "--seeds", default = 100, type = int,
        help = "number of seed genomes to use [default 100]")
    parser_core.add_argument("-p", "--core_prefilter", default = 0.90,
        type = float,
        help = "minimum relative frequency of single-copy presence to be used "
            "for initial database construction [default 0.90]")
    parser_core.add_argument("-f", "--core_filter", default = 0.95,
        type = float,
        help = "minimum relative frequency of single-copy presence to be used "
            "for final database selection [default 0.95]")
    parser_core.add_argument("-m", "--max_cores", default = 0,
        type = int,
        help = "maximum number of core genes to retrieve (0 = no maximum) " 
            "[default 0]")
    parser_core.add_argument("-t", "--threads", default = 8,
        type = int, help = "number of threads to use [default 8]")
    parser_core.add_argument("-c", "--cont", action = "store_true",
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
    logging.info("welcome to SCARAP")

    if "cont" in args and args.cont and os.path.exists(args.outfolder):
        logging.info("continuing in existing output folder")
        check_outfile(os.path.join(args.outfolder, "SCARAP.log"))
    else:
        logging.info("creating output folder and log file")
        check_outdir(args.outfolder)
        os.makedirs(args.outfolder, exist_ok = True)
        logging.info(f"output folder '{args.outfolder}' created")

    handler = logging.FileHandler(
        filename = os.path.join(args.outfolder, "SCARAP.log"),
        mode = 'w'
    )
    handler.setFormatter(logging.getLogger().handlers[0].formatter)
    logging.getLogger().addHandler(handler)
    logging.info("log file created!")

    args.func(args)

    logging.info("SCARAP out")
