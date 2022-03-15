#! /usr/bin/env python3

# This is the main script of SCARAP; it only contains the commandline
# interface.

__author__ = "Stijn Wittouck"
__version__ = "0.3.2"

import argparse
import logging
import sys

from utils import *
from module_wrappers import *

def print_help():

    message = '''\
VERSION
    {0}
AUTHORS
    Stijn Wittouck (development)
    Sarah Lebeer (supervision)
USAGE
    scarap [-h] <task> <task-specific arguments>
MODULES
    pan           --> infer a pangenome from a set of faa files
    core          --> infer a core genome from a set of faa files
    build         --> build a profile database for a core/pangenome
    search        --> search query genes in a profile database
    checkgenomes  --> assess the quality of genomes in a core genome
    checkgroups   --> assess the quality of orthogroups in a core genome
    filter        --> filter the genomes/orthogroups in a pangenome
    concat        --> construct a concatenated core orthogroup alignment from a
                      core genome
    sample        --> sample a subset of representative genomes
    fetch         --> fetch sequences and store in fasta per orthogroup
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
    
    # help messages for positional arguments (alphabetically)
    h_coregenome = "input file with a core genome"
    h_db = "input database containing alignments and profile score cutoffs"
    h_faa_files = "folder with one fasta file with amino acid sequences (faa "\
        "file) per genome, or file with paths to faa files"
    h_fasta_files = "folder with one fasta file with amino acid sequences "\
        "(faa file) or one fasta file with nucleic acid sequences (ffn file) "\
        "per genome, or a file with paths to faa/ffn files"
    h_outfolder = "output folder"
    h_pangenome = "input file with pangenome"
    
    # help messages for "optional" arguments (alphabetically)
    h_core_prefilter = "minimum relative frequency of single-copy presence to "\
        "be used for initial database construction"
    h_cont = "continue in existing output folder [default False]"
    h_core_filter = "minimum relative frequency of single-copy presence to be "\
        " used for"
    h_exact = "perform full alignments [default False]"
    h_ffn_files = "file with paths to or folder with ffn files of genomes; if "\
        "given, a nucleotide supermatrix will be constructed in addition "\
        "to the amino acid suprmatrix"
    h_genomes = "input file with genomes to extract from pangenome"
    h_identity = "maximum sequence identity between sampled genomes [default 1]"
    h_max_core_genes = "maximum number of core genes to retrieve (0 = no "\
        "maximum) [default 0]"
    h_max_genomes = "maximum number of genomes to sample (0 = no maximum) "\
        "[default 0]"
    h_method = "pangenome inference method [default: FH]"
    h_method_sample = "genome-genome comparison method [default: mean]"
    h_orthogroups = "input file with orthogroups to extract from pangenome"
    h_seeds = "number of seed genomes to use [default 100]"
    h_species = "input file with species of genomes; if given, a hierarchical "\
        "pangenome strategy will be used for final database selection"
    h_threads = "number of threads to use"
    
    # other interface components
    method_choices = ["H", "FH", "FT", "H-nl", "T-nl", "P", "O-B", "O-D", "S"]
    
    parser_pan = subparsers.add_parser('pan')
    parser_pan.set_defaults(func = run_pan_withchecks)
    parser_pan.add_argument("faa_files", metavar = "faa-files", 
        help = h_faa_files)
    parser_pan.add_argument("outfolder", help = h_outfolder)
    parser_pan.add_argument("-d", "--method", default = "FH", 
        choices = method_choices, help = h_method)
    parser_pan.add_argument("-s", "--species", help = h_species)
    parser_pan.add_argument("-t", "--threads", default = 8, type = int, 
        help = h_threads)
    parser_pan.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_core = subparsers.add_parser('core')
    parser_core.set_defaults(func = run_core_withchecks)
    parser_core.add_argument("faa_files", metavar = "faa-files", 
        help = h_faa_files)
    parser_core.add_argument("outfolder", help = h_outfolder)
    parser_core.add_argument("-d", "--method", default = "FH", 
        choices = method_choices, help = h_method)
    parser_core.add_argument("-e", "--seeds", default = 100, type = int,
        help = h_seeds)
    parser_core.add_argument("-p", "--core-prefilter", default = 0.90,
        type = float, help = h_core_prefilter + " [default 0.90]")
    parser_core.add_argument("-f", "--core-filter", default = 0.95,
        type = float, help = h_core_filter + " [default 0.95]")
    parser_core.add_argument("-m", "--max-core-genes", default = 0,
        type = int, help = h_max_core_genes)
    parser_core.add_argument("-t", "--threads", default = 8,
        type = int, help = h_threads)
    parser_core.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_build = subparsers.add_parser('build')
    parser_build.set_defaults(func = run_build_withchecks)
    parser_build.add_argument("faa_files", metavar = "faa-files", 
        help = h_faa_files)
    parser_build.add_argument("pangenome", help = h_pangenome)
    parser_build.add_argument("outfolder", help = h_outfolder)
    parser_build.add_argument("-p", "--core-prefilter", default = 0,
        type = float, help = h_core_prefilter + " [default 0]")
    parser_build.add_argument("-f", "--core-filter", default = 0,
        type = float, help = h_core_filter + " [default 0]")
    parser_build.add_argument("-m", "--max-core-genes", default = 0,
        type = int, help = h_max_core_genes)
    parser_build.add_argument("-t", "--threads", default = 8,
        type = int, help = h_threads)
    parser_build.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_search = subparsers.add_parser('search')
    parser_search.set_defaults(func = run_search_withchecks)
    parser_search.add_argument("faa_files", metavar = "faa-files", 
        help = h_faa_files)
    parser_search.add_argument("db", help = h_db)
    parser_search.add_argument("outfolder", help = h_outfolder)
    parser_search.add_argument("-t", "--threads", default = 8,
        type = int, help = h_threads)
    parser_search.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_concat = subparsers.add_parser('concat')
    parser_concat.set_defaults(func = run_concat_withchecks)
    parser_concat.add_argument("faa_files", metavar = "faa-files", 
        help = h_faa_files)
    parser_concat.add_argument("pangenome", help = h_pangenome)
    parser_concat.add_argument("outfolder", help = h_outfolder)
    parser_concat.add_argument("-f", "--core-filter", default = 0,
        type = float, help = h_core_filter + " [default 0]")
    parser_concat.add_argument("-m", "--max-core-genes", default = 0,
        type = int, help = h_max_core_genes)
    parser_concat.add_argument("-n", "--ffn-files", help = h_ffn_files)
    parser_concat.add_argument("-t", "--threads", default = 8,
        type = int, help = h_threads)
    parser_concat.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_sample = subparsers.add_parser('sample')
    parser_sample.set_defaults(func = run_sample_withchecks)
    parser_sample.add_argument("fasta_files", metavar = "fasta-files", 
        help = h_fasta_files)
    parser_sample.add_argument("pangenome", help = h_pangenome)
    parser_sample.add_argument("outfolder", help = h_outfolder)
    parser_sample.add_argument("-m", "--max-genomes", default = 0, type = int,
        help = h_max_genomes)
    parser_sample.add_argument("-i", "--identity", default = 1, type = float, 
        help = h_identity)
    parser_sample.add_argument("-d", "--method", default = "mean",
        help = h_method_sample)
    parser_sample.add_argument("-x", "--exact", action = "store_true",
        help = h_exact)
    parser_sample.add_argument("-t", "--threads", default = 8,
        type = int, help = h_threads)
    parser_sample.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_checkgenomes = subparsers.add_parser("checkgenomes")
    parser_checkgenomes.set_defaults(func = run_checkgenomes_withchecks)
    parser_checkgenomes.add_argument("coregenome", help = h_coregenome)
    parser_checkgenomes.add_argument("outfolder", help = h_outfolder)
    parser_checkgenomes.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_checkgroups = subparsers.add_parser("checkgroups")
    parser_checkgroups.set_defaults(func = run_checkgroups_withchecks)
    parser_checkgroups.add_argument("coregenome", help = h_coregenome)
    parser_checkgroups.add_argument("outfolder", help = h_outfolder)
    parser_checkgroups.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_filter = subparsers.add_parser('filter')
    parser_filter.set_defaults(func = run_filter_withchecks)
    parser_filter.add_argument("pangenome", help = h_pangenome)
    parser_filter.add_argument("outfolder", help = h_outfolder)
    parser_filter.add_argument("-g", "--genomes", help = h_genomes)
    parser_filter.add_argument("-o", "--orthogroups", help = h_orthogroups)
    parser_filter.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

    parser_fetch = subparsers.add_parser('fetch')
    parser_fetch.set_defaults(func = run_fetch_withchecks)
    parser_fetch.add_argument("fasta_files", metavar = "fasta-files", 
        help = h_fasta_files)
    parser_fetch.add_argument("pangenome", help = h_pangenome)
    parser_fetch.add_argument("outfolder", help = h_outfolder)
    parser_fetch.add_argument("-c", "--cont", action = "store_true",
        help = h_cont)

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
