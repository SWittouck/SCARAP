import logging
import os
import sys

from scarap.utils import *
from scarap.checkers import *
from scarap.modules import *

####################
# helper functions #
####################

def correct_freq(freq, name):
    if freq > 100:
        logging.error(f"{name} should be between 0 and 1")
        sys.exit(1)
    elif freq > 1:
        freq = freq / 100
        logging.info(f"corrected {name} value to {str(freq)}")
    return(freq)

def process_reps(args):

    if args.min_reps < 2:
        logging.error("at least two representative sequences are needed for "
            "family splitting")

    if (args.max_reps == 0):
        args.max_reps = round(args.min_reps * 1.25)
        logging.info(f"maximum number of representatives set to "
            f"{args.max_reps}")

    if (args.max_align == 0):
        args.max_align = round(args.max_reps * 16)
        logging.info(f"maximum number of sequences to align set to "
            f"{args.max_align}")

    if args.min_reps > args.max_reps:
        logging.error("the minimum number of representatives must be smaller "
        "than or equal to the maximum number")

    if args.max_reps > args.max_align:
        logging.error("the maximum number of representatives must be smaller "
        "than or equal to the maximum number of sequences to align")

    return(args)

###################
# module wrappers #
###################

def run_pan_withchecks(args):

    logging.info("welcome to the pan task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faa_files)
    if not args.species is None:
        check_infile(args.species)

    logging.info("checking dependencies")
    if args.method in ["O-B", "O-D"]:
        check_tool("orthofinder")
    elif args.method == "S":
        check_mmseqs()
    else:
        check_mmseqs()
        check_mafft()

    args = process_reps(args)

    run_pan(args)

def run_build_withchecks(args):

    logging.info("welcome to the build task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faa_files)
    check_infile(args.pangenome)
    faapaths = read_fastapaths(args.faa_files)
    args.core_prefilter = correct_freq(args.core_prefilter, "core prefilter")
    args.core_filter = correct_freq(args.core_filter, "core filter")
    
    # give warning if core filter or max core genes is given without prefilter
    if args.core_filter != 0 or args.max_core_genes != 0:
        if args.core_prefilter == 0:
            logging.warning("a core filter or maximum number of core genes was "
                "given without a core prefilter - this can be needlessly slow")

    logging.info("checking dependencies")
    check_mmseqs()
    check_mafft()

    run_build(args)

def run_search_withchecks(args):

    logging.info("welcome to the search task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faa_files)
    check_db(args.db)

    logging.info("checking dependencies")
    check_mmseqs()

    run_search(args)

def run_checkgenomes_withchecks(args):

    logging.info("welcome to the checkgenomes task")

    logging.info("checking arguments other than output folder")
    check_infile(args.coregenome)

    run_checkgenomes(args)

def run_checkgroups_withchecks(args):

    logging.info("welcome to the checkgroups task")

    logging.info("checking arguments other than output folder")
    check_infile(args.coregenome)

    run_checkgroups(args)

def run_filter_withchecks(args):

    logging.info("welcome to the filter task")

    logging.info("checking arguments other than output folder")
    check_infile(args.pangenome)
    if not args.genomes is None:
        check_infile(args.genomes)
    if not args.orthogroups is None:
        check_infile(args.orthogroups)

    run_filter(args)

def run_concat_withchecks(args):

    logging.info("welcome to the concat task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faa_files)
    check_infile(args.pangenome)
    if not args.ffn_files is None:
        check_fastas(args.ffn_files)
    args.core_filter = correct_freq(args.core_filter, "core filter")

    logging.info("checking dependencies")
    check_mafft() 
    
    run_concat(args)

def run_sample_withchecks(args):

    logging.info("welcome to the sample task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.fasta_files)
    fastapaths = read_fastapaths(args.fasta_files)
    n_genomes = len(fastapaths)
    if args.max_genomes > n_genomes:
        args.max_genomes = n_genomes
        logging.info(f"max_genomes reduced to {args.max_genomes}, since "
            "that's the total number of genomes")
    check_infile(args.pangenome)
    if args.identity > 100:
        logging.error("identity should be between 0 and 1")
        sys.exit(1)
    elif args.identity > 1:
        args.identity = args.identity / 100
        logging.info(f"corrected identity value to {str(args.identity)}")
    if args.method == "median":
        logging.info("whole-genome ANI will be calculated as the median of all "
            "per-gene identity values")
    elif args.method == "mean":
        logging.info("whole-genome ANI will be calculated as the mean of all "
            "per-gene identity values")
    elif args.method[:4] == "mean":
        p = args.method[4:]
        try:
            p = int(p)
            if p > 100 or p <= 0: raise ValueError
        except ValueError:
            logging.error("method unknown")
            sys.exit(1)
        logging.info(f"whole-genome ANI will be calculated as the mean of the "
            f"middle {p}% of per-gene identity values")
    else:
        logging.error("method unknown")
        sys.exit(1)

    logging.info("checking dependencies")
    check_mmseqs()
    
    run_sample(args)

def run_fetch_withchecks(args):

    logging.info("welcome to the fetch task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.fasta_files)
    check_infile(args.pangenome)
    
    run_fetch(args)

def run_core_withchecks(args):

    logging.info("welcome to the core pipeline")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faa_files)
    fastapaths = read_fastapaths(args.faa_files)
    n_genomes = len(fastapaths)
    if args.seeds > n_genomes:
        args.seeds = n_genomes
        logging.info(f"number of seeds reduced to {args.seeds}, since that's "
            "the number of genomes")
    args.core_prefilter = correct_freq(args.core_prefilter, "core prefilter")
    args.core_filter = correct_freq(args.core_filter, "core filter")

    logging.info("checking dependencies")
    if args.method in ["O-B", "O-D"]:
        check_tool("orthofinder")
    else:
        check_mmseqs()
    check_mafft()

    args = process_reps(args)
    
    run_core(args)
