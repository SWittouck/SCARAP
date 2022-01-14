import logging
import os
import sys

from utils import *
from checkers import *
from tasks_composite import *
from tasks_core import *

def run_pan_withchecks(args):

    logging.info("welcome to the pan task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faapaths)
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

    run_pan(args)

def run_build_withchecks(args):

    logging.info("welcome to the build task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faapaths)
    check_infile(args.pangenome)
    faapaths = read_fastapaths(args.faapaths)
    n_genomes = len(faapaths)
    if args.min_genomes > n_genomes:
        args.min_genomes = n_genomes
        logging.info(f"min_genomes reduced to {args.min_genomes}, since that's "
            "the total number of genomes")

    logging.info("checking dependencies")
    check_tool("hmmbuild", ["-h"])
    check_mafft()

    run_build(args)

def run_search_withchecks(args):

    logging.info("welcome to the search task")

    logging.info("checking arguments other than output folder")
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

def run_supermatrix_withchecks(args):

    logging.info("welcome to the supermatrix task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faapaths)
    check_infile(args.coregenome)
    if not args.ffnpaths is None:
        check_fastas(args.ffnpaths)

    logging.info("checking dependencies")
    check_mafft() 
    
    run_supermatrix(args)

def run_clust_withchecks(args):

    logging.info("welcome to the clust task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.fastapaths)
    fastapaths = read_fastapaths(args.fastapaths)
    n_genomes = len(fastapaths)
    if args.max_clusters > n_genomes:
        args.max_clusters = n_genomes
        logging.info(f"max_clusters reduced to {args.max_clusters}, since "
            "that's the total number of genomes")
    check_infile(args.coregenome)
    if args.identity > 100:
        logging.error("identity should be between 0 and 1")
        sys.exit(1)
    elif args.identity > 1:
        args.identity = args.identity / 100
        logging.info(f"corrected identity value to {str(args.identity)}")
    if args.method == "median":
        logging.info("per-gene identity values will be aggregated using the "
            "median of all values")
    elif args.method == "mean":
        logging.info("per-gene identity values will be aggregated using "
            "the mean of all values")
    elif args.method[:4] == "mean":
        p = args.method[4:]
        try:
            p = int(p)
            if p > 100 or p <= 0: raise ValueError
        except ValueError:
            logging.error("method unknown")
            sys.exit(1)
        logging.info(f"per-gene identity values will be aggregated using "
            f"the mean of the middle {p}% of values")
    else:
        logging.error("method unknown")
        sys.exit(1)

    logging.info("checking dependencies")
    check_mmseqs()
    
    run_clust(args)

def run_fetch_withchecks(args):

    logging.info("welcome to the fetch task")

    logging.info("checking arguments other than output folder")
    check_fastas(args.fastapaths)
    check_infile(args.genes)
    
    run_fetch(args)

def run_pan_pipeline_withchecks(args):

    logging.info("welcome to the pan pipeline")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faapaths)
    if not args.species is None:
        check_infile(args.species)

    logging.info("checking dependencies")
    if args.method in ["O-B", "O-D"]:
        check_tool("orthofinder")
    else:
        check_mmseqs()
    check_mafft()
    check_tool("hmmbuild", ["-h"])
    check_tool("hmmsearch", ["-h"])
    
    run_pan_pipeline(args)

def run_core_pipeline_withchecks(args):

    logging.info("welcome to the core pipeline")

    logging.info("checking arguments other than output folder")
    check_fastas(args.faapaths)
    n_genomes = sum([1 for line in open(args.faapaths, "r")])
    if args.seeds > n_genomes:
        seedfilter_new = n_genomes * args.seedfilter // args.seeds # result is still integer
        args.seeds = n_genomes
        logging.info(f"number of seeds reduced to {args.seeds}, since that's "
            "the number of genomes")
        args.seedfilter = seedfilter_new
        logging.info(f"seedfilter set to {args.seedfilter}")
    if args.allfilter > 100:
        logging.error("allfilter should be between 0 and 1")
        sys.exit(1)
    elif args.allfilter > 1:
        args.allfilter = args.allfilter / 100
        logging.info(f"corrected allfilter value to {str(args.allfilter)}")

    logging.info("checking dependencies")
    if args.method in ["O-B", "O-D"]:
        check_tool("orthofinder")
    else:
        check_mmseqs()
    check_mafft()
    check_tool("hmmbuild", ["-h"])
    check_tool("hmmsearch", ["-h"])
    
    run_core_pipeline(args)
