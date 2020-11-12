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
    check_infile(args.faapaths)
    if args.method in ["O-B", "O-D"]:
        check_fastas(args.faapaths) # mmseqs can handle zipped fastas
    if not args.species is None:
        check_infile(args.species)

    logging.info("checking dependencies")
    if args.method in ["O-B", "O-D"]:
        check_tool("orthofinder")
    elif args.method == "S":
        check_tool("mmseqs", ["-h"])
    else:
        check_tool("mmseqs", ["-h"])
        check_tool("mafft", ["--help"]) # --help avoids mafft interactive mode

    run_pan(args)

def run_build_withchecks(args):

    logging.info("welcome to the build task")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_infile(args.pangenome)

    logging.info("checking dependencies")
    check_tool("hmmbuild", ["-h"])
    check_tool("mafft", ["--help"]) # --help avoids mafft interactive mode

    run_build(args)

def run_search_withchecks(args):

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
    check_infile(args.faapaths)
    check_infile(args.coregenome)
    if not args.ffnpaths is None:
        check_infile(args.ffnpaths)

    logging.info("checking dependencies")
    check_tool("mafft", ["--help"]) # --help avoids mafft interactive mode
    
    run_supermatrix(args)

def run_clust_withchecks(args):

    logging.info("welcome to the clust task")

    logging.info("checking arguments other than output folder")
    check_infile(args.fastapaths)
    check_infile(args.coregenome)

    logging.info("checking dependencies")
    check_tool("mmseqs", ["-h"])
    
    run_clust(args)

def run_pan_pipeline_withchecks(args):

    logging.info("welcome to the pan pipeline")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_fastas(args.faapaths)
    if not args.species is None:
        check_infile(args.species)

    logging.info("checking dependencies")
    check_tool("orthofinder")
    check_tool("hmmbuild", ["-h"])
    check_tool("mafft", ["--help"])
    check_tool("hmmsearch", ["-h"])
    
    run_pan_pipeline(args)

def run_core_pipeline_withchecks(args):

    logging.info("welcome to the core pipeline")

    logging.info("checking arguments other than output folder")
    check_infile(args.faapaths)
    check_fastas(args.faapaths)

    logging.info("checking dependencies")
    check_tool("orthofinder")
    check_tool("hmmbuild", ["-h"])
    check_tool("mafft", ["--help"])
    check_tool("hmmsearch", ["-h"])
    
    run_core_pipeline(args)
