import logging
import os
import subprocess
import sys

from Bio import SeqIO
from pathlib import Path
from utils import *

def check_tool(tool, arguments = []):
    devnull = open(os.devnull, 'w')
    try:
        subprocess.call([tool] + arguments, stdout = devnull, stderr = devnull)
    except FileNotFoundError:
        logging.error(f"{tool} not found")
        sys.exit(1)
    logging.info(f"{tool} found")

def check_infile(infile):
    if not os.path.isfile(infile):
        logging.error(f"input file '{infile}' not found")
        sys.exit(1)

def check_outfile(outfile):
    if not os.path.isfile(outfile):
        return()
    i = 0
    outfile_prev = outfile + f".prev{i}"
    while os.path.isfile(outfile_prev):
        i = i + 1
        outfile_prev = outfile + f".prev{i}"
    logging.info(f"output file '{outfile}' already exists; moving it to "
        f"'{outfile_prev}'")
    os.rename(outfile, outfile_prev)

def check_outdir(outdir, move_existing = True):
    if not os.path.exists(outdir):
        return()
    elif move_existing:
        i = 0
        outdir_prev = outdir + f".prev{i}"
        while os.path.exists(outdir_prev):
            i = i + 1
            outdir_prev = outdir + f".prev{i}"
        logging.info(f"output folder '{outdir}' already exists; moving it to "
            f"'{outdir_prev}'")
        os.rename(outdir, outdir_prev)
    else:
        logging.info(f"output folder '{outdir}' already exists")
        sys.exit(1)

def check_fasta(path):
    if not os.path.isfile(path):
        logging.error("one or more fasta files not found")
        sys.exit(1)
    with open_smart(path) as handle:
        try:
            fasta = SeqIO.parse(handle, "fasta")
            if not any(fasta):
                raise Exception()
        except Exception:
            logging.error("one or more fasta files are not in fasta format")
            sys.exit(1)

def check_fastas(path):
    filenames = []
    for fastapath in open(path):
        fastapath = fastapath.strip()
        filenames.append(filename_from_path(fastapath))
        check_fasta(fastapath)
    if not len(filenames) == len(set(filenames)):
        logging.error("names of fasta files are not unique")
        sys.exit(1)

def check_db(path):
    exts = [".h3f", ".h3i", ".h3m", ".h3p"]
    for db_file in [os.path.join(path, "hmm_db" + ext) for ext in exts]:
        if not os.path.isfile(db_file):
            logging.error("profile hmm database not found")
            sys.exit(1)
