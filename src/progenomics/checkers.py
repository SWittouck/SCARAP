import logging
import os
import re
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
    
def check_mafft():
    try:
        res = subprocess.run(["mafft", "--version"], stderr = subprocess.PIPE)
    except FileNotFoundError:
        logging.error(f"MAFFT not found")
        sys.exit(1)
    r = re.compile("[0-9]+\\.[0-9]+")
    version = r.search(res.stderr.decode()).group()
    logging.info(f"detected MAFFT v{version}")
    if float(version) < 7.310:
        logging.warning("progenomics has been tested with MAFFT v7.310 or newer")
    
def check_mmseqs():
    try:
        res = subprocess.run(["mmseqs", "-h"], stdout = subprocess.PIPE)
    except FileNotFoundError:
        logging.error(f"MMseqs2 not found")
        sys.exit(1)
    r = re.compile("MMseqs2 Version: ([a-f0-9]{5})")
    version = r.search(res.stdout.decode()).group(1)
    releases_tested = {"e1a1c": "11", "113e3": "12"}
    release = releases_tested.get(version, "unknown")
    logging.info(f"detected MMseqs2 version {version} (release {release})")
    if not version in releases_tested.keys():
        logging.warning("progenomics has only been tested with MMseqs2 "
            "releases 11 and 12")

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
    if os.stat(path).st_size == 0:
        logging.error("fastapaths file is empty")
        sys.exit(1)
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
