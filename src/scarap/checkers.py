import logging
import os
import re
import subprocess
import sys
from collections import Counter

from Bio import SeqIO
from pathlib import Path
from scarap.utils import *

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
        logging.warning("SCARAP has been tested with MAFFT v7.310 or newer")
    
def check_mmseqs():
    try:
        res = subprocess.run(["mmseqs", "-h"], stdout = subprocess.PIPE)
    except FileNotFoundError:
        logging.error(f"MMseqs2 not found")
        sys.exit(1)
    r = re.compile("MMseqs2 Version: ([a-f0-9.]{5})")
    try:
        version = r.search(res.stdout.decode()).group(1)
    except AttributeError:
        version = "unknown"
    releases_tested = {"e1a1c": "11", "113e3": "12", "45111": "13"}
    if version.split(".")[0] in releases_tested.values():
    # new version shows human readable version
        release = version.split(".")[0]
    else:
    # Old version showed the first 5 characters of the commit sha for the version
        release = releases_tested.get(version, "unknown")
    logging.info(f"detected MMseqs2 version {version} (release {release})")
    if not release in releases_tested.values():
        logging.warning("SCARAP has only been tested with MMseqs2 "
            "releases 11, 12 and 13")

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
            logging.error(f"file {path} is not in fasta format")
            sys.exit(1)

# works for dir with fastas or file with fastapaths
def check_fastas(path):
    # when path is a file
    if os.path.isfile(path):
        # error if file with fastapaths is empty
        if os.stat(path).st_size == 0:
            logging.error("fastapaths file is empty")
            sys.exit(1)
        # store fastapaths in list
        fastapaths = [fp.strip() for fp in open(path)]
        # error if fasta filenames not unique
        filenames = [filename_from_path(fp) for fp in fastapaths]
        if not len(filenames) == len(set(filenames)):
            # Check for the culprit
            offenders = [x for x,c in Counter(filenames).items() if c > 1]
            logging.error(f"names of fasta files are not unique. Offending values are: {offenders}")
            sys.exit(1)
    # when path is a folder
    elif os.path.isdir(path):
        # store fastapaths in list
        fastapaths = [os.path.join(path, file) for file in os.listdir(path)]
        fastapaths = [fp for fp in fastapaths if os.path.isfile(fp)]
        extensions = ("fasta", "fa", "faa", "ffn", "fasta.gz", "fa.gz", 
            "faa.gz", "ffn.gz")
        fastapaths = [fp for fp in fastapaths if fp.endswith(extensions)]
    # error when path doesn't exist
    else: 
        logging.error(f"input file/folder '{path}' not found")
        sys.exit(1)
    # check individual fasta files
    for fastapath in fastapaths:
        check_fasta(fastapath)

def check_db(path):
    din_alis = os.path.join(path, "alignments")
    fin_cutoffs = os.path.join(path, "cutoffs.csv")
    if not os.path.isdir(din_alis):
        logging.error("no folder with alignments found")
        sys.exit(1)
    if not os.path.isfile(fin_cutoffs):
        logging.error("no file with score cutoffs found")
        sys.exit(1)
