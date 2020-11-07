import concurrent.futures
import logging
import os
import shutil
import subprocess
import sys

from Bio.Align.Applications import MafftCommandline
from pathlib import Path

from utils import *

def run_iqtree(fin_aln, dout_tree, threads, options):
    makedirs_smart(dout_tree)
    result = subprocess.call(["iqtree", "-te", "BIONJ", "-fixbr", "-s", 
        fin_aln, "-pre", f"{dout_tree}/tree", "-nt", str(threads)] + options, 
        stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    if result != 0:
        logging.error("something went wrong with iqtree; see log file "
            f"{dout_tree}/tree.log")
        sys.exit(1)

def run_mmseqs(arguments, logfout, skip_if_exists = "", threads = 1):
    if skip_if_exists != "":
        if os.path.exists(skip_if_exists):
            if os.path.getsize(skip_if_exists) > 0:
                logging.info("existing output detected - moving on")
                return()
    args = ["mmseqs"] + arguments
    if not arguments[0] in ["createdb", "convertmsa", "tsv2db"]:
        args = args + ["--threads", str(threads)]
    with open(logfout, "w") as loghout:
        result = subprocess.call(args, stdout = loghout, stderr = loghout)
    if result != 0:
        logging.error("something went wrong with mmseqs; see log file "
            f"{logfout}")
        sys.exit(1)

def run_orthofinder(faafins, dout, logfout, threads, engine = "blast"):
    for faafin in faafins:
        faaname = os.path.basename(faafin)
        shutil.copyfile(faafin, os.path.join(dout, faaname))
    with open(logfout, 'w') as loghout:
        with open(os.devnull, 'w') as devnull:
            result = subprocess.call(["orthofinder", "-og", "-S", engine, "-t",
                str(threads), "-f", dout], stdout = loghout, stderr = loghout)
    if result != 0:
        logging.error("something went wrong with orthofinder; see log file "
            "'orthofinder.log'")
        sys.exit(1)

def run_mafft(fin_aa_seqs, fout_aa_seqs_aligned, threads = 1):
    logging.debug(f"aligning {Path(fin_aa_seqs).stem}")
    mafft_cline = MafftCommandline(input = fin_aa_seqs, thread = threads)
    stdout, stderr = mafft_cline()
    with open(fout_aa_seqs_aligned, "w") as hout_aa_seqs_aligned:
        hout_aa_seqs_aligned.write(stdout)

def run_mafft_parallel(fins_aa_seqs, fouts_aa_seqs_aligned):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_mafft, fins_aa_seqs, fouts_aa_seqs_aligned)

def run_hmmbuild(fin_aa_seqs_aligned, fout_hmm):
    logging.debug(f"building profile for {Path(fin_aa_seqs_aligned).stem}")
    subprocess.run(["hmmbuild", fout_hmm, fin_aa_seqs_aligned],
        stdout = subprocess.PIPE)

def run_hmmbuild_parallel(fins_aa_seqs_aligned, fouts_hmms):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_hmmbuild,fins_aa_seqs_aligned, fouts_hmms)

def run_hmmpress(fins_hmms, dout_hmm_db):
    fout_hmm_db = dout_hmm_db + "/hmm_db"
    with open(fout_hmm_db, 'w') as hout_hmm_db:
        for fin_hmm in fins_hmms:
            with open(fin_hmm) as hin_hmm:
                hout_hmm_db.write(hin_hmm.read())
    subprocess.run(["hmmpress", fout_hmm_db], stdout = subprocess.PIPE)
    subprocess.run(["rm", fout_hmm_db])

def run_hmmsearch(din_hmm_db, fins_genomes, fout_domtbl):
    fin_hmm_db = os.path.join(din_hmm_db, "hmm_db")
    fout_temp_genomes = din_hmm_db + "/genomes.temp"
    with open(fout_temp_genomes, 'w') as hout_temp_genomes:
        for fin_genome in fins_genomes:
            with open(fin_genome) as hin_genome:
                hout_temp_genomes.write(hin_genome.read())
    subprocess.run(["hmmsearch", "-o", "/dev/null", "--domtblout", fout_domtbl,
        fin_hmm_db, fout_temp_genomes])
    subprocess.run(["rm", fout_temp_genomes])
