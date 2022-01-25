import logging
import os
import shutil

from utils import *
from readerswriters import *
from computers import *
from callers import *

# used by the build and search modules
def run_profilesearch(fins_faas, fins_alis, fout_hits, dout_tmp, threads):
  
    # define tmp paths 
    dout_mmseqs = os.path.join(dout_tmp, "mmseqs2")
    dout_logs = os.path.join(dout_tmp, "mmseqs2_logs")
    dout_rubbish = os.path.join(dout_tmp, "rubbish")
    fout_sto = os.path.join(dout_tmp, "alis.sto")

    # create tmp subfolders
    for dir in [dout_mmseqs, dout_logs, dout_rubbish]:
        os.makedirs(dir, exist_ok = True)

    logging.info("converting alignments from fasta to stockholm format")
    fastas2stockholm(fins_alis, fout_sto)
    
    logging.info("creating mmseqs2 database from faa files")
    run_mmseqs(["createdb"] + fins_faas + [f"{dout_mmseqs}/faadb"], 
        f"{dout_logs}/create_faadb.log")
    
    logging.info("creating mmseqs2 database with profiles of orthogroups")
    run_mmseqs(["convertmsa", fout_sto, f"{dout_mmseqs}/msadb",
        "--identifier-field", "0"], f"{dout_logs}/create_msadb.log")
    run_mmseqs(["msa2profile", f"{dout_mmseqs}/msadb", 
        f"{dout_mmseqs}/profiledb", "--match-mode", "1"], 
        f"{dout_logs}/create_msadb.log")
    
    logging.info("searching faas against profiles") 
    run_mmseqs(["search", f"{dout_mmseqs}/faadb", f"{dout_mmseqs}/profiledb", 
        f"{dout_mmseqs}/resultdb", f"{dout_rubbish}"], 
        f"{dout_logs}/search.log", threads = threads)
    run_mmseqs(["createtsv", f"{dout_mmseqs}/faadb", f"{dout_mmseqs}/profiledb",
        f"{dout_mmseqs}/resultdb", fout_hits, "--full-header"], 
        f"{dout_logs}/create_tsv.log")
        
    # remove tmp
    shutil.rmtree(dout_tmp)
