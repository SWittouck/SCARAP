import gzip
import os
import pandas as pd
import shutil

def create_prefdb(TDB, PDB):
    """Initialize a fake prefilter database.
    
    Create a fake prefilter database to be used for searching all queries
    against all targets.
    
    See https://github.com/soedinglab/MMseqs2/wiki#how-to-create-a-fake-prefiltering-for-all-vs-all-alignments
    
    Args:
        TDB (str): The path to the target database, but relative to the prefilter 
            database!!
        PDB (str): The path to the prefilter database. 
    """
    os.symlink(f"{TDB}.index", PDB)
    with open(f"{PDB}.dbtype", "w") as hout_PDB:
        chars = [chr(7)] + [chr(0)] * 3
        for char in chars:
            hout_PDB.write(char)
    
def update_prefdb(QDB, TDB, PDB):
    """Setup the queries for a fake prefilter database. 
    
    Set the queries of a fake prefilter database to be used for searching all
    queries against all targets.
    
    See https://github.com/soedinglab/MMseqs2/wiki#how-to-create-a-fake-prefiltering-for-all-vs-all-alignments
    
    Args: 
        QDB (str): The path to the query database.
        TDB (str): The path to the target database.
        PDB (str): The path to the prefilter database.
    """
    index_size = os.path.getsize(f"{TDB}.index")
    open(f"{PDB}.index", "w").close()
    with open(f"{QDB}.index", "r") as hin_QDB:
        with open(f"{PDB}.index", "a") as hout_PDB:
            for query_line in hin_QDB:
                values = query_line.strip().split("\t")
                towrite = "\t".join([values[0], "0", str(index_size)])
                hout_PDB.write(f"{towrite}\n")
        
def padded_counts(n):
    d = len(str(n))
    counts = [f"{str(i + 1).zfill(d)}" for i in range(n)]
    return(counts)

def makedirs_smart(dout):
    if os.path.exists(dout):
        shutil.rmtree(dout)
    os.makedirs(dout, exist_ok = True)

def open_smart(filename, mode = "rt"):
    if filename.endswith(".gz"):
        return(gzip.open(filename, mode))
    else:
        return(open(filename, mode))

def make_paths(filenames, folder, extension):
    paths = [os.path.join(folder, filename + extension) for filename in
        filenames]
    return(paths)

def filename_from_path(path):
    if (path.endswith(".gz")):
        path = path[:-3]
    filename = os.path.basename(path)
    filename = os.path.splitext(filename)[0]
    return(filename)

def write_tsv(genes, fout):
    genes.to_csv(fout, sep = "\t", index = False, header = False)

def write_lines(lines, fout):
    lines = pd.DataFrame({"line": lines})
    lines.to_csv(fout, index = False, header = False)

def read_lines(fin):
    with open(fin) as hin:
        lines = [line.strip() for line in hin.readlines()]
    return(lines)
  
def listpaths(folder):
    paths = [os.path.join(folder, f) for f in os.listdir(folder)]
    return(paths)
