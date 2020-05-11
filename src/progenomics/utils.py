import gzip
import os
import pandas as pd
import shutil

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
