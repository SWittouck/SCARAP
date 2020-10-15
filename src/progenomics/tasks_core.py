import logging
import os

from argparse import Namespace
from utils import *
from readerswriters import *
from computers import *
from callers import *
from pan import *

def run_pan(args):
    if "species" in args and not args.species is None:
        run_pan_hier(args)
    else:
        run_pan_nonhier(args)

def run_pan_nonhier(args):

    pangenomefout = os.path.join(args.outfolder, "pangenome.tsv")
    if os.path.isfile(pangenomefout):
        logging.info("existing pangenome detected - moving on")
        return()
    
    faafins = read_lines(args.faapaths)
    
    if args.method in ["O-B", "O-D"]:

        logging.info("creating orthofinder directory")
        dir_orthofinder = os.path.join(args.outfolder, "orthofinder")
        os.makedirs(dir_orthofinder, exist_ok = True)

        logging.info("running orthofinder")
        logfile = os.path.join(args.outfolder, "orthofinder.log")
        
        if args.method == "O-B":
            engine = "blast"
        else:
            engine = "diamond"
        run_orthofinder(faafins, dir_orthofinder, logfile, args.threads, 
            engine)
    
        logging.info("creating tidy pangenome file")
        pangenome = read_pangenome_orthofinder(dir_orthofinder)
        write_tsv(pangenome, pangenomefout)
        
    else:
      
        logging.info(f"constructing pangenome with {args.method} strategy")
        infer_pangenome(faafins, args.method, args.outfolder, args.threads)

def run_pan_hier(args):

    logging.info("processing species file")
    genomes_species = read_species(args.species)
    genomes_species.species = [species.replace(" ", "_") for species in
        genomes_species.species]
    speciesdict = {}
    for ix, row in genomes_species.iterrows():
        speciesdict.setdefault(row.species, []).append(row.genome)

    logging.info("processing faapaths")
    faafins = read_lines(args.faapaths)
    genomedict = {}
    for faafin in faafins:
        genome = filename_from_path(faafin)
        genomedict[genome] = faafin

    logging.info("started building pangenomes of species")
    speciespansdio = os.path.join(args.outfolder, "speciespangenomes")
    os.makedirs(speciespansdio, exist_ok = True)
    reprpaths = []
    speciespanfios = []
    for species, genomes in speciesdict.items():
        logging.info(f"inferring pangenome of {species}")
        faapaths = [genomedict[genome] for genome in genomes]
        dout = os.path.join(speciespansdio, species)
        os.makedirs(dout, exist_ok = True)
        faapathsfio = os.path.join(dout, "faapaths.txt")
        write_lines(faapaths, faapathsfio)
        run_pan_nonhier(Namespace(faapaths = faapathsfio, outfolder = dout,
            threads = args.threads, method = "FH"))
        speciespanfio = os.path.join(dout, "pangenome.tsv")
        speciespanfios.append(speciespanfio)
        reprfio = os.path.join(dout, species + ".faa")
        reprpaths.append(reprfio)
        speciespan = read_genes(speciespanfio)
        collapse_pangenome(speciespan, faapathsfio, reprfio, species,
            os.path.join(dout, "temp"))
    reprpathsfio = os.path.join(args.outfolder, "reprpaths.txt")
    write_lines(reprpaths, reprpathsfio)

    logging.info("started building metapangenome using representatives")
    metapandio = os.path.join(args.outfolder, "metapangenome")
    os.makedirs(metapandio, exist_ok = True)
    run_pan_nonhier(Namespace(faapaths = reprpathsfio,
        outfolder = metapandio, threads = args.threads,
        method = "FH"))

    logging.info("started inflating metapangenome with species pangenomes")
    speciespans = [read_genes(panfin) for panfin in speciespanfios]
    speciespans = pd.concat(speciespans)
    speciespans = speciespans.rename(columns = {"orthogroup": "speciesfam"})
    speciespans = pd.merge(speciespans, genomes_species)
    speciespans.speciesfam = [species + "-" + speciesfam for species,
        speciesfam in zip(speciespans.species, speciespans.speciesfam)]
    metapan = read_genes(os.path.join(metapandio, "pangenome.tsv"))
    metapan = metapan.rename(columns = {"gene": "speciesfam",
        "genome": "species"})
    pangenome = pd.merge(speciespans, metapan)
    pangenome = pangenome[["gene", "genome", "orthogroup"]]
    write_tsv(pangenome, os.path.join(args.outfolder, "pangenome.tsv"))

def run_build(args):

    if os.path.isfile(os.path.join(args.outfolder, "hmm_db.h3f")):
        logging.info("existing database detected - moving on")
        return()

    logging.info("creating output subfolders")
    orthogroupsdio = os.path.join(args.outfolder, "orthogroups")
    alignmentsdio = os.path.join(args.outfolder, "alignments")
    profilesdio = os.path.join(args.outfolder, "profiles")
    os.makedirs(orthogroupsdio, exist_ok = True)
    os.makedirs(alignmentsdio, exist_ok = True)
    os.makedirs(profilesdio, exist_ok = True)

    logging.info("gathering sequences of orthogroups")
    pangenome = read_genes(args.pangenome)
    faafins = read_lines(args.faapaths)
    gather_orthogroup_sequences(pangenome, faafins, orthogroupsdio,
        args.min_genomes)
    logging.info(f"gathered sequences for {len(os.listdir(orthogroupsdio))} "
        f"orthogroups occurring in at least {args.min_genomes} genome(s)")

    logging.info("aligning orthogroups")
    orthogroups = [os.path.splitext(file)[0] for file in
        os.listdir(orthogroupsdio)]
    orthogroupfouts = make_paths(orthogroups, orthogroupsdio, ".fasta")
    alifouts = make_paths(orthogroups, alignmentsdio, ".aln")
    run_mafft_parallel(orthogroupfouts, alifouts)

    logging.info("building profile hmms")
    profilefouts = make_paths(orthogroups, profilesdio, ".hmm")
    run_hmmbuild_parallel(alifouts, profilefouts)

    logging.info("pressing profile hmm database")
    run_hmmpress(profilefouts, args.outfolder)

def run_search(args):

    genesfout = os.path.join(args.outfolder, "genes.tsv")
    if os.path.isfile(genesfout):
        logging.info("existing search results detected - moving on")
        return()

    cutoffsfio = os.path.join(args.db, "orthogroups.tsv")

    logging.info("performing hmmsearch")
    queryfins = read_lines(args.qpaths)
    domtblfio = os.path.join(args.outfolder, "hmmer_domtbl.tmp")
    if os.path.isfile(domtblfio):
        logging.info("existing hmmer domtbl detected - skipping hmmsearch")
    else:
        run_hmmsearch(args.db, queryfins, domtblfio)
    hits = read_domtbl(domtblfio)
    # write_tsv(hits, os.path.join(args.outfolder, "hits.tsv"))
    genes_genomes = extract_genes(queryfins)

    if args.trainstrategy == "pan":

        logging.info("traninig profile-specific hmmer cutoffs with the pan "
            "strategy")
        pangenome = read_genes(args.pangenome)
        cutoffs = train_cutoffs_pan(hits, pangenome)
        write_tsv(cutoffs, cutoffsfio)

    elif args.trainstrategy == "core":

        logging.info("traninig profile-specific hmmer cutoffs with the core "
            "strategy")
        cutoffs = train_cutoffs_core(hits, genes_genomes)
        write_tsv(cutoffs, cutoffsfio)

    logging.info("applying hmmer score cutoffs")
    orthogroups = read_orthogroups(cutoffsfio)
    genes = process_scores(hits, orthogroups)
    genes = pd.merge(genes, genes_genomes, how = "left")
    genes = genes[["gene", "genome", "orthogroup"]]
    write_tsv(genes, genesfout)

    logging.info("removing temporary files")
    os.remove(domtblfio)

def run_checkgenomes(args):

    logging.info("checking genomes")
    coregenome = read_genes(args.coregenome)
    genomes = checkgenomes(coregenome)
    write_tsv(genomes, os.path.join(args.outfolder, "genomes.tsv"))

def run_checkgroups(args):

    logging.info("checking core orthogroups")
    coregenome = read_genes(args.coregenome)
    orthogroups = checkgroups(coregenome)
    write_tsv(orthogroups, os.path.join(args.outfolder, "orthogroups.tsv"))

def run_filter(args):

    logging.info("reading pangenome")
    pangenome = read_genes(args.pangenome)
    if not args.genomes is None:
        logging.info("filtering genomes")
        genomes = read_lines(args.genomes)
        pangenome = filter_genomes(pangenome, genomes)
    if not args.orthogroups is None:
        logging.info("filtering orthogroups")
        orthogroups = read_lines(args.orthogroups)
        pangenome = filter_groups(pangenome, orthogroups)
    write_tsv(pangenome, os.path.join(args.outfolder, "pangenome.tsv"))

def run_supermatrix(args):

    sm_aas_fout = os.path.join(args.outfolder, "supermatrix_aas.fasta")
    sm_nucs_fout = os.path.join(args.outfolder, "supermatrix_nucs.fasta")
    seqs_aas_dio = os.path.join(args.outfolder, "seqs_aas")
    seqs_nucs_dio = os.path.join(args.outfolder, "seqs_nucs")
    alis_aas_dio = os.path.join(args.outfolder, "alis_aas")
    alis_nucs_dio = os.path.join(args.outfolder, "alis_nucs")

    logging.info("creating output subfolders")
    os.makedirs(seqs_aas_dio, exist_ok = True)
    os.makedirs(alis_aas_dio, exist_ok = True)

    if os.path.isfile(sm_aas_fout):

        logging.info("existing amino acid supermatrix detected - moving on")

    else:

        logging.info("gathering amino acid sequences of orthogroups")
        coregenome = read_genes(args.coregenome)
        faa_fins = read_lines(args.faapaths)
        gather_orthogroup_sequences(coregenome, faa_fins, seqs_aas_dio)
        logging.info(f"gathered sequences for {len(os.listdir(seqs_aas_dio))} "
            f"core orthogroups")

        logging.info("aligning orthogroups on the amino acid level")
        orthogroups = [os.path.splitext(file)[0] for file in
            os.listdir(seqs_aas_dio)]
        seqs_aas_fios = make_paths(orthogroups, seqs_aas_dio, ".fasta")
        alis_aas_fios = make_paths(orthogroups, alis_aas_dio, ".aln")
        run_mafft_parallel(seqs_aas_fios, alis_aas_fios)

        logging.info("concatenating amino acid alignments")
        construct_supermatrix(coregenome, alis_aas_fios, sm_aas_fout)

    if not args.ffnpaths is None and os.path.isfile(sm_nucs_fout):

        logging.info("existing nucleotide supermatrix detected - moving on")

    if not args.ffnpaths is None and not os.path.isfile(sm_nucs_fout):

        logging.info("creating output subfolders")
        os.makedirs(seqs_nucs_dio, exist_ok = True)
        os.makedirs(alis_nucs_dio, exist_ok = True)

        logging.info("gathering nucleotide sequences of orthogroups")
        coregenome = read_genes(args.coregenome)
        ffn_fins = read_lines(args.ffnpaths)
        gather_orthogroup_sequences(coregenome, ffn_fins, seqs_nucs_dio)

        logging.info("aligning orthogroups on the nucleotide level")
        orthogroups = [os.path.splitext(file)[0] for file in
            os.listdir(seqs_aas_dio)]
        seqs_nucs_fios = make_paths(orthogroups, seqs_nucs_dio, ".fasta")
        alis_aas_fios = make_paths(orthogroups, alis_aas_dio, ".aln")
        alis_nucs_fios = make_paths(orthogroups, alis_nucs_dio, ".aln")
        reverse_align_parallel(seqs_nucs_fios, alis_aas_fios, alis_nucs_fios)

        logging.info("concatenating nucleotide alignments")
        construct_supermatrix(coregenome, alis_nucs_fios, sm_nucs_fout)
