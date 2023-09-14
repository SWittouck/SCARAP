#!/usr/bin/env bats

setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    din_concat=testdata/data/06_concatenate
    din_faas=testdata/data/01_pan/faas
    din_faapath=tempdir
    dout=testresults
    mkdir -p $din_faapath
    fin_faas=$din_faapath/faapaths.txt
    ls $din_faas/*.faa.gz | head -3 > $din_faapath/faapaths.txt
}

@test "Can concatenate" {
    scarap concat $fin_faas $din_concat/pangenome.tsv $dout
    [ -s $dout/supermatrix_aas.fasta ]
}

teardown() {
    rm -r $dout
    rm -r $din_faapath 
}