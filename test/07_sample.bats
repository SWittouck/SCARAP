#!/usr/bin/env bats

setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    din_faas=testdata/data/01_pan/faas
    din_concat=testdata/data/07_sample
    din_faapath=tempdir
    dout=testresults
    mkdir -p $din_faapath
    fin_faas=$din_faapath/faapaths.txt
    ls $din_faas/*.faa.gz | head -3 > $din_faapath/faapaths.txt
    ls $din_faas/extra/*.faa.gz >> $din_faapath/faapaths.txt
}

@test "Can sample results" {
    scarap sample $fin_faas $din_concat/pangenome.tsv $dout \
      --max-genomes 5 --method mean50
    #[ -s $dout/clusters.tsv ]
    [ "$(wc -l < $dout/seeds.txt)" -eq 5 ]
}

teardown() {
    rm -r $dout
    rm -r $din_faapath
}