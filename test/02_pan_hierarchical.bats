#!/usr/bin/env bats

setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    din=testdata/data/02_hierarchical
    dout=testresults    # Prepare data
    mkdir -p $dout
    ls $din/faas_hiertest/*.faa.gz > $dout/faapaths_hiertest.txt
}

@test "Can run hierarchical with faapaths" {
    scarap pan $dout/faapaths_hiertest.txt $dout/pan_hierarchical \
    -s $din/species_hiertest.tsv -t 16  
    [ -s $dout/pan_hierarchical/pangenome.tsv ]
}

@test "Can run hierarchical with faafolder" {
    scarap pan $din/faas_hiertest $dout/pan_hierarchical \
    -s $din/species_hiertest.tsv -t 16
    [ -s $dout/pan_hierarchical/pangenome.tsv ]
}

teardown() {
    rm -r $dout 
}