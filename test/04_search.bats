#!/usr/bin/env bats

setup() {
    din=testdata/data
    dout=testresults
    mkdir -p $dout
    ls $din/01_pan/faas/*.faa.gz | head -3 > $dout/faapaths.txt
    cat $dout/faapaths.txt | head -4 | tail -1 > $dout/qpath.txt
}

@test "Can run a search" {
    scarap search $dout/faapaths.txt $din/04_search $dout/search
    [ -s $dout/search/genes.tsv ]
}

teardown() {
    rm -r $dout
}
