#!/usr/bin/env bats

setup() {
    din=testdata/data/01_pan/faas
    din_build=testdata/data/03_build
    dout=testresults
    mkdir -p $dout
    ls $din/*.faa.gz | head -3 > $dout/faapaths.txt
}

@test "Build" {
    scarap build $dout/faapaths.txt $din_build/pangenome.tsv $dout/build \
      -p 90 -f 95 -m 100
    [ -s $dout/build/cutoffs.csv ]
}

teardown() {
    rm -r $dout
}
