#!/usr/bin/env bats

setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    din=testdata/data/01_pan/faas
    dout=testresults
    # Prepare data
    mkdir -p $dout
    ls $din/*.faa.gz | head -5 > $dout/faapaths_coretest.txt
}

@test "Can run core module" {
    scarap core $dout/faapaths_coretest.txt $dout/core -t 16
}

teardown() {
    rm -r $dout 
}