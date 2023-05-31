#!/usr/bin/env bats

setup() {
    load 'test_helper/bats-support/load'
    load 'test_helper/bats-assert/load'
    din=testdata/data/01_pan/faas
    dout=testresults
    # Prepare data
    mkdir -p $dout
    ls $din/*.faa.gz | head -3 > $dout/faapaths.txt
}
    
@test "Can run with faapaths" {
    scarap pan $dout/faapaths.txt $dout/pan -t 16 # with faapaths
}

@test "Can run with faafolder" {
    scarap pan $din $dout/pan -t 16 # with faafolder
}

@test "Running with redundant genes should give error" {
    ls $din/{extra,redundant}/*.faa.gz >> $dout/faapaths.txt
    redundant_run() {
        scarap pan $dout/faapaths.txt $dout/pan -t 16
    }
    run redundant_run
    assert_failure
}

teardown() {
    rm -r $dout 
}
