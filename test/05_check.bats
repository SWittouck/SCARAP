#!/usr/bin/env bats

setup() {
    din=testdata/data/03_build
    dout=testresults
    dgenomes=$dout/checkgenomes
    dgroups=$dout/checkgroups
    mkdir -p $dout
}

@test "Can run check genomes" {
    scarap checkgenomes $din/pangenome.tsv $dgenomes
    [ -s $dgenomes/genomes.tsv ]
}

@test "Can run check groups" {
    scarap checkgroups $din/pangenome.tsv $dgroups
    [ -s $dgroups/orthogroups.tsv ]
}

@test "Can filter" {
    scarap checkgroups $din/pangenome.tsv $dgroups
    awk '$2 == 1.0 {print $1}' $dgroups/orthogroups.tsv \
     > $dout/orthogroups_core.txt

    cat $dout/orthogroups_core.txt | head -n 10 > $dout/orthogroups_core_sub.txt

    scarap filter -c -o $dout/orthogroups_core_sub.txt \
     $din/pangenome.tsv $dout
    [ -s $dout/pangenome.tsv ]
}

teardown() {
    rm -r $dout
}
