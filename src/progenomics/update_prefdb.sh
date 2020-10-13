#!/bin/bash 

QDB="$1"
TDB="$2"
RES="$3"

# create new index repeatedly pointing to same entry
INDEX_SIZE="$(echo $(wc -c < "${TDB}.index"))"
awk -v size=$INDEX_SIZE '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"