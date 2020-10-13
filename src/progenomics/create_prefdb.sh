#!/bin/bash 

TDB="$1"
RES="$2"

# create link to data file which contains a list of all targets that should be aligned
cp "${TDB}.index" "${RES}"

# create dbtype (7)
awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"