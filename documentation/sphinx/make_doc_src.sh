#!/bin/bash

cd ../../
rm -rf doc_src
for f in $(find src/ -type f -print); do
    newf="doc_${f}"
    mkdir -p $(dirname $newf)
    LANG=C sed -e 's%LIB_LOCAL%%g' -e 's%LIB_PUBLIC%%g' $f > $newf
done
