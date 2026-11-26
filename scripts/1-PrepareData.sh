#!/bin/bash
mkdir -p working_sequences
for outlet in data/sequences_by_outlet/outlet_*; do
    name=$(basename $outlet)
    for file in $outlet/*.fasta; do
        sed "s/^>/>${name}|/" $file > working_sequences/$(basename $file)
    done
done