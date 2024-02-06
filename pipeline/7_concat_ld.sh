#!/bin/bash

for CHROM in {1..22}
do
    zcat chroms/${CHROM}/chr${CHROM}.exonic.pol.ld.genematch.annot.gz \
        | awk 'NF==14'
done \
    | awk 'NR==1 {print; next} { if ($0!~"CHR") print }' \
    | gzip -c \
    > genome.exonic.pol.ld.genematch.annot.gz

