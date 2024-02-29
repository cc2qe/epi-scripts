#!/bin/bash

for CHROM in {1..22}
do
    zcat chroms/${CHROM}/chr${CHROM}.exonic.pol.ld.genematch.annot.csq.gz \
        | awk 'NF==16'
done \
    | awk 'NR==1 {print; next} { if ($0!~"CHR") print }' \
    | gzip -c \
    > genome.exonic.pol.ld.genematch.annot.csq.gz

