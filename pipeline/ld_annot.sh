#!/bin/bash

CHROM=$1

CHROM_DIR=chroms/${CHROM}
LD_FILE=${CHROM_DIR}/chr${CHROM}.exonic.pol.ld.genematch.gz

zcat ${LD_FILE} \
    | /net/home/colby/projects/epistasis/code/extract_LD_class.py \
        -a ${CHROM_DIR}/impact.txt.gz \
    | groupBy -header -i stdin -g 3,7 -full -c 11,12 -o distinct,distinct \
    | awk 'NR==1 {print; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$15,$16,$13,$14 }' OFS="\t" \
    | gzip -c \
    > ${CHROM_DIR}/chr${CHROM}.exonic.pol.ld.genematch.annot.gz

