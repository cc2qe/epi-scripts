#!/bin/bash

CHROM=$1

CHROM_DIR=chroms/${CHROM}
VCF_PREFIX=chroms/${CHROM}/chr${CHROM}.af.cadd.snpEff.polarized
PLINK_PREFIX=chroms/${CHROM}/chr${CHROM}.exonic.pol
OUT_PREFIX=chroms/${CHROM}

# # make impact file, tab delimited file of functional impact by variant
# cat ${VCF_PREFIX}.noGT.vcf \
#     | vawk --header 'I$AF>0' \
#     | vawk '{
#         split(I$ANN,x,",");
#         for(idx in x) { split(x[idx],y,"|");
#             print $1,$2,$3,y[4],y[2],I$CADDPHRED }
#         }' \
#     | zapdups -u \
#     | bgzip -c \
#     > ${OUT_PREFIX}/impact.txt.gz &&

# # From the LD file, extract only variant pairs affecting the same gene
# zcat ${PLINK_PREFIX}.ld.gz | tr -s ' ' '\t' | sed 's/^\t//g' \
#     | zjoin -p "CHR" -a stdin -b <(zcat ${OUT_PREFIX}/impact.txt.gz) -1 3 -2 3 \
#     | cut -f -10,15 \
#     | zjoin -p "CHR" -a stdin -b <(zcat ${OUT_PREFIX}/impact.txt.gz) -1 7 -2 3 \
#     | cut -f -11,15 \
#     | awk 'NR==1 { OFS="\t"; print $0,"GENE_A","GENE_B"; next} { if ($(NF-1)==$NF) print }' \
#     | zapdups \
#     | bgzip -c \
#     > ${PLINK_PREFIX}.ld.genematch.gz

# join the LD file with annotations (CADD score, gene) and output data
zcat ${PLINK_PREFIX}.ld.genematch.gz \
    | /net/home/colby/projects/epistasis/code/extract_LD_class.py \
        -a ${CHROM_DIR}/impact.txt.gz \
    | /net/home/colby/projects/epistasis/code/filter_impact.py \
    | gzip -c \
    > ${CHROM_DIR}/chr${CHROM}.exonic.pol.ld.genematch.annot.csq.gz

