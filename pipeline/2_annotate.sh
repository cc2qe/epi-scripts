#!/bin/bash

CHROM=$1

VCF_PREFIX=chroms/${CHROM}/chr${CHROM}
FASTA=/net/home/colby/projects/epistasis/annotations/GRCh37/ancestral/homo_sapiens_ancestor_GRCh37_e71/chr_prefix/chr${CHROM}.fa

# Note: must use my forked version of vcfdo (https://github.com/cc2qe/vcfdo)
# since the main release throws errors due to numpy deprecations
zcat ${VCF_PREFIX}.vcf.gz \
    | bcftools plugin fill-tags -- -t AF \
    | vawk --header 'I$AF>0' \
    | vcfdo polarize -v \
    -f $FASTA \
    | vcfanno -p 1 cadd.toml /dev/stdin \
        | bgzip -c \
        > ${VCF_PREFIX}.af.cadd.vcf.gz &&

tabix -p vcf -f ${VCF_PREFIX}.af.cadd.vcf.gz

