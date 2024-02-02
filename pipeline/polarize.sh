#!/bin/bash

CHROM=$1

VCF_PREFIX=chroms/${CHROM}/chr${CHROM}.af.cadd.snpEff
FASTA=/net/home/colby/projects/epistasis/annotations/GRCh37/ancestral/homo_sapiens_ancestor_GRCh37_e71/chr_prefix/chr${CHROM}.fa

zcat ${VCF_PREFIX}.vcf.gz \
    | vawk --header 'I$AA!="X"' \
    | bcftools norm \
        -m-any \
        -Ov \
        --check-ref ws \
        -f $FASTA \
    | bcftools plugin fill-tags -- -t AF \
    | bcftools view --min-ac 1 \
    | bgzip -c \
    > ${VCF_PREFIX}.polarized.vcf.gz &&

tabix -p vcf -f ${VCF_PREFIX}.polarized.vcf.gz

