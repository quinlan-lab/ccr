#!/bin/bash

mkdir -p data; mkdir -p results

if [ ! -s data/Homo_sapiens.GRCh37.75.gtf.gz ]; then
    wget -P data ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
    cd data; tabix Homo_sapiens.GRCh37.75.gtf.gz; cd -
fi

if [ ! -s data/grch37.fa ]; then
    wget -P data https://s3.us-east-2.amazonaws.com/pathoscore-data/fastas/grch37.fa
    wget -P data https://s3.us-east-2.amazonaws.com/pathoscore-data/fastas/grch37.fa.fai
fi

if [ ! -s data/segmental.bed.gz ]; then
    wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz # segmental duplications
    # grab and make self-chain and segdup tracks
    zcat < genomicSuperDups.txt.gz | cut -f 2-4 | sed 's/^chr//g' | sort -k1,1 -k2,2n | bedtools merge | bgzip -c > data/segmental.bed.gz; cd data; tabix -f segmental.bed.gz; cd -
    rm genomicSuperDups.txt.gz
fi

if [ ! -s data/self-chains.id90.bed.gz ]; then
    python get-chain.py | sed 's/^chr//g' | sort -k1,1 -k2,2n | bedtools merge | bgzip -c > data/self-chains.id90.bed.gz; cd data; tabix -f self-chains.id90.bed.gz; cd -
fi

# get gnomAD VCF, coverage files
#gnomAD v2.0.1 VCF
if [ `ls -1 data/gnomad.exomes.r2.0.1.sites.vcf.gz* 2>/dev/null | wc -l ` -eq 0 ]; then
    wget -P data https://storage.googleapis.com/gnomad-public/release/2.0.1/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz
    wget -P data https://storage.googleapis.com/gnomad-public/release/2.0.1/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz.tbi
fi
#gnomAD v2.0.2 VCF
#if [ `ls -1 data/gnomad.exomes.r2.0.2.sites.vcf.bgz* 2>/dev/null | wc -l ` -eq 0 ]; then
#    wget -P data https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz
#    wget -P data https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz.tbi
#fi
#gnomAD v2.0.2 coverage (same as v2.0.1 coverage)
if [ `ls -1 data/gnomad.exomes.r2.0.2.chr* 2>/dev/null | wc -l ` -eq 0 ]; then
    for i in {1..22} X Y; do wget -P data https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/exomes/gnomad.exomes.r2.0.2.chr$i.coverage.txt.gz; done
    for i in {1..22} X Y; do wget -P data https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/exomes/gnomad.exomes.r2.0.2.chr$i.coverage.txt.gz.tbi; done
fi

# if and only if you want ExAC v1 for some reason:
if [ `ls -1 data/ExAC.r1.sites.vep.vcf.gz* 2>/dev/null | wc -l ` -eq 0 ]; then
    wget -P data ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz
    wget -P data ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz.tbi
fi
if [ `ls -1 data/Panel.chr* 2>/dev/null | wc -l ` -eq 0 ]; then
    for i in {1..22} X Y; do wget -P data ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/coverage/Panel.chr$i.coverage.txt.gz; done
    for i in {1..22} X Y; do wget -P data ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/coverage/Panel.chr$i.coverage.txt.gz.tbi; done
fi

#gnomAD
#v2.0.1
if [ ! -s data/gnomad-vep-vt.vcf.gz ]; then
    bash varmake.sh data/gnomad.exomes.r2.0.1.sites.vcf.gz
fi
#ExAC
#v1 (final)
if [ ! -s data/ExAC-vep-vt.vcf.gz ]; then
    bash varmake.sh data/ExAC.r1.sites.vep.vcf.gz
fi
