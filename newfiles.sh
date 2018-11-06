#!/bin/bash

mkdir -p data; mkdir -p results

if [ ! -s data/Homo_sapiens.GRCh37.85.gtf.gz ]; then
    wget -P data ftp://ftp.ensembl.org/pub/grch37/release-85/gtf/homo_sapiens/Homo_sapiens.GRCh37.85.gtf.gz 
    cd data; gunzip Homo_sapiens.GRCh37.85.gtf.gz; sort -k1,1 -k4,4n Homo_sapiens.GRCh37.85.gtf | bgzip -c > Homo_sapiens.GRCh37.85.gtf.gz; tabix Homo_sapiens.GRCh37.85.gtf.gz; cd -
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
#gnomAD v2.1 VCF
if [ `ls -1 data/gnomad.exomes.r2.1.sites.vcf.bgz* 2>/dev/null | wc -l ` -eq 0 ]; then
    wget -P data https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz
    wget -P data https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.bgz.tbi
fi
#gnomAD v2.1 coverage is same as 2.0.2
if [ `ls -1 data/nomad.exomes.coverage.summary.tsv.bgz | wc -l ` -eq 0 ]; then
    wget -P data https://storage.googleapis.com/gnomad-public/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz
    tabix -b 2 -e 2 -c c data/gnomad.exomes.coverage.summary.tsv.bgz # gnomAD 2.1 has an uncommented first line but no line should start with c other than that
fi

#gnomAD
#v2.1
if [ ! -s data/gnomad2.1-vep-vt.vcf.gz ]; then
    bash varmake.sh data/gnomad.exomes.r2.1.sites.vcf.bgz 2.1
fi
