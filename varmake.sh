set -eo pipefail
original=$1
fileName="${original##*/}"
fileExt=${fileName#*.}
FILE=${fileName%*.$fileExt}
OTHERNAME=$2
if [ ! -s data/${FILE}${OTHERNAME}_dec.vcf ]; then
    vt decompose $original -o data/${FILE}${OTHERNAME}_dec.vcf -s 
fi

if [ ! -s data/${FILE}${OTHERNAME}_vt.vcf ]; then
    vt normalize data/${FILE}${OTHERNAME}_dec.vcf -o data/${FILE}${OTHERNAME}_vt.vcf -r data/grch37.fa
fi

rm data/${FILE}${OTHERNAME}_dec.vcf # files can get QUITE large

# add variant_effect_predictor to path, contains line at the top #!/usr/bin/env perl to allow for execution
variant_effect_predictor.pl -i data/${FILE}${OTHERNAME}_vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o data/${FILE}${OTHERNAME}-vep-vt.vcf --vcf --fields ALLELE,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM,cDNA_position --offline --fork 12 --force_overwrite

rm data/${FILE}${OTHERNAME}_vt.vcf

cat <(grep "^#" data/${FILE}${OTHERNAME}-vep-vt.vcf) <(grep -v "^#" data/${FILE}${OTHERNAME}-vep-vt.vcf | sort -k1,1 -k2,2n) | bgzip -c > data/${FILE}${OTHERNAME}-vep-vt.vcf.gz; tabix -f data/${FILE}${OTHERNAME}-vep-vt.vcf.gz

rm data/${FILE}${OTHERNAME}-vep-vt.vcf
