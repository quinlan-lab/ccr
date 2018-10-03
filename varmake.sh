original=$1
fileName="${original##*/}"
fileExt=${fileName#*.}
FILE=${fileName%*.$fileExt}
if [ ! -s data/${FILE}_dec.vcf ]; then
    vt decompose $original -o data/${FILE}_dec.vcf -s 
fi

if [ ! -s data/${FILE}_vt.vcf ]; then
    vt normalize data/${FILE}_dec.vcf -o data/${FILE}_vt.vcf -r data/grch37.fa
fi

perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i data/${FILE}_vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o data/${FILE}-vep-vt.vcf --vcf --fields ALLELE,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM,cDNA_position --offline --fork 12 --force_overwrite

cat <(grep "^#" data/${FILE}-vep-vt.vcf) <(grep -v "^#" data/${FILE}-vep-vt.vcf | sort -k1,1 -k2,2n) | bgzip -c > data/${FILE}-vep-vt.vcf.gz; tabix -f data/${FILE}-vep-vt.vcf.gz
