A map of constrained coding regions (CCRs) in the human genome.
-------------
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=for-the-badge)](https://quinlan-lab.github.io/ccr/)

Click on the badge above to read the docs and learn more about how to run the model.

Go to INSTALL to see all packages and software necessary to run the model.
 
--------------

Overview
========
This repository is linked to our [manuscript](https://www.nature.com/articles/s41588-018-0294-6) describing constrained coding regions in the human genome. If you would like to view CCRs throughout the genome or download the model in its most current form, go to the [CCR Browser](https://www.rebrand.ly/ccrregions).  The version used in the paper, which is the version currently available for download, utilizes the hg19/GRCh37 reference genome.

The constrained coding regions model (CCR) uses the Genome Aggregation Database (gnomAD, version 2.0.1 in the paper) to reveal regions of protein coding genes that are likely to be under potentially purifiying selection. We used protein-altering variation from across 123,136 ostensibly healthy individuals' exomes to reveal coding regions that are completely devoid of any protein-coding variation. We infer such coding regions to be constrained; the higher the constraint percentile, the more constrained we predict the region to be. 

The most constrained regions (&ge;90th percentile, and especially at or above the &ge;99th percentile) have been shown to be extremely enriched for pathogenic variation in ClinVar, _de novo_ dominant mutations in patients with severe developmental disorders, and critical Pfam domains exome-wide.  Even more exciting, 72% of genes harboring a CCR in the 99th percentile or higher have no known pathogenic variants.  There is great opportunity for discovery of function in these understudied genes as well as their role in disease phenotypes or potentially in embryonic lethality when altered.

Citation
========
If you use this model in any way, please cite the paper:

[Havrilla, J.M., Pedersen, B.S., Layer, R.M. & Quinlan, A.R. A map of constrained coding regions in the human genome. Nature Genetics (2018). doi:10.1038/s41588-018-0294-6](https://www.nature.com/articles/s41588-018-0294-6)

CCR BED Files
=============
- [Autosomal CCR regions in extended BED format (gnomAD v2.0.1)](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz)
- [X chromosome CCR regions in extended BED format (gnomAD v2.0.1)](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz)
- [Autosomal CCR regions at 90th percentile and higher in extended BED format (gnomAD v2.0.1)](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.90orhigher.v2.20180420.bed.gz)
- [X chromosome CCR regions at 90th percentile and higher in extended BED format (gnomAD v2.0.1)](https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.90orhigher.v2.20180420.bed.gz)

Each column in the above CCR BED files is described below:

#### BED file columns
Column              | Description |
--------            | ----------- |
chrom               | Chromosome ID  
start               | Start coordinate (0-based, may be part of a multi-exon CCR)
end                 | End coordinate (1-based, may be part of a multi-exon CCR)
ccr_pct             | CCR percentile.  0 represents gnomAD variants and is total non-constraint.  100 represents complete constraint, the highest constrained region in the model. 
gene                | HGNC gene name.
ranges              | The range of coordinates that represent the CCR.  For multi-exon spanning CCRs, this will be a comma-separated list of ranges.
varflag             | VARTRUE = 0th percentile CCR, and thus an ExAC variant coordinate (or several ExAC deletions merged into one CCR).  VARFALSE = Anything that is not a 0th percentile CCR. 
syn_density         | A calculation of the synonymous variant density of the CCR region.  Used variants that were SNPs and did not change amino acids or stop/start codons.  Allowed multiple alleles at same bp.
cpg                 | CpG dinucleotide density of the whole CCR region. 
cov_score           | The score of length scaled by coverage proportion at 10x for each base pair.  
resid               | Raw residual value from the linear regression model. 
resid_pctile        | Raw residual percentile, not weighted by proportion of exome represented.
unique_key          | A unique key ID for each CCR.
