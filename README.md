A map of constrained coding regions (CCRs) in the human genome.
-------------
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=for-the-badge)](https://quinlan-lab.github.io/ccr/)

Click on the badge above to read the docs and learn more about how to run the model.

Go to INSTALL to see all packages and software necessary to run the model.
 
--------------

If you use this model in any way, please cite the paper:

[Havrilla, J.M., Pedersen, B.S., Layer, R.M. & Quinlan, A.R. A map of constrained coding regions in the human genome. Nature Genetics (2018). doi:10.1038/s41588-018-0294-6](https://www.nature.com/articles/s41588-018-0294-6)

This repository is linked to our Nature Genetics manuscript on constrained coding regions in the human genome. If you would like to view CCRs throughout the genome or download the model in its most current form, go to the [CCR Browser](https://rebrand.ly/ccrregions).  The version used in the paper, which is the version currently available for download utilizes the hg19/GRCh37 reference genome.

The constrained coding regions model (CCR) uses the power of the Genome Aggregation Database (gnomAD) to discover areas of the human exome under potentially purifiying selection.  Using protein-altering variation from across 123,136 ostensibly healthy individuals' exomes as the boundaries of each region of constraint, we create regions completely devoid of any protein-coding variation. These regions' size is variable and dependent upon where the variants lie. The model then ranks these regions using a linear regression of coverage-weighted length against the CpG content of those regions.

The most constrained regions (&ge;90th percentile) have been shown to be extremely enriched for pathogenic variation in ClinVar, _de novo_ mutations in patients with severe developmental disorders, and critical Pfam domains exome-wide.  Even more exciting, 72% of genes harboring a CCR in the 99th percentile or higher have no known pathogenic variants.  There is great opportunity for discovery of function in these understudied genes as well as their role in disease phenotypes or potentially even in embryonic lethality when altered.

After downloading the CCR files at the [CCR Browser](https://rebrand.ly/ccrregions), note the meaning of the columns in the file is as described at the [browser repo](https://github.com/quinlan-lab/ccrhtml):

#### BED file columns
Column              | Description |
--------            | ----------- |
chrom               | Chromosome ID  
start               | Start coordinate (may be part of a multi-exon CCR)
end                 | End coordinate (may be part of a multi-exon CCR)
ccr_pct             | CCR percentile.  0 represents ExAC variants and is total non-constraint.  100 represents complete constraint, the highest constrained region in the model. 
gene                | HGNC gene name.
ranges              | The range of coordinates that represent the CCR.  For multi-exon spanning CCRs, this will be a comma-separated list of ranges.
varflag             | VARTRUE = 0th percentile CCR, and thus an ExAC variant coordinate (or several ExAC deletions merged into one CCR).  VARFALSE = Anything that is not a 0th percentile CCR. 
syn_density         | A calculation of the synonymous variant density of the CCR region.  Used variants that were SNPs and did not change amino acids or stop/start codons.  Allowed multiple alleles at same bp.
cpg                 | CpG dinucleotide density of the whole CCR region. 
cov_score           | The score of length scaled by coverage proportion at 10x for each base pair.  
resid               | Raw residual value from the linear regression model. 
resid_pctile        | Raw residual percentile, not weighted by proportion of exome represented.
unique_key          | A unique key ID for each CCR.
