A map of constrained coding regions (CCRs) in the human genome.
==============

If you use this model in any way, please cite the paper:

[Havrilla, J.M., Pedersen, B.S., Layer, R.M. & Quinlan, A.R. A map of constrained coding regions in the human genome. Nature Genetics (2018). doi:10.1038/s41588-018-0294-6](https://www.nature.com/articles/s41588-018-0294-6)

This repository is linked to our Nature Genetics manuscript on constrained coding regions in the human genome, which is still in press.  Currently the paper is on [bioRxiv](https://www.biorxiv.org/content/early/2017/11/22/220814).  If you would like to view CCRs throughout the genome or download the model in its most current form, go to the [CCR Browser](https://rebrand.ly/ccrregions).

The constrained coding regions model (CCR) uses the power of the Genome Aggregation Database (gnomAD) to discover areas of the human exome under potentially purifiying selection.  Using protein-altering variation from across 123,136 ostensibly healthy individuals' exomes as the boundaries of each region of constraint, we create regions completely devoid of any protein-coding variation. These regions' size is variable and dependent upon where the variants lie. The model then ranks these regions using a linear regression of coverage-weighted length against the CpG content of those regions.

The most constrained regions (&ge;90th percentile) have been shown to be extremely enriched for pathogenic variation in ClinVar, _de novo_ mutations in patients with severe developmental disorders, and critical Pfam domains exome-wide.  Even more exciting, 72% of genes harboring a CCR in the 99th percentile or higher have no known pathogenic variants.  There is great opportunity for discovery of function in these understudied genes as well as their role in disease phenotypes or potentially even in embryonic lethality when altered.


Go to INSTALL to see all packages and software necessary to run the model.

Click on the badge below to read the docs and learn more about how to run the model.

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=for-the-badge)](https://quinlan-lab.github.io/ccr/)
