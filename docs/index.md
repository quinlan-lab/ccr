CCR Scripts
==============

This repo will allow you to create the model as we did.  If you want to run the model as in the paper, go to [Examples --> gnomAD](examples/gnomad.md) above.  New versions are under [Examples --> New Versions](examples/updates.md).

exac-regions.py
---

This script generates the regions devoid of genetic variation.  It utilizes scripts in `utils.py` to do so.  Regions are actually calculated using the functions in `utils.py` which will be explained in the next section.

+ `-n` will prevent singletons from being used in the model.  If you do not trust variants with an allele count of 1 in your VCF file, this may be an option you are interested in.  Alternatively, if you are trying to create a recessive model you may not want these either.  In the future, this may be replaced by an allele count option (e.g., AC=2,Hom=1).
+ `-w` is a necessary option as of now.  It is the version we used in the manuscript where variants are separate from CCRs.  Originally, the downstream variant was included in the region definition.  Now, that is no longer the case.  Variants are merged together if overlapping INDELs and treated as one large region of non-constraint.  If you want to incorporate variants into the regions, be aware the code to do so is outdated, and likely will be bugged.  We recommend that you use this version.
+ `-x` tells the code where your vt decomposed, normalized and VEP annotated variant file is.  You definitely want the VEP fields used in varmake.sh.  You can just run `bash varmake.sh` on whatever your VCF file is and look for it with the partial name of the file added to "-vep-vt.vcf.gz".
+ `-e` tells the code where your Ensembl formatted exons file is.  You want this in GTF gzipped format.
+ `-c` tells the code where the coverage files are.  We recommend you create a data folder and store them in there.  There should be one for each chromosome.  The actual input is formatted like `data/gnomad.exomes.r2.0.2.chr{chrom}.coverage.txt.gz` so it will pull the file for each chromosome internally using the literal string `{chrom}` it just needs to know where to look.
+ `-d` allows you to pass a coverage depth into the model.  The default is 10 for 10x coverage, as in the manuscript.
+ `-l` allows you to provide a proportional cutoff for coverage into the model.  As in the manuscript the default is 0.5, or 50% of individuals in your dataset that have at least the coverage depth you provided in the -d option.
+ `-f` will use a FASTA file of your choice to calculate CpG density.  The default is the GRCh37/hg19 FASTA downloaded from our s3 server by getfiles.sh.  You do want to make sure it does NOT start with chr.
+ `-s` tells the code where your exclude files are.  There are no defaults, but in the manuscript we used gzipped BED files of segmental duplications and self-chains of 90% identity as short reads are unreliable in those regions of the genome.  You can use any exclusionary BED.gz files here.

Make sure you save the output to a file, by default it is printed to stdout.
`exac-regions.py` runs using the multiprocessing package, so debugging may be difficult for those less experienced.

resid-plot.py
---

This script is used to take the output of `exac-regions.py` and creates the model for the regions.  It also calculates synonymous density. Synonymous density calculation currently takes a very long time, so it may be in the interest of the user to forgo the calculation if iterating through several different model versions.  In the future it should be trivial to add more aspects to the model or switch from linear regression to a logistic regression model or some other ML model.  This script does use one function from `utils.py`, which is `issynonymous` to aid with determination of whether variants are synonymous or not.

+ `-c` triggers CpG density input into the model.  This is used in the manuscript's model.  It will always be calculated no matter what because it is already calculated in the output of `exac-regions.py`.
+ `-s` tells the model to calculate synonymous variant density, which you may not want as the model takes a lot longer to do so. If you do this without the -c option for CpG, it will calculate a model based on coverage-scaled length and synonymous density instead of CpG content.  If you want to calculate the density, but not incorporate it into the model, use -r with -s.
+ `-f` tells the code where your input file is.  This is the output from exac-regions.py.
+ `-n` tells the script not to use singletons for synonymous density calculation.  You don't need this option unless you invoke -s.  Also, if you feel confident in your variant calls, there is no reason to use this option.
+ `-w` is a necessary option as of now.  It assumes there is a "varflag" field.  It is the version we used in the manuscript where variants are separate from true CCRs.  Originally, the downstream variant was included in the region definition, but we no longer do so.  If you want to incorporate variants into the regions it must have been done already at the `exac-regions.py` step, and you should be very aware the code to do so is outdated, and likely will be bugged.  We recommend that you use this flag.
+ `-p` allows you to only run the model on one chromosome or a few of them.  We used this option to run a CCR model based solely on the X chromosome, as in the manuscript. If using multiple chromosomes, e.g. X Y, separate inputs with spaces.
+ `-x` will let you exclude certain chromosomes from your analysis.  The default is Y but as in the manuscript, we recommend using this to get autosomal CCRs and excluding X and Y.  But if you want everything just run with a blank.  Separate inputs with spaces.
+ `-q` will create a version of the full model where the X chromosome is waited differently based on the allele count max in areas not in pseudoautosomal regions.  Only relevant if using the autosomes and at least the X chromosome.  It could be used for an X-only model, but would have little point.
+ `-r` removes synonymous density from the model and just uses the calculation invoked by -s as a field in the BED file.  We did this for our manuscript, where the model uses CpG content only.  No need to use this option without -s.

This script will produce its output in stdout.  Capture it in a file, and make sure that the file is sorted by percentile (the last column) and the header as the first line before passing the file into `weightpercentile.py`.

weightpercentile.py
---

This only takes one input, the file from `resid-plot.py` sorted by percentile, keeping the header.

This script recalculates the percentile scale using the percentiles from `resid-plot.py`.  It starts at 100 (the most constrained region) and scales it down proportionally for each subsequent region by the fraction of the length of protein-coding exome already covered by the preceding region.  It is highly recommended that you use this script to create the final percentiles as the percentiles created by `resid-plot.py` are very imbalanced.  The output of this file contains the final CCRs for your input variant set and exons of choice.  This is the end of the pipeline.

utils.py
---

This script contains all the functions necessary for creating the regions as well as a few other useful functions.  There are several doctests for each function that check if they work as expected.

+ `overlaps` is used to help check if exclusionary ranges overlap generated regions in the function `split_ranges`.
+ `split_ranges` takes the output of the functions `get_ranges` and `read_exons` and splits generated regions on exclude files provided by the user.
+ `get_ranges` is the star of the program.  It takes variant start coordinate data, end coordinate data (if INDEL, it is longer), and the exon\_starts and exon\_ends for all transcripts for a gene acquired from `read_exons`, as well as the relevant chromosome for multiprocessing purposes.  It also takes in the last coordinate used as either a variant or the end of the region.  It can create long chains of INDELs as one region of non-constraint and regions of true constraint with no variation.  This is also how multi-exon regions are generated.  If there are several VARFALSE flagged regions in the output, it is a multi-exon region.  There are numerous doctests for this function for the many corner cases.  If -w is used in `exac-regions.py` this is the function used.  
+ `get_ranges_w_variant` is the old function for generating regions including the variant on the edge, but is outdated.  Do not use this function without checking it first.
+ `path` just expands the variables in the path.
+ `floatfmt` is mostly for formatting the float value of CpG density into a shorter string.
+ `read_coverage` uses the relevant chromosome's coverage file, for the depth provided by the user, and creates an array of the proportion of individuals at that depth into a numpy array for the appropriate coordinates.
+ `read_exons` uses the Ensembl GTF file and the coverage cutoff proportion as well as the exclude files provided by the user to get all the exon starts and ends flattened across all transcripts in an IntervalSet object for use by `get_ranges`.  It also outputs the Interval objects to split the ranges on in `split_ranges` for later.
+ `get_cdna_start_end` gets cDNA coordinate information from VEP annotations on the variant entries.  However, this is no longer used for any particular purpose in `exac-regions.py`.
+ `isfunctional` is the function used to determine if any variant is protein-altering in the variant input set which is what is used to break regions in `exac-regions.py`.
+ `ismissense` works similarly to `isfunctional` but only checks if a variant alters one amino acid.
+ `issynonymous` is actually used in `resid-plot.py` to check if a variant alters a codon but does not change the protein product.  It is used for synonymous density calculation.
+ `cg_content` is used to calculate CpG dinucleotide density for `exac-regions.py` to be passed into the model generated in `resid-plot.py`.
