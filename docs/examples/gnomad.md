Creating CCRs as in the paper
==============

This repo will allow you to create the model as we did.  This is a great way to get a quickstart on how the model is run.  Run these two simple bash commands and it will run all the python scripts and bash code you need to grab all files and create the CCR model files.

First, run:
```
bash getfiles.sh
```
to grab everything you need.  This will also run `bash varmake.sh` and create vt decomposed and normalized files for the gnomAD VCF annotated with appropriate VEP fields.  It may take a while.  This also uses `get-chain.py` to create the self-chain file and creates the `results` and `data` folders in this folder.  The output of this script goes entirely into the `data` folder.  The segmental duplication and self-chain files `segmental.bed.gz` and `self-chains.id90.bed.gz` are placed into this folder, the GRCh37 FASTA file `grch37.fa`, the Ensembl exons GTF file `Homo_sapiens.GRCh37.75.gtf.gz`, the first version gnomAD VCF `gnomad.exomes.r2.0.1.sites.vcf.gz`, and lastly the coverage files for gnomAD labeled as `gnomad.exomes.r2.0.2.chr$i.coverage.txt.gz` with `$i` replaced with each chromosome name.  The `varmake.sh` script takes the gnomAD VCF file and outputs a decomposed, normalized, VEP-annotated file called `gnomad-vep-vt.vcf.gz` in the same folder.

Then run:
```
bash sbatch.sh
```
to run `regions.sh` for all versions of the model found in the manuscript.  It contains examples of how the code is run to obtain the versions used in the manuscript.  It creates output using `regions.sh` in three folders.  `results/gnomAD10.x5syn` is the full autosome-based model used in the paper.  `results/ExACv1syn` is the version for ExAC v1 which only uses 60,706 exomes which is also referenced in the paper's supplementary figures.  Lastly, `results/Xchromonly` is the X chromosome only version of CCR, dubbed X-CCR, which is also referenced in the paper's supplement.

In each results folder, `weightedresiduals-cpg-synonymous-novariant.txt` is the final file containing all CCRs with appropriately adjusted percentiles.  The `resids-cpg-synonymous-novariant.txt` is the file before size-weighted percentile adjustment, and `exac-regions-novariant.txt` is the file of all the regions generated with no model applied to ranking them at all.

Lastly, to create the bgzipped, sorted and tabixed bedfiles, as seen on the CCR browser run:
```
bash makebeds.sh
```
which will create the appropriately versioned bgzipped bedfiles in the corresponding results folders.

### regions.sh

So to run the full autosomal model as we did in the manuscript, we used the following command:
```
bash regions.sh -c -s -w -v gnomAD10x.5syn -d 10 -d 0.5 -x data/segmental.bed.gz -x data/self-chains.id90.bed.gz -q X -q Y -g -u
```
So the options for running `regions.sh` are as follows:

+ -v allows you to choose the name of the results folder you store the results in. This will be in the directory results inside of this folder. It can be a name, the date you ran it, whatever you like.
+ -x allows you to pass files of choice (you don't have to pass any) to exclude regions from the genome when generating your model.  We recommend excluding segmental duplications and self-chains of 90% sequence identity, as these regions are unreliable when short reads are used to call variants.
+ -s tells the model to calculate synonymous variant density, which you may not want as the model takes a lot longer to do so. If you do this without the -c option for CpG, it will calculate a model based on coverage-scaled length and synonymous density instead of CpG content.
+ -u removes synonymous density from the model and just puts it as a field in the BED file.  We did this for our model version, CpG content only.
+ -e will use the version of regions.sh that incorporates ExAC v1 (assuming you downloaded it in getfiles.sh).
+ -g will use the version of regions.sh that incorporates gnomAD (or ExAC v2), downloaded by default in getfiles.sh.
+ -c triggers CpG density input into the model.  This is the version used in the manuscript.  It will always be calculated.
+ -n will prevent singletons from being used in the model.  If you do not believe in variants that occur only once in your VCF input file, this may be an option you are interested in.
+ -f will use a FASTA file of your choice to calculate CpG density.  The default is the GRCh37/hg19 FASTA downloaded from our s3 server by getfiles.sh. 
+ -w is the version we used in the manuscript where variants are separate from CCRs.  Variants are merged together if overlapping INDELs and treated as one large region of non-constraint.  It may be you want to incorporate variants into the regions, but the code to do so is quite old.  We recommend that you use this version.
+ -d allows you to pass a coverage depth and a coverage cutoff into the model.  As in the manuscript the default is 10x at 50%.  So, 10 and 0.5.  Separate inputs with spaces.
+ -r allows you to only run the model on one chromosome or a couple if you like.  We used this option to run CCRs just for the X chromosome quickly, as in the manuscript. Like -d, separate inputs with spaces.
+ -q will let you exclude certain chromosomes from your analysis.  The default is Y but as in the manuscript, we recommend using this to get autosomal CCRs and excluding X and Y.  But if you want everything just run with a blank.  Like -d, separate inputs with spaces.

### About the model code

`regions.sh` will run `exac-regions.py` and `resid-plot.py` (it also runs `weightpercentile.py`) using the above.  If you want more information on the inputs to those individual functions, just run the python scripts with -h to get the descriptions of all the inputs.  You do not have to use regions.sh to run the model.

`exac-regions.py` is the python code that generates the regions.  It utilizes the functions in the script `utils.py` to do so.  The functions in `utils.py` would need to be modified if you wanted to change the way regions are calculated.  `exac-regions.py` runs using parallel processing, so debugging may be difficult for those less experienced.

`resid-plot.py` is the code to create the model and calculate synonymous density.  Synonymous density calculation currently takes a very long time, so it may be in the interest of the user to forgo the calculation if iterating through several different model versions.

`weightpercentile.py` calculates the length-scaled percentiles off of the far-less balanced ones generated by `resid-plot.py`.  It is highly recommened you use this script to create the final percentiles. 
