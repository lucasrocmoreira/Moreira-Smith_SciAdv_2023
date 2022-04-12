# Code from: Parallel genomic signatures of local adaptation across a continental-scale environmental gradient

Code used for analyses in Moreira et al. 2022:

### Read alignment, variant calling, and filtering

* Code available [here](https://github.com/lucasrocmoreira/Moreira-et-al-2022).

### Genotype-environment association analysis

* `LFMM.R`: performs genotype-environment association analysis using [LFMM 2](https://bcm-uga.github.io/lfmm/index.html).

* `LFMM_intersection.R`: finds a set of (parallel) SNPs intersecting two sets of candidate SNPs and estimates probability of sharing.

* `LFMM-PCAdapt_intersection.R`: finds a set of SNPs intersecting three different methods.

### Signatures of selective sweep

* `H-scan_outlier.R`: detects the top outliers from [H-scan](https://messerlab.org/resources/) analysis.

### Population structure outliers

* `PCAdapt.R`: detects population struction outliers using [PCAdapt](https://bcm-uga.github.io/pcadapt/index.html).

### FST-outlier analysis

* `ANGSD_SAF.sh`: estimates allele frequencies directly from genotype likelihoods using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

* `ANGSD_realSFS.sh`: generates a folded 2d site frequency spectrum (SFS) using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

* `ANGSD_FST.sh`: estimates FST across sliding windows using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

* `Global_FST.sh`: estimates genome-wise FST using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

* `Count_FSToutliers.R`: counts the number of outlier FST windows across several pairwise population comparisons. 

* `Shared_FST_outliers.R`: identifies outlier windows shared between species. 

### Gene ontology enrichment

* `Manipulating_annotation.R`: transforms a gff file into a easier-to-manipulate table. 

* `SNP_annotation.R`: merges vcf data with gene annotation.

* `Extracting_GO_names_from_GO_ids.R`: goes from GO IDs to GO term name.

* `Extracting_annotation_info_from_outliers.R`: extracts annotation information from FST outlier windows.

* `Gene_ontology_enrichment.R`: performs gene ontology enrichment analysis.

### Vizualization

* `Fst_Manhattan_plot.R`: produces Manhattan plots from FST estimates.

* `Genotype_plot.R`: plots genotypes across a genomic segment of choice.
