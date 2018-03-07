# Single Variant Association -- Linear/logistic/firth

## Description 

This workflow implements a single variant association analysis of genotype data with a single trait using linear/logistic/firth. The primary code is written in R using the SeqVarTools package for association testing.

### Authors

This workflow is produced and maintained by the [Manning Lab](https://manning-lab.github.io/). Contributing authors include:

* Tim Majarian (tmajaria@broadinstitute.org)

## Dependencies

### Workflow execution

* [WDL](https://software.broadinstitute.org/wdl/documentation/quickstart)
* [Cromwell](http://cromwell.readthedocs.io/en/develop/)

### R packages

* [Biobase](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html)
* [SeqVarTools](https://www.bioconductor.org/packages/release/bioc/html/SeqVarTools.html)
* [dplyr](http://dplyr.tidyverse.org/)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [qqman](https://cran.r-project.org/web/packages/qqman/index.html)

## Workflow elements

### association_glm.R

This function performs an association test to generate p-values for each variant included.

Inputs:
* gds.file : a genotype file containing data for all samples are variants to be tested (.gds)
* phenotype.file : tsv or csv of phenotypes for samples contained in the gds.file (.csv, .tsv)
* outcome.name : outcome to be tested (string)
* covariate.string : comma separated list of covariates to include in the model (string)
* id.col : column name in phenotype file with sample ids (string)
* label : prefix for output filename (string)
* test : statistical test (linear, logistic, firth)
* sample.file : text file with one sample id per line to include in analysis (optional, .txt)
* mac : minimum minor allele count for variants to be included in analysis (int, default = 5)
* variant.range : comma separated list of variant ranges in format chr#:pos_start-pos_end

Outputs:
assoc : an RData file of associations results (.RData)

### summary.R

Generate a summary of association results including quantile-quantile and manhattan plots for variants subseted by minor allele frequency (all variants, maf < 5%, maf >= 5%). Also generates CSV files of all variants and variants with P < pval_threshold.

Inputs:
* test : statistical test (linear, logistic, firth)
* pval_threshold : p-value threshold for the returning top associations, top association output will include only variants with a p-value less than the threshold (float, default = 0.0001)
* label : prefix for output filename (string)
* assoc : output of assocTest (Array[.RData])

## Other workflow inputs

* this_memory : amount of memory in GB for each execution of a task (int)
* this_disk : amount of disk space in GB to allot for each execution of a task (int)



