### association_glm.R
# Description: This function performs an association test to generate p-values for each variant included.
# Inputs:
# gds.file : a genotype file containing data for all samples are variants to be tested (.gds)
# phenotype.file : tsv or csv of phenotypes for samples contained in the gds.file (.csv, .tsv)
# outcome.name : outcome to be tested (string)
# covariate.string : comma separated list of covariates to include in the model (string)
# id.col : column name in phenotype file with sample ids (string)
# label : prefix for output filename (string)
# test : statistical test (linear, logistic, firth)
# sample.file : text file with one sample id per line to include in analysis (optional, .txt)
# mac : minimum minor allele count for variants to be included in analysis (int, default = 5)
# variant.range : comma separated list of variant ranges in format <chr#>:pos_start-pos_end
# Outputs:
# assoc : an RData file of associations results (.RData)

####### Testing inputs ########
# gds.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/gds_files/freeze.5b.chr10.pass_and_fail.gtonly.minDP10.chunk1.gds"
# phenotype.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/Pooled_MIXED_TM_19JAN18_T2D_freeze5b_harmonized_ancestry.csv"
# outcome.name <- "t2d_ctrl"
# covariate.string <- "last_exam_age,study,sex"
# id.col <- "topmedid"
# label <- "linear_testing"
# test <- "linear"
# sample.file <- "/Users/tmajaria/Documents/projects/topmed/data/test_inputs/phenotypes/first_500_samples.txt"
# mac <- 5
# variant.range <- "NA"
##############################

# Load packages
library(data.table)
library(Biobase)
library(SeqVarTools)
library(dplyr)

# Parse input arguments
input_args <- commandArgs(trailingOnly=T)
gds.file <- input_args[1] 
phenotype.file <- input_args[2]
outcome.name <- input_args[3]
covariate.string <- input_args[4]
id.col <- input_args[5]
label <- input_args[6]
test <- input_args[7]
sample.file <- input_args[8]
mac <- as.numeric(input_args[9])
variant.range <- input_args[10]

# get covariates
if (!(covariate.string == "NA")){
  covariates <- unlist(strsplit(covariate.string, ","))
} else {
  covariates <- c()
}

# check if test is valid of Score, wald, linear, logistic, firth
test.vals <- c("linear", "logistic", "firth")
test.ind <- match(tolower(test), tolower(test.vals))
if (!(is.na(test.ind))){
  test <- test.vals[test.ind]
} else {
  stop("Unrecognized test argument, must be one of: linear, logistic, firth")
}

# Print input arguments
print_ <- function(to_print) {
  print(paste(to_print, collapse=" "))
}

print_(c("gds.file", gds.file))
print_(c("phenotype.file", phenotype.file))
print_(c("outcome.name", outcome.name))
print_(c("covariate.string", covariate.string))
print_(c("id.col", id.col))
print_(c("label", label))
print_(c("test", test))
print_(c("sample.file", sample.file))
print_(c("mac", mac))
print_(c("variant.range", variant.range))

# load phenotype data
phenotype.data <- fread(phenotype.file,header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)

# make sure there are no duplicated ids
phenotype.data <- phenotype.data[!duplicated(phenotype.data[,id.col]),]

# get vector of all fields needed in phenotype file
all_vals <- c(id.col,outcome.name,covariates)

# subset by only these fields
phenotype.slim <- na.omit(as.data.frame(phenotype.data[,all_vals,drop=F]))

# subset by sample ids if they are given
if (!(sample.file == "NA")){
  sample.ids <- unique(readLines(sample.file))
  phenotype.slim <- phenotype.slim[phenotype.slim[,id.col] %in% sample.ids,na.omit(all_vals,drop=F)]
}

# get the right naming convention for the phenotype file
names(phenotype.slim)[1]<- "sample.id"

# put in annotated data frame format for regression function
phenotype.anno <- AnnotatedDataFrame(data = phenotype.slim)

# Make sure that the data frame has the right fields
gds.data <- seqOpen(gds.file)
sample.id <- seqGetData(gds.data, "sample.id")
sample.id <- data.frame(sample.id, stringsAsFactors = F)
pData(phenotype.anno) <- left_join(sample.id, pData(phenotype.anno), by="sample.id")

# Filter by desired MAC
seqSetFilter(gds.data,sample.id=phenotype.slim$sample.id, action="intersect", verbose=TRUE)
gds.freq <- seqAlleleFreq(gds.data, .progress=TRUE)
gds.maf <- pmin(gds.freq, 1-gds.freq)
gds.mac.filt <- 2 * gds.maf * (1-gds.maf) * length(phenotype.slim$sample.id) >= mac

# Filter to snps with mac greater than threshold
seqSetFilter(gds.data, variant.sel=gds.mac.filt, action="intersect", verbose=TRUE)

# Organize data for output
id <- seqGetData(gds.data,"variant.id")
chr <- seqGetData(gds.data,"chromosome")
pos <- seqGetData(gds.data,"position")
ref <- as.character(ref(gds.data))
alt <- as.character(unlist(alt(gds.data)))
snps.pos <- data.frame(id,chr,pos,ref,alt)

# If variant range is input, subset by var range
if (!(variant.range == "NA")){
  variant.range = as.list(unlist(strsplit(variant.range,",")))
  variant.range.list <- lapply(lapply(variant.range, function(x) unlist(strsplit(x,":"))), function(y) c(y[1],unlist(strsplit(y[2],"-"))))
  
  var.tokeep.id <- c()
  for (rng in variant.range.list){
    cur.var <- subset(snps.pos, chr == rng[1] & pos >= as.numeric(rng[2]) & pos <= as.numeric(rng[3]))
    var.tokeep.id <- c(var.tokeep.id, cur.var$id)
  }
  var.tokeep.id <- unique(var.tokeep.id)
  
  seqSetFilter(gds.data, variant.id=var.tokeep.id, action="intersect", verbose=TRUE)
  snps.pos <- snps.pos[snps.pos$id %in% var.tokeep.id,]
}

seqClose(gds.data)
reg.in <- SeqVarData(gds.file, sampleData = phenotype.anno)

# why is this taking so long?
reg.out <- regression(reg.in,
                      outcome = outcome.name, 
                      covar=covariates, 
                      model.type=test)

assoc <- cbind(snps.pos, reg.out)
save(assoc, file=paste(label, ".assoc.RData", sep=""))

