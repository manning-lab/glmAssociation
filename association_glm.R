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

# opend gds file
gds.data <- seqOpen(gds.file)

# get the sample ids
sample.id <- seqGetData(gds.data, "sample.id")

# make sure we dont have any extra sample ids
phenotype.slim <- phenotype.slim[phenotype.slim[,id.col] %in% sample.id,]

# need to assign NA rows for sample ids that are in gds but not pheno
# get the ids
no.pheno <- sample.id[!(sample.id %in% phenotype.slim[,id.col])]

if (length(no.pheno) > 0){
  # make an na df to merge
  no.pheno.rows <- data.frame(sample.id=no.pheno, stringsAsFactors=F)

  # ensure that we merge on the right column
  names(no.pheno.rows) <- id.col

  # add the rows to pheno
  phenotype.slim[(nrow(phenotype.slim) + 1):(nrow(phenotype.slim) + nrow(no.pheno.rows)), names(no.pheno.rows)] <- no.pheno.rows

}

# order the rows to match gds
phenotype.slim <- phenotype.slim[match(sample.id,phenotype.slim[,id.col]),,drop=F]

# change col name, needed for seqvardata
all_vals[1] <- "sample.id"
names(phenotype.slim) <- all_vals

# metadata for df
meta <- data.frame(labelDescription=all_vals, row.names=names(phenotype.slim))

# make final sample df
pheno.anno <- AnnotatedDataFrame(phenotype.slim, meta)
reg.in <- SeqVarData(gds.data, pheno.anno)

# subset to some ids if given
if (!(sample.file == "NA")){
  sample.ids <- unique(readLines(sample.file))
  seqSetFilter(gds.data, sample.id=sample.ids, action="intersect", verbose=TRUE)
}

# Filter by desired MAC
gds.freq <- seqAlleleFreq(gds.data, .progress=TRUE)
gds.maf <- pmin(gds.freq, 1-gds.freq)
gds.mac.filt <- 2 * gds.maf * (1-gds.maf) * length(seqGetData(gds.data,"sample.id")) >= mac

# Filter to snps with mac greater than threshold
seqSetFilter(gds.data, variant.sel=gds.mac.filt, action="intersect", verbose=TRUE)

# Organize data for output
id <- seqGetData(gds.data,"variant.id")
chr <- seqGetData(gds.data,"chromosome")
pos <- seqGetData(gds.data,"position")
ref <- as.character(ref(gds.data))
alt <- as.character(unlist(alt(gds.data)))
maf <- gds.maf[gds.mac.filt]
snps.pos <- data.frame(id,chr,pos,ref,alt,maf)

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

# run regression
reg.out <- regression(reg.in,
                      outcome = outcome.name, 
                      covar=covariates, 
                      model.type=test)
# close gds
seqClose(gds.data)

# merge results with snp data
assoc <- merge(reg.out, snps.pos, by.x = "variant.id", by.y = "id")

# save results
save(assoc, file=paste(label, ".assoc.RData", sep=""))

