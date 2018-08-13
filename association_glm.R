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

# Load packages
library(data.table)
library(Biobase)
library(SeqVarTools)
library(dplyr)
library(tidyr)

.variantDF <- function(gds) {
  data.frame(variant.id=seqGetData(gds, "variant.id"),
             chromosome=seqGetData(gds, "chromosome"),
             position=seqGetData(gds, "position"),
             ref=refChar(gds),
             alt=altChar(gds),
             nAlleles=seqNumAllele(gds),
             stringsAsFactors=FALSE)
}
.expandAlleles <- function(gds) {
  .variantDF(gds) %>%
    separate_rows_("alt", sep=",") %>%
    rename_(allele="alt") %>%
    group_by_("variant.id") %>%
    mutate_(allele.index=~1:n()) %>%
    as.data.frame()
}


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

if (sum(gds.mac.filt) == 0){
  print("No SNPs pass MAC filter. Finished Association Step")
  assoc <- NA
  
} else {
  # Filter to snps with mac greater than threshold
  seqSetFilter(gds.data, variant.sel=gds.mac.filt, action="intersect", verbose=TRUE)

  # Filter to only passing variants
  var.ids <- seqGetData(gds.data,"variant.id")
  filt <- seqGetData(gds.data, "annotation/filter")
  var.ids <- var.ids[filt == "PASS"]

  if (length(var.ids) == 0){
    print("No SNPs pass MAC filter. Finished Association Step")
    assoc <- NA

  } else {
    seqSetFilter(gds.data, variant.id = var.ids, action="intersect", verbose=TRUE)  

    # Organize data for output
    snps.pos <- .expandAlleles(gds.data)[,c(1,2,3,4,5)]
    names(snps.pos) <- c("id","chr","pos","ref","alt")
    
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
    
    # merge results with snp data
    assoc <- merge(reg.out, snps.pos, by.x = "variant.id", by.y = "id")

    # get case/control
    pheno <- data.frame(sample = as.character(reg.in@sampleData@data$sample.id), outcome = reg.in@sampleData@data[,outcome.name], stringsAsFactors = F)
    
    # fix assoc format
    #MarkerName	chr	pos	ref	alt	minor.allele	maf	pvalue	n	Score.Stat	homref	het	homalt
    assoc$MarkerName <- paste(assoc$chr, assoc$pos, assoc$ref, assoc$alt, sep = "-")
    assoc$minor.allele <- "alt"
    
    if (test == "linear"){
      assoc$maf <- pmin(assoc$freq, 1-assoc$freq)
      assoc <- assoc[,c("MarkerName","chr","pos","ref","alt","minor.allele","maf",names(assoc)[endsWith(tolower(names(assoc)),"pval")],"n","Est","SE","Wald.Stat")]
      names(assoc) <- c("MarkerName","chr","pos","ref","alt","minor.allele","maf","pvalue","n","beta","stderr","test.stat")
    } else {
      if (test == "logistic"){
        assoc$maf <- (assoc$n0*assoc$freq0 + assoc$n1*assoc$freq1)/(assoc$n0+assoc$n1)
        assoc$maf <- pmin(assoc$maf, 1-assoc$maf)
        assoc$n <- assoc$n0 + assoc$n1
        assoc <- assoc[,c("variant.id","MarkerName","chr","pos","ref","alt","minor.allele","maf",names(assoc)[endsWith(tolower(names(assoc)),"pval")],"n","Est","SE","Wald.Stat")]
      } else {
        assoc$maf <- (assoc$n0*assoc$freq0 + assoc$n1*assoc$freq1)/(assoc$n0+assoc$n1)
        assoc$maf <- pmin(assoc$maf, 1-assoc$maf)
        assoc$n <- assoc$n0 + assoc$n1
        assoc <- assoc[,c("variant.id","MarkerName","chr","pos","ref","alt","minor.allele","maf",names(assoc)[endsWith(tolower(names(assoc)),"pval")],"n","Est","SE","PPL.Stat")]
      }
      
      assoc$Est <- exp(assoc$Est)
      names(assoc) <- c("variant.id","MarkerName","chr","pos","ref","alt","minor.allele","maf","pvalue","n","or","stderr","test.stat")
    
      # get the variants that pass both maf and pval threshold
      assoc.top_var <- assoc[(assoc$maf < 0.05 & assoc$pvalue < 0.01), "variant.id"]
      
      # set filter
      seqSetFilter(gds.data, variant.id = assoc.top_var, sample.id = pheno$sample)
      
      # get genotypes
      geno <- altDosage(gds.data)
      geno.ctrl <- geno[row.names(geno) %in% pheno[pheno$outcome == 0, "sample"],]
      geno.case <- geno[row.names(geno) %in% pheno[pheno$outcome == 1, "sample"],]
      rm(geno)
      
      # get counts per geno
      geno.ctrl.counts <- apply(geno.ctrl, 2, function(x) sum(x == 0, na.rm = T))
      geno.ctrl.counts <- data.frame(variant.id = names(geno.ctrl.counts), homref = as.numeric(as.character(geno.ctrl.counts)), stringsAsFactors = F)
      geno.ctrl.counts$het <- as.numeric(as.character(apply(geno.ctrl, 2, function(x) sum(x == 1, na.rm = T))))
      geno.ctrl.counts$homalt <- as.numeric(as.character(apply(geno.ctrl, 2, function(x) sum(x == 2, na.rm = T))))
      
      geno.case.counts <- apply(geno.case, 2, function(x) sum(x == 0, na.rm = T))
      geno.case.counts <- data.frame(variant.id = names(geno.case.counts), homref = as.numeric(as.character(geno.case.counts)), stringsAsFactors = F)
      geno.case.counts$het <- as.numeric(as.character(apply(geno.case, 2, function(x) sum(x == 1, na.rm = T))))
      geno.case.counts$homalt <- as.numeric(as.character(apply(geno.case, 2, function(x) sum(x == 2, na.rm = T))))
      
      # get to right format
      geno.counts <- data.frame(
        variant.id = geno.ctrl.counts$variant.id, 
        homref = paste0(geno.case.counts$homref, "/", geno.ctrl.counts$homref),
        het = paste0(geno.case.counts$het, "/", geno.ctrl.counts$het),
        homalt = paste0(geno.case.counts$homalt, "/", geno.ctrl.counts$homalt),
        stringsAsFactors = F
      )
  
      geno.counts <- data.frame(
        variant.id = geno.ctrl.counts$variant.id, 
        homref.case = geno.case.counts$homref,
        homref.control = geno.ctrl.counts$homref,
        het.case = geno.case.counts$het,
        het.control = geno.ctrl.counts$het,
        homalt.case = geno.case.counts$homalt,
        homalt.control = geno.ctrl.counts$homalt,
        stringsAsFactors = F
      )
      
      assoc <- merge(assoc, geno.counts, by.x = "variant.id", by.y = "variant.id", all.x = T)
      assoc[is.na(assoc)] <- ""
      assoc <- assoc[,c("MarkerName","chr","pos","ref","alt","minor.allele","maf","pvalue","n","or","stderr","test.stat","homref.case","homref.control","het.case","het.control","homalt.case","homalt.control")]
    } else {
      
    }
    # close gds
    seqClose(gds.data)


    # save results
    save(assoc, file=paste(label, ".assoc.RData", sep=""))

    }
  }