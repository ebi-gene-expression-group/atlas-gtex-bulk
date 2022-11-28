library(Rsamtools)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

bamPath <- args[1]

testPairedEndBam(bamPath)
