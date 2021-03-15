# Example of generating simulated data of transcripts using polyester packages
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

fasta_file=as.character(args[1])
outdir = as.character(args[2])


if(!require(polyester))
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("polyester")
}
if(!require(Biostrings))
{
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("Biostrings")
}

library(polyester)
library(Biostrings)

fasta = readDNAStringSet(fasta_file)
# generate reads with coverage of 2
unif.countmat=20*width(fasta)
# generate at least 10000 reads
unif.countmat[unif.countmat<10000]=10000
unif.countmat=as.matrix(unif.countmat)
simulate_experiment_countmat(fasta_file, readmat=unif.countmat, outdir=outdir,error_rate=0, strand_specific=FALSE) 
rm(list=ls())

