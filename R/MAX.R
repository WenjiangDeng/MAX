## Take the arguments
workdir=NULL
design.matrix="X_matrix.RData"
#core = 16 #default
hg="hg38"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="workdir") workdir=as.character(res[2])
	if (res[1]=="mut") mut=as.character(res[2])
	if (res[1]=="gtf") gtf=as.character(res[2])
	if (res[1]=="ref") ref=as.character(res[2])
	if (res[1]=="hg") hg=as.character(res[2])
}

options(stringsAsFactors=FALSE)
setwd(workdir)

if (!suppressMessages(suppressWarnings(require(GenomicFeatures)))){ 
if (!suppressMessages(suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))))
install.packages("BiocManager")
BiocManager::install("GenomicFeatures")}
require(GenomicFeatures)

if(hg=='hg38'){
if (!suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg38)))) {
if (!suppressMessages(suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")}
suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg38))) }#library(BSgenome.Hsapiens.UCSC.hg38)}

if(hg=='hg19'){
if (!suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19)))){
 if (!suppressMessages(suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19)))
}#library(BSgenome.Hsapiens.UCSC.hg19)}

if (!suppressMessages(suppressWarnings(require(foreach)))) install.packages('foreach')
if (!suppressMessages(suppressWarnings(require(doParallel)))) install.packages('doParallel')
suppressMessages(suppressWarnings(require(foreach)))
suppressMessages(suppressWarnings(require(doParallel)))
suppressMessages(suppressWarnings(require(Biostrings)))
##input files 1.mutation list: mutation-list.txt
library(data.table)
mut.list = fread(mut)
mut.list = as.data.frame(mut.list)
#mut.list = read.table(mut,header=TRUE)
index = paste(mut.list[,1],mut.list[,2],mut.list[,3],mut.list[,4],mut.list[,5],mut.list[,6])
mut.list = mut.list[!duplicated(index),]
mut.list$MutID=c(1:nrow(mut.list))
gene.list = unique(mut.list[,6])

#Nghia /05Feb2022:
write.table(mut.list,file="MAX_mut_list_keepMutID.txt",col.names=TRUE,quote=FALSE, row.names=FALSE,sep="\t")

## generate sqlite from gtf file

gtfFile=gtf
gtfSqliteFn="tmp_gtf.sqlite"

gtfTxdb <- makeTxDbFromGFF(file=gtfFile,
                 format="gtf",
                 dataSource=paste("Link to the source",sep=""),
                 organism="Homo sapiens")
saveDb(gtfTxdb,file=gtfSqliteFn)
cat("\nSqlite file generated\n")
###########
gtfSqlite="tmp_gtf.sqlite"

anntxdb <- loadDb(gtfSqlite)
genes.all = genes(anntxdb, single.strand.genes.only = FALSE )
tx.all = transcripts(anntxdb)
exon_all=exons(anntxdb)

genes.tx.all = suppressMessages(suppressWarnings(select(anntxdb, keys=names(genes.all), columns=c("GENEID", "TXNAME","EXONID","EXONSTART","EXONEND","EXONCHROM","EXONSTRAND"), keytype = "GENEID")))

## get the wild-type sequences for gene.list
if(hg=='hg38')fasta_genome<-BSgenome.Hsapiens.UCSC.hg38
if(hg=='hg19')fasta_genome<-BSgenome.Hsapiens.UCSC.hg19

g=gene.list
mygene=genes.tx.all[genes.tx.all$GENEID %in% g,] ##mygene strand info
mygene_txlist=unique(mygene$TXNAME)
save(mygene,file='gene_tx_tmp.RData')


mySeq=NULL
for (i in 1: length(mygene_txlist)){
  tx=mygene_txlist[i]
  mytx=genes.tx.all[genes.tx.all$TXNAME==tx,]
  mytx=mytx[order(mytx$EXONSTART),]
  mytx$len=mytx$EXONEND-mytx$EXONSTART+1

  chrID=mytx$EXONCHROM
  startR=mytx$EXONSTART
  endR=mytx$EXONEND
  gr = GRanges(seqnames=chrID, ranges=IRanges(start=startR, end=endR))
  exonseq=as.character(getSeq(fasta_genome, gr))
  txseq=paste0(exonseq,collapse = "")
  mySeq=c(mySeq,txseq)
}

mySeq = DNAStringSet(mySeq)
mySeq_name=paste0(mygene_txlist," ",g)
## ####################### ########## g here not correct
names(mySeq) = mySeq_name

writeXStringSet(mySeq, "wild_type_subset_genes.fa") ########## folder
### wild_type_subset_genes.fa will always be forward strand
cmd = as.character("awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' < wild_type_subset_genes.fa > tmp.fa;
tail -n +2 tmp.fa > wild_type_subset_genes.fa;
rm tmp.fa")
system(cmd)


#wild-type
mywt=readDNAStringSet("wild_type_subset_genes.fa")
## mygene has strand info
tempIdx=0
for(i in 1:length(gene.list))
{
	gene1 = gene.list[i]
	mut1 = mut.list[mut.list$Gene==gene1,]
	strand1 = unique(mygene[mygene$GENEID==gene1,'EXONSTRAND'])
	chrID = unique(mygene[mygene$GENEID==gene1,'EXONCHROM'])
	
	for(j in 1:nrow(mut1)){	
		pos_start = mut1[j,2]
		if((mut1[j,3]-mut1[j,2]+1)<20)
		{gr = GRanges(seqnames=chrID, ranges=IRanges(start=pos_start, end=pos_start+20))
		t1=as.character(getSeq(fasta_genome, gr))}
		if((mut1[j,3]-mut1[j,2]+1)>=20)
		{gr = GRanges(seqnames=chrID, ranges=IRanges(start=pos_start, end=mut1[j,3]+20))
		t1=as.character(getSeq(fasta_genome, gr))}
		
		t2 = mut1[j,5]
		gr = GRanges(seqnames=chrID, ranges=IRanges(start=mut1[j,3]+1, end=mut1[j,3]+20))
		t2.1=as.character(getSeq(fasta_genome, gr))
		t2 = paste0(t2,t2.1)# 
		
		tempIdx=tempIdx+1
		tempFile=paste0("tmp_wild_type_subset_genes_",tempIdx,".fa")		
		system(paste0("cp wild_type_subset_genes.fa ",tempFile))
		cmd1 = paste0("sed -i \"s/",t1,"/",t2,"/g\" ",tempFile,";")
		system(cmd1)
		#system("cat wild_type_subset_genes_tmp.fa >> Mutant.fa;rm wild_type_subset_genes_tmp.fa;")
		## Nghia /05Feb2022:
		mymut=readDNAStringSet(tempFile)
#		names(mymut)[1]="hello"
#		cat("\n",which(mymut!=mywt))
		mymut2=mymut
		mymut2_txname=sapply(names(mymut2),function(x) strsplit(x," ")[[1]][1])
		mymut2_genename=sapply(names(mymut2),function(x) strsplit(x," ")[[1]][2])
		mymut2_txname=paste0(mymut2_txname,"_mut_",mut1$MutID[j], " ",mymut2_genename)
		names(mymut2)=mymut2_txname
		writeXStringSet(mymut2,tempFile)
		system(paste0("cat ",tempFile," >> Mutant.fa;"))
	}  #sed -i 's/\..*$//' Mutant.fa.tmp 	
}
## rename mutant isoform
#cmd = as.character("awk '{if($1~/^>/)printf(\"%s_mut%s %s\\n\",$1,NR,$2);if($1!~/^>/)print $0;}' Mutant.fa > Mutant.fa.1;")
#system(cmd)
## Nghia /05Feb2022: no need the rename step
system("cat Mutant.fa > Mutant.fa.1")

#add wide type isoforms
#system("cat wild_type_subset_genes.fa Mutant.fa.1  > Mutant2.fa
#echo ' ' >>Mutant2.fa")
system("cat wild_type_subset_genes.fa Mutant.fa.1  > Mutant2.fa")

# put sequences into 1 line
cmd = as.character("awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' < Mutant2.fa > tmp.fa;
tail -n +2 tmp.fa > Mutant2.fa;
rm tmp.fa")
system(cmd)

# remove duplicated sequences
tt=readLines('Mutant2.fa');tt2 = unlist(lapply(tt,function(x)if(nchar(x)>3) return(x)))#remove duplicated sequences
index = duplicated(tt2);tt3=NULL
for(i in seq(1,length(tt2),2))if(!index[i+1])tt3=c(tt3,tt2[i],tt2[i+1])
writeLines(tt3,'Mutant3.fa')

