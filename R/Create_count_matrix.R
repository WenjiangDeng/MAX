## Take the workdir and core arguments
workdir=NULL
design.matrix="X_matrix.RData"
core = 8 #default

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="workdir") workdir=as.character(res[2])
	if (res[1]=="core") core=as.numeric(res[2])
  if (res[1]=="design.matrix") design.matrix=as.character(res[2])
}

cat("\n Create_count_matrix.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n workdir: ",workdir)
cat("\n core: ",core)
cat("\n design.matrix: ",design.matrix)
cat("\n ----------------------------------------------------- ")

source("/path/to/R/Rsource.R")

options(stringsAsFactors=FALSE)
setwd(workdir)
load(design.matrix)

flist = list.files(workdir,pattern="eqClass.txt",recursive=TRUE,full.names = TRUE)

for(j in 1:length(flist))
{
	doc1 = flist[j]
	a = read.table(doc1,header=TRUE)
	doc2 = gsub('eqClass.txt','',doc1)
	a_merge = NULL


for(i in unique(a$eqClass))
	{ 
	 b = NULL
	 b = subset(a,eqClass==i)	 
	 
	 tt = mutant_isoforms
	 
	for(isoform in paste0(tt,'_mut'))
	 {
	 index119 = grep(isoform,b$Transcript)
	 if(length(index119)>1)
		{
		 b119 = colSums(b[index119,2:6])
		 b119[2:5] = b119[2:5]/length(index119)
		 b119 = matrix(b119,nrow=1)
		 b119 = data.frame(Transcript=paste0(isoform,'_all'),b119)
		 colnames(b119) = names(b)
		 b = rbind(b119,b[grep(isoform,b$Transcript,invert=T),])
		}
	 }
	
		a_merge = rbind(a_merge,b)
		
	}
	if(is.null(a_merge)){next}	
	rownames(a_merge) = c(1:nrow(a_merge))
	x = a_merge
#	x$Transcript = gsub("_mut[0-9]*[0-9]",'_mut_all',x$Transcript)
  #####
  #Nghia 05Feb2022:
  x$Transcript = gsub("_mut_[[:print:]]+$",'_mut_all',x$Transcript)#letters, numbers, punctuation, and whitespace.
	#need to merge the same eqclasses after reduction
  x_reduced=NULL
  for(i in unique(x$eqClass)){
  x1 = NULL
  x1 = subset(x,eqClass==i)
  x1 = x1[!duplicated(x1$Transcript),]
  x1 = x1[order(x1$Transcript),]
  x_reduced=rbind(x_reduced,x1)
  }
  x=x_reduced
  x_reduced=NULL
	#####	

x_merge = NULL

for(i in unique(x$eqClass))
	{ 
	 x1 = NULL
	 x1 = subset(x,eqClass==i)
	 x1 = x1[order(x1$Transcript),]
	 x1 = cbind(x1,Isoform=paste(x1$Transcript,collapse='_'))
	 x_merge = rbind(x_merge,x1)
	}

z_merge = NULL

for(i in 1:length(unique(x_merge$Isoform)))
	{ 
	 b = NULL
	 b = subset(x_merge,Isoform==unique(x_merge$Isoform)[i])[,-7]
     b1 = rbind(tapply(b$Weight,b$Transcript,sum),tapply(b$Count,b$Transcript,sum),tapply(b$EffLength,b$Transcript,mean),tapply(b$RefLength,b$Transcript,mean))
	 b1 = t(b1)
	 b2 = data.frame(Transcript=rownames(b1),b1,eqClass=i)
	 colnames(b2) = colnames(b)
	 z_merge = rbind(z_merge,b2)
	}
	rownames(z_merge) = c(1:nrow(z_merge))

	write.table(z_merge,file=paste0(doc2,'Mutated_Combined_eqclass.txt'),quote=F,row.names=F,sep='\t')
### Merge eqClass

	
}

setwd(workdir)

if(!dir.exists("Ycount")) dir.create("Ycount")
setwd(paste(workdir,"/Ycount",sep=""))
flist = list.files(workdir,pattern="Mutated_Combined_eqclass.txt",recursive=TRUE,full.names = TRUE)
tx_length = flist[1]
#initialization for parallel computing
library(foreach)
library(doParallel)
registerDoParallel(cores=core)

res=foreach(id = 1:length(flist),.combine=c) %dopar% {
  y = NULL
  y = crpcount(flist[id])
  samplename = paste(y$samplename,".RData",sep="")

  save(y,file=samplename)
  return(flist[id])
}
res=NULL

flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
npat = sapply(CRP,nrow)   # number of occupancy patterns per clusterloc2 = which(npat>1) 
loc2 = which(npat>1) 
CCRP1 = CCRP[loc2]
Y=NULL
for(id in 1:length(flist)){
 cat("Merging results from sample ",flist[id],' ...\n')
 load(flist[id])
 if(id==1){
  Y = y[[1]][loc2]
 }
 if(id>1){
  y1 = y[[1]][loc2]
  for(i in 1:length(Y)){
   Y[[i]] = cbind(Y[[i]],sample1=y1[[i]][,'sample1'])
   }
  }
 }



samplename1 = NULL
for(id in 1:length(flist))
 {
  s.1 = strsplit(flist[id],"/")[[1]]
  s = s.1[length(s.1)]
  s = gsub(".RData","",s)
  samplename1 = c(samplename1,s)
 }

 setwd(workdir)

 for(i in 1:length(Y))
 {
 y2 = Y[[i]]
 xloc = which(colnames(y2) != "sample1")
 y2.1 = y2[,-xloc]
 y3 = cbind(CCRP1[[i]],y2.1)
 Y[[i]] = as.matrix(y3)
 }

save(Y,samplename1,tx_length,file='Ycount.RData')

cat("\n...Done...\n")
