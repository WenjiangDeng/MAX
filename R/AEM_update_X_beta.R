## 27 Jan. 2021 / Wenjiang:
# refine the codes to save final result
## 27 Nov 2019 / Nghia:
# add standard error for the estimates
## 04 June 2019 / Nghia:
# add parameter: isoform.method=average/total to report the expression of the individual members of a paralog i) average (default) or ii) total from the paralog set
## 01 Apr 2019 / Wenjiang:
# add "merge.paralogs" parameter to turn on/off the paralog merging in XAEM. The default is off, which will generate the same set of isoforms between different projects. To turn it on, just add "merge.paralogs=TRUE"
# Example of command: Rscript buildCRP.R in=eqClass.txt isoform.out=X_matrix.RData merge.paralogs=TRUE

## Take the workdir and core arguments
workdir=NULL
core = 8 #default
design.matrix="X_matrix.RData"
#isoform.method="average" #  "average" or "total"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")
remove.ycount=TRUE

for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="workdir") workdir=as.character(res[2])
	if (res[1]=="core") core=as.numeric(res[2])
  if (res[1]=="design.matrix") design.matrix=as.character(res[2])
#  if (res[1]=="isoform.out") fout=as.character(res[2])
  if (res[1]=="max.out") foutr=as.character(res[2])
#	if (res[1]=="merge.paralogs") merge.paralogs=as.logical(res[2])
#  if (res[1]=="isoform.method") isoform.method=as.character(res[2])
  if (res[1]=="remove.ycount") remove.ycount=as.logical(res[2])
}

cat("\n AEM_update_X_beta.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n workdir: ",workdir)
cat("\n core: ",core)
cat("\n design.matrix: ",design.matrix)
#cat("\n isoform.out: ",fout)
cat("\n max.out: ",foutr)
#cat("\n isoform.method: ",isoform.method)
cat("\n remove.ycount: ",remove.ycount)
cat("\n ----------------------------------------------------- ")


source("/path/to/R/Rsource.R")

options(stringsAsFactors=FALSE)
if(!dir.exists(workdir))dir.create(workdir)
setwd(workdir)

#load input data
load(design.matrix)
load("Ycount.RData")
#set parallel
library(foreach)
library(doParallel)
ncores = detectCores()
nc = min(ncores,core)     # use 8 or 16 as needed!!
cl <- makePSOCKcluster(nc)   #
registerDoParallel(cl)

##### start from here
X.y= Y


fun = function(crpdat, maxiter.X=1, modify=FALSE){
  #xloc = grep('N', colnames(crpdat))
  xloc = which(colnames(crpdat) != "sample1")
  X0 = matrix(crpdat[,xloc], ncol=length(xloc))
  Ymat = crpdat[,-xloc]
  est = AEM(X0, Ymat, maxiter.X=maxiter.X, modify=FALSE)
  return(est)
}
X.y= Y

cat("\n...Estimation using AEM algorithm...\n")
EST <- foreach(i= 1:length(X.y)) %dopar% fun(X.y[[i]], maxiter.X=1)
names(EST) = names(X.y)

x.all = list()
beta.all = list()
X.y=Y
for(i in 1:length(X.y))
{
 x1=b1=NULL
 x1 = EST[[i]]$X
 b1=EST[[i]]$BETA
 x.y =  X.y[[i]]
 xloc = which(colnames(x.y) != "sample1")
 colnames(x1) = colnames(x.y)[xloc]
 colnames(b1) = colnames(x.y)[xloc]
 rownames(x1) = rownames(x.y)
 rownames(b1) = samplename1
 x.all[[i]] = x1
 beta.all[[i]] = b1
}
names(x.all) = names(X.y)

updated_X = x.all
#save(updated_X,file='Updated_X.RData')

# beta.all
seq1 = Reduce(cbind,beta.all);seq1=t(seq1)
XAEM_count=NULL

##### process singletons
## singletons
estfun = function(mat){
  CRP.y = mat$crpcount
  #sum(!(names(CRP.y)==names(CRP))) 
  ## cluster info
  npat = sapply(CRP,nrow)   # number of occupancy patterns per cluster
  table(npat)
  loc1 = which(npat==1)  # clusters with 1 pattern

  TC1 = sapply(CRP.y[loc1],function(x) x[1, ncol(x)])
   names(TC1)=names(CRP.y[loc1])
  ## single tx
  est.all =  c(TC1)
  return(est.all)
}

flist = list.files(paste(workdir,"/Ycount/",sep=""),pattern="RData",recursive=TRUE,full.names = TRUE)
result_est=NULL
for(id in 1:length(flist)){ # call crpcount()
  load(flist[id])
  est1 = estfun(y)# estimation step
  result_est = cbind(result_est,est1)
  #cat("sample ",id,'\n')
}
if(nrow(result_est)>0)
{
colnames(result_est)=samplename1
XAEM_count = rbind(result_est,seq1)
}else{XAEM_count=seq1}

##### done with the estimation

### collect the list of txnames from count data, more than 1 transcripts if the isoform is a paralog
txList = sapply(rownames(XAEM_count),function(x){return(unlist(strsplit(x," ")))})
txNum=sapply(txList,length)
### get median txlength for paralog
txLen=sapply(txList, function(x){
  pick=names(txlength) %in% x
  return(median(txlength[pick]))
})
### compute TPM
#normalize count to length
isoform_lenNorm=apply(XAEM_count,2,function(x)return(x/txLen))
libsize_lenNorm=apply(isoform_lenNorm,2,sum)
XAEM_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
XAEM_tpm=t(XAEM_tpm)
#keep information of raw output of XAEM
MAX_count=XAEM_count
MAX_tpm=XAEM_tpm
attr(MAX_count,'conv')=NULL
isoform_count=MAX_count
isoform_tmp=MAX_tpm
save(isoform_count, isoform_tmp,file=foutr)

### expand isoforms can not separated from CRP 
### isoform.method (value=average/total)
#
#paralogID=which(txNum >1)
#paralog_count=XAEM_count[paralogID,]
#paralog_tpm=XAEM_tpm[paralogID,]
#
##get mapping ID
#if(length(paralogID)>0)
#{
#matchID=sapply(c(1:length(paralogID)), function(x) rep(paralogID[x],txNum[paralogID][x]))
#matchID=unlist(matchID)
##get isoform names
#matchNames=sapply(c(1:length(paralogID)), function(x) unlist(strsplit(names(paralogID[x])," ")))
#matchNames=unlist(matchNames)
##update to data XAEM_count
#expandDat=matrix(0,ncol = ncol(XAEM_count), nrow = sum(txNum[paralogID]))
#expandDat=XAEM_count[matchID,] #if (isoform.method=="total"):the counts of isoform members are equal to the count of paralog
#if (isoform.method=="average"){ #the counts of isoform members are equal to the average count of paralog
#  expandDat_txnum=txNum[match(names(matchID),names(txNum))]
#  expandDat2=apply(cbind(expandDat_txnum,expandDat),1,function(x) x[-1]/x[1])
#  expandDat=t(expandDat2)
#}
#rownames(expandDat)=matchNames
#isoform_count=XAEM_count
#isoform_count=rbind(isoform_count,expandDat)
#isoform_count=isoform_count[-paralogID,]
###recalculate TPM
#txLen2=txlength[match(rownames(isoform_count),names(txlength))]
#isoform_lenNorm=apply(isoform_count,2,function(x)return(x/txLen2))
#libsize_lenNorm=apply(isoform_lenNorm,2,sum)
#isoform_tpm=apply(isoform_lenNorm,1,function(x) return(x*1e6/libsize_lenNorm))
#isoform_tpm=t(isoform_tpm)
#
##export to file
#pick=names(txlength) %in% rownames(isoform_count)
#pick=which(!pick)
#if(length(pick)>0){
#  newDat=matrix(0,ncol = ncol(isoform_count), nrow = length(pick))
#  rownames(newDat)=names(txlength)[pick]
#  colnames(newDat)=colnames(isoform_count)
#  isoform_count=rbind(isoform_count,newDat)
#  isoform_tpm=rbind(isoform_tpm,newDat)
#}
#isoform_count=isoform_count[names(txlength),]
#isoform_tpm=isoform_tpm[names(txlength),]
#
##save(isoform_count,isoform_tpm,file=fout)
#}

### clean Ycount
if (remove.ycount){
  system("rm -rf Ycount")
  system("rm -rf Ycount.RData")
}

cat("\n...Done...\n")
