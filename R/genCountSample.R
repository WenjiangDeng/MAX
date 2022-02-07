## Rscript genCountSample.R xmatEq=$xmatEqFn sampleMut=$sampleMutFn sampleID=$sampleFn sampleEq=$inputFn YcountDir=$outdir

# default settings
xmatEqFn="Xmatrix/eqClass.txt"
sampleMutFn="mutDat_TP53.txt"
inputFn="sample_01_1.fasta/Mutated_eqClass.txt"
sampleFn="sample_01_1.fasta"
outdir="Ycount_MAX2"

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="xmatEq") xmatEqFn=res[2]
  if (res[1]=="sampleMut") sampleMutFn=res[2]
  if (res[1]=="sampleID") sampleFn=res[2]
  if (res[1]=="sampleEq") inputFn=res[2]
  if (res[1]=="YcountDir") outdir=res[2]
}

if(!dir.exists(outdir)) dir.create(outdir)

source("/path/to/R/Rsource.R")

#Nghia /05Feb2022:
mut.list=read.csv(sampleMutFn,header=TRUE, sep="\t",stringsAsFactors = FALSE)
s=mut.list$Samples
names(s)=mut.list$MutID
sMap=sapply(s,function(x) sort(unique(trimws(strsplit(x,",")[[1]]))))
s=unique(unlist(sMap))
myMut=NULL
for (i in 1:length(sMap)) if (sampleFn %in% sMap[[i]]){
  myMut=c(myMut,names(sMap)[i])
}


eqDat=read.table(xmatEqFn, sep="\t",header=TRUE, stringsAsFactors = FALSE)
txAll=unique(eqDat$Transcript)
pick=grep("mut",txAll)
txMut=txAll[pick]
txWt=txAll[-pick]
mutList=sapply(txMut,function(x) strsplit(x,"_mut_")[[1]][2])
mutList=unique(mutList)
names(mutList)=NULL
#myMut=mutList[c(5,6)]
selectedTxMut=sapply(myMut,function(x) txMut[grep(x,txMut)])
selectedTxMut=unique(unlist(c(selectedTxMut)))
txSelected=c(txWt,selectedTxMut)
eqList=lapply(unique(eqDat$eqClass), function(e){
  eqDat[eqDat$eqClass==e,]
})  
eqList2=lapply(eqList,function(eq){
  pick=eq$Transcript %in% txSelected
  eq1=eq[pick,drop=FALSE,]
  nCount=sum(eq1$Weight)
  if (nrow(eq1)>0) eq1$Count=nCount
  return(eq1)
})
s=sapply(eqList2,nrow)
eqList3=eqList2[s>0]
s=sapply(eqList3,nrow)
eqList3=lapply(eqList3, function(eq){
  eq=eq[order(eq$Transcript),drop=FALSE,]
  return(eq)
})
eqList3.names=sapply(eqList3, function(eq) paste0(eq$Transcript, collapse = " "))

eqList4=tapply(eqList3,eqList3.names,function(eqs){
  eq=eqs[[1]]
  res=lapply(eqs,function(eq) eq$Weight)
  res=do.call(rbind,res)
  eq$Weight=colSums(res)
  eq$Count=sum(eq$Weight)
  return(eq)
})
rawmat_eqmut=do.call(rbind,eqList4)
rownames(rawmat_eqmut)=NULL
eqIdx_old=unique(rawmat_eqmut$eqClass)
eqIdx=seq(length(eqIdx_old))
rawmat_eqmut$eqClass=eqIdx[match(rawmat_eqmut$eqClass,eqIdx_old)]
## now need to combine them into mut_all
rawmat=mergeMut(rawmat_eqmut,txWt)
rawmat$Transcript=as.character(rawmat$Transcript)
## get CRP and CCRP
res=getCRP(rawmat,H_thres=0.025,ncore=1)
CRP=res$CRP
CCRP=res$CRP

y = crpcount(inputFn)
npat = sapply(CRP,nrow)   # number of occupancy patterns per clusterloc2 = which(npat>1) 
loc2 = which(npat>1) 
Y = y[[1]][loc2]
samplename1=y$samplename
save(Y,samplename1,file=paste0(outdir,"/",sampleFn,'_Ycount.RData'))

cat("\n...Done...\n")