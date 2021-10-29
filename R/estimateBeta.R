### estimation
# if mutation-sample mapping is available, example command:
# - Rscript estimateBeta.R workdir=Ycount sampleMut=mutDat.txt out=MAX_isoform_expression.RData
# if mutation-sample mapping is NOT available, example command:
# - Rscript estimateBeta.R workdir=Ycount out=MAX_isoform_expression.RData

# default settings
sampleMutFn=NULL
workdir="Ycount"
fout="MAX2_isoform_expression.RData"

maxiter.X=1
maxiter.b=100
modify=FALSE

args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="workdir") workdir=res[2]
  if (res[1]=="sampleMut") sampleMutFn=res[2]
  if (res[1]=="out") fout=res[2] 
  if (res[1]=="maxiterX") maxiter.X=as.integer(res[2])
  if (res[1]=="maxiterb") maxiter.b=as.integer(res[2])
  if (res[1]=="modify") modify=as.logical(res[2])
}

source("Rsource.R")

flist = list.files(workdir,pattern="_Ycount.RData",recursive=TRUE,full.names = TRUE)

if (!is.null(sampleMutFn)){
  cat("\n Estimation of isoform expression (AEM)")
### using AEM
  flist2 = list.files(workdir,pattern="_Ycount.RData",recursive=TRUE,full.names = FALSE)
  sNameID=gsub("_Ycount.RData","",flist2)
  fileMutDat=read.csv(sampleMutFn,header=FALSE, sep="\t",stringsAsFactors = FALSE)
  sGroup=tapply(fileMutDat[,1],fileMutDat[,2],c)
  maxList=list()
  txAll=NULL
  for (g in 1:length(sGroup)){
    pick=which(sNameID %in% sGroup[[g]])
    flist3=flist[pick]
    gSnames=sNameID[pick]
    fn=flist3[1]
    load(fn)
    crpdat=Y[[1]]
    xloc = which(colnames(crpdat) != "sample1")
    X0 = crpdat[,xloc]
    txAll=unique(c(txAll,colnames(X0)))
    Ymat=NULL
    for (fn in flist3){
      load(fn)
      crpdat=Y[[1]]
      xloc = which(colnames(crpdat) != "sample1")
      X0 = crpdat[,xloc]
      Ymat = cbind(Ymat,crpdat[,-xloc,drop=FALSE])
    }
    colnames(Ymat)=gSnames

    estRes = AEM(X0, Ymat, maxiter.X=maxiter.X, maxiter.b=maxiter.b, modify=modify)
  #  convergeVec=c(convergeVec,estRes$beta.conv)
    th2=estRes$BETA
    colnames(th2)=colnames(X0)
    rownames(th2)=gSnames
    maxList[[g]]=th2
  }
  #create the final expression matrix
  maxMat=matrix(0,length(txAll),length(sNameID))
  dim(maxMat)
  rownames(maxMat)=txAll
  colnames(maxMat)=sNameID
  for (g in 1:length(maxList)){
    thMat=maxList[[g]]
    for (i in 1:nrow(thMat)){
      x=thMat[i,drop=FALSE,]
      maxMat[colnames(x),rownames(x)]=x
    }
  }
}else{
  ### Using EM
  cat("\n Estimation of isoform expression (EM)")
  maxList=list()
  txAll=NULL
  convergeVec=NULL
  for (fn in flist){
    load(fn)
    crpdat=Y[[1]]
    xloc = which(colnames(crpdat) != "sample1")
    X0 = crpdat[,xloc]
    Ymat = crpdat[,-xloc,drop=FALSE]
    estRes = AEM(X0, Ymat, maxiter.X=1, maxiter.b=maxiter.b, modify=modify)
    convergeVec=c(convergeVec,estRes$beta.conv)
    th2=estRes$BETA[1,]
    names(th2)=colnames(X0)
    txAll=unique(c(txAll,colnames(X0)))
    maxList[[samplename1]]=th2
  }
  maxMat=matrix(0,length(txAll),length(maxList))
  dim(maxMat)
  colnames(maxMat)=names(maxList)
  rownames(maxMat)=txAll
  for (i in 1:ncol(maxMat)){
    x=maxList[[colnames(maxMat)[i]]]
    maxMat[names(x),i]=x
  }
}

isoformCount=maxMat
save(isoformCount, file=fout)
