## 15 Otc 2019 / Nghia:
# - fix a bug when separe a CRP into more than 1 CRP
## 30 Apr 2019 / Nghia:
# - compute CRP from CRPCOUNT: fix the bug when colSums(mycrpCount)==0 
## 03 Apr 2019 / Nghia:
# - add H_thres for filtering connections with low counts between two transcripts
# Example of command: Rscript in=eqClass.txt out=X_matrix.RData H=0.025
## 02 Nov 2018 / Nghia:
# - grow fully TC, 
# - improve speed

### Build transcript cluster and CRP from simulated data
# Before running this script: 
# 1) generate simulated data of the whole transcriptome using genPolyesterSimulation.R
# 2) run GenTC for the simulated data
# Input: eqClass.txt from results of GenTC


### default settings
#eqClassFn="eqClass.txt"
eqClassFn=NA
fout='X_matrix.RData'
#H_thres=0.0
H_thres=0.025
workdir=NULL
args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args,"\n")

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1]=="in") eqClassFn=res[2]
  if (res[1]=="H") H_thres=as.double(res[2])
  if (res[1]=="out") fout=res[2]
  if (res[1]=="workdir") workdir=res[2]
}

setwd(workdir)

cat("\n buildCRP.R will run with the following parameter setting: ")
cat("\n ----------------------------------------------------- ")
cat("\n in: ",eqClassFn)
cat("\n out: ",fout)
cat("\n workdir: ",workdir)
cat("\n ----------------------------------------------------- ")

source("/path/to/R/Rsource.R")

########## Merge eqClass
setwd(workdir)
load('gene_tx_tmp.RData') 
doc1 = 'eqClass.txt'
a = read.table(doc1,header=TRUE)

a_merge = NULL

for(i in unique(a$eqClass))
	{ 
	 b = NULL
	 b = subset(a,eqClass==i)	 
	 
	 tt = unique(mygene$TXNAME)
	 
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
	#x$Transcript = gsub("_mut[0-9]*[0-9]",'_mut_all',x$Transcript)
  #####
  #Nghia 05Feb2022:
  x$Transcript = gsub("_mut_[[:print:]]+$",'_mut_all',x$Transcript)#letters, numbers, punctuation, and whitespace.
#	#need to merge the same eqclasses after reduction
#  x_reduced=NULL
#  for(i in unique(x$eqClass)){
#  x1 = NULL
#  x1 = subset(x,eqClass==i)
#  x1 = x1[!duplicated(x1$Transcript),]
#  x1 = x1[order(x1$Transcript),]
#  x_reduced=rbind(x_reduced,x1)
#  }
#  x=x_reduced
#  x_reduced=NULL
#	#####

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

	write.table(z_merge,file='Mutated_Combined_eqclass.txt',quote=F,row.names=F,sep='\t')
### Merge eqClass


library(foreach)
library(doParallel)
ncores = detectCores()
#nc = min(ncores,core)     # use 8 or 16 as needed!!
nc=min(ncores,8)
cl <- makePSOCKcluster(nc)   #
registerDoParallel(cl)


#get input from eqClass.txt in the output directory of GenTC
rawmat  = read.table(eqClassFn, header=TRUE, as.is=TRUE, sep='\t')

tx2eqc = tapply(rawmat$eqClass,rawmat$Transcript,c)  ## map: tx --> eqClass 
a = sapply(tx2eqc, length)
#table(a)
eqc2tx = tapply(rawmat$Transcript, rawmat$eqClass,c) ## map: eqClass --> tx

# neighbors of each tx
fn = function(i) {eqc=as.character(tx2eqc[[i]]);
unique(unlist(eqc2tx[eqc]))
}


system.time(NB <- foreach(i=1:length(tx2eqc)) %dopar% fn(i) )  ## about 25sec 
names(NB) = names(tx2eqc)

### build transcript clusters - new codes

#convert name of tx to index
txi=seq(length(NB))
NBi=NB
names(NBi)=txi
#map from tx to crp
t2c_map=unlist(NBi,use.names=FALSE)
t2c_mapi=match(t2c_map,names(NB))
t2c_mapi_group=rep(txi, lengths(NBi))
NBi2=tapply(t2c_mapi,t2c_mapi_group,c)
NBi=NBi2

f1=sapply(NBi,function(x) min(x))
f2=sapply(NBi,function(x) min(f1[x]))

growTimes=1
isOK=FALSE
repeat{  
  f2=sapply(NBi,function(x) min(f1[x]))  
  difNum=sum(f2!=f1)  
  growTimes=growTimes+1
  f1=f2
  if (difNum==0) isOK=TRUE
  if (isOK) break();
}

f1=names(NB)[f1]
names(f1)=names(NB)

NB2=tapply(names(f1),f1,c)
NB2=lapply(NB2,function(x) sort(x))
#create TC3
TC3=NB
TC3[names(f1)]=NB2[f1]


OTC <- sapply(TC3, paste, collapse=' ') # pasted version
names(OTC) = names(NB)
otcmap = as.list(OTC)
names(otcmap) = names(NB)

#output: list of clusters and tx->cluster map 
clust = names(table(OTC))  # clusters
OTC = strsplit(clust,split=' ')
names(OTC) = clust

save.image(file="tmp1.RData")

#get CRP count
system.time(
CRPCOUNT <- foreach(i=1:length(OTC)) %dopar%{
  myOTC=unlist(OTC[i])
  names(myOTC)=NULL
  #get binary codes of eqc
  myeqc=unique(unlist(sapply(myOTC,function(x) tx2eqc[x])))
  #eqc2tx[myeqc]
  bcode=lapply(eqc2tx[myeqc],function(x) as.integer(!is.na(match(myOTC,x))))
  #bcode

  #get corresponding count - new codes
  pick=which(rawmat$eqClass %in% myeqc)
  countname=paste(rawmat$Transcript[pick],rawmat$eqClass[pick],sep="__")
  count=rawmat$Weight[pick] 
  myname=expand.grid(myOTC,myeqc)
  myname=paste(myname[,1],myname[,2],sep="__")
  myCount=rep(0,length(myname))
  matchID=match(countname,myname)
  myCount[matchID]=count

  #create TC count matrix
  mycrpCount=matrix(unlist(myCount),nrow=length(bcode),ncol=length(myOTC),byrow=TRUE)
  colnames(mycrpCount)=myOTC
  rownames(mycrpCount)=sapply(bcode, function(x) paste(x,collapse=""))
 
  if (sum(mycrpCount)>0 & nrow(mycrpCount)>1){ # if there are more than 1 row
     ## sum(mycrpCount)>0 because some tx may be too short,<100 bp
	tmp.colnames = colnames(mycrpCount)
	tmp.rownames = rownames(mycrpCount)
    mycrpCount=as.matrix(mycrpCount[order(rownames(mycrpCount)),])
	colnames(mycrpCount) = tmp.colnames
	rownames(mycrpCount) = tmp.rownames[order(tmp.rownames)]
    #check if duplicated rownames
    repID=table(rownames(mycrpCount))
    repID=repID[which(repID>1)]
    if (length(repID)>0){
      rmID=which(rownames(mycrpCount) %in% names(repID))
      sumcrpCount=NULL
      for (j in 1:length(repID)){
	  
	    tmp.crp = as.matrix(mycrpCount[rownames(mycrpCount) %in% names(repID)[j],])
		
        sumcrpCount=rbind(sumcrpCount,colSums(tmp.crp))
      }
      rownames(sumcrpCount)=names(repID)
      mycrpCount=mycrpCount[-rmID,]
      mycrpCount=rbind(mycrpCount,sumcrpCount)
	  
	  if(nrow(mycrpCount)>1)
		{
		mycrpCount=mycrpCount[order(rownames(mycrpCount)),]
		mycrpCount=as.matrix(mycrpCount)
		colnames(mycrpCount) = tmp.colnames
		}
    }
  }
  
  return(mycrpCount)  
#  #extract CRP matrix
#  mycrp=t(t(mycrpCount)/colSums(mycrpCount))
#  mycrp
}
)

save.image(file="tmp2.RData")

#### now compute CRP, use H_thres (default=0) to filter out too low proportion sharing between two transcripts

CRP=list()
for (k in 1:length(CRPCOUNT)){
  mycrpCount=CRPCOUNT[[k]]

  #fix the bug when colSums(mycrpCount)==0 
  txSum=colSums(mycrpCount)
  mycrpCount=mycrpCount[,which(txSum>0),drop=FALSE]
  if(ncol(mycrpCount)==0) next(); 
  
  #x=mycrpCount
  y=t(t(mycrpCount)/colSums(mycrpCount))
  z=apply(y,1,max)
  pick=z>H_thres
  x1=mycrpCount[pick,drop=FALSE,]
  #decode eq in x1
  myclust=c(1:ncol(x1))
  for (i in 1:nrow(x1)){
    x2=which(x1[i,] >0)
    x3=which(myclust %in% myclust[x2])
    myclust[x3]=min(myclust[x3])
  }
  #generate new CRP
  clustID=as.integer(names(table(myclust)))
  #newx=list()
  for (i in 1:length(clustID)){  
    pick=which(myclust == clustID[i])
    X=mycrpCount[,pick,drop=FALSE]
    X=X[rowSums(X)>0,drop=FALSE,]
    rownames(X)=NULL
    bcode=apply(X,1,function(x) paste(as.integer(x>0),collapse=""))
    ubcode=unique(bcode)
    cbcode=match(bcode,ubcode) #bcode clustering   
    #redistribute values in row
    for (cID in unique(cbcode)){
      pick=which(cbcode==cID)
      if (length(pick)>1){
        for (j in 1:ncol(X)){
          X[pick,j]=sum(X[pick,j])
        }
      }
    }
    pick=which(!duplicated(cbcode))
    X=X[pick,drop=FALSE,]
    rownames(X)=bcode[pick]
    X_names=paste(colnames(X),collapse=" ")
    #normalise to get crp
    X=t(t(X)/colSums(X))    
    #add up results
    newx=list()
    newx[[X_names]]=X
    CRP=c(CRP,newx)
  }
}

CCRP <- lapply(CRP, ccrpfun)

for(j in 1:length(CRP))
{
X = CRP[[j]]
index = grep('_mut_all',colnames(X))
index_wt = grep('_mut_all',colnames(X),invert=TRUE)
if(length(index)==0)next;
if(length(index)>1&length(index_wt)>1)
{
Xwt = X[,-index];Xwt=ccrpfun(Xwt, clim=50)
Xmut= X[,index];Xmut=ccrpfun(Xmut, clim=50)
X = cbind(Xwt, Xmut)
CCRP[[j]] = X
next;
}
if(length(index)==1&length(index_wt)>1)
{
Xwt = X[,-index];Xwt=ccrpfun(Xwt, clim=50)
Xmut= X[,index]
X1 = cbind(Xwt, X[,index])
colnames(X1)=c(colnames(Xwt),colnames(X)[index])
CCRP[[j]] = X1
next;
}
if(length(index)==1&length(index_wt)==1)
{
CCRP[[j]] = X
}
}



#get txlength
txlength=as.integer(rawmat$RefLength)
names(txlength)=as.character(rawmat$Transcript)
pick=!duplicated(names(txlength))
txlength=txlength[pick]

mutant_isoforms=unique(mygene$TXNAME)

#export CRP to file
save(CRP,CCRP,txlength,mutant_isoforms,file=fout)
cat('\nX matrix is generated\n\n')