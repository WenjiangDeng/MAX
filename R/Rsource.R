## 06 Feb 2022/ Nghia:
# add extra functions for MAX2
## 28 Jan 2021/ YP
## revise AEM() and estim.BETA() to 
##   (i) allow starting value beta0 for AEM.
##   (ii) report the convergence codes (TRUE/FALSE) of the X and beta estimation
## BUT: Old codes still run!
##
## 16 Nov 2020 / Lu:
# Revised crpcount()
# Added new function hamming_correction()
## 25 Jul 2020 / Nghia:
# revise function ccrpfun() to obtain stable CCRPs by improve the convergence of k-mean 
## 11 Jul 2020 / Nghia:
# revised old and slow crpcount function for 10x genomics data, will improve speed later
## 27 Nov 2019 / Nghia:
# add function getSE() to compute the standard error for the estimates
## Nghia 05 May 2019: 
# - crpcount: vectorisation for a extremely fast computation
## YP 23 April 2018!!
##    revise tcfun(.) not separating the NB1 and removing from the rest,
##    which generated wrong singletons    
##
## YP 25 April 2016: add EMM = modified EM
## 13 Dec 2016: adding estimation by LSP
## 7 Dec 2016: marking the paralogs
## 5 Dec 2016: parallel version of the 4Dec version

## Processing transcript-response profile (TRP), to make
## transcription cluster (TC) and cluster response profile (CRP)
##
## tcfun function:
## input TRP  = the list of trp matrices (must be matrix!!)
## output = transcription clusters
##
## crpfun function:
## input TRP 
##   computes transcription clusters using tcfun
## output= cluster-response profile


## ................................................
## ................................................ tcfun: cluster growing function
## ................................................
tcfun = function(TRP, lim=0.025)
{
# setup neighborhood (NB) list from TRP
# drop likely false neighbors from the TRP if tx-count leakage is too low
NB.n= lapply(TRP,apply,2,sum)  # neighbors with strength=n count
fn = function(x, lim=lim) {  ## keep neighbors with leakage rate >lim
     x = x/max(x); 
     return(unique(names(x[x>lim])))
  }
NB = lapply(NB.n,fn,lim=lim)

## checking singletons, then remove them from the rest
#NBvec= sapply(NB,paste, collapse=' ')
#tst1 = names(NB) == NBvec
#  NB1 = NB[tst1]
#NB2.tmp = NB[!tst1] 
#  fn = function(i) setdiff(NB2.tmp[[i]],NB1)
#  NB2 = foreach(i= 1:length(NB2.tmp)) %dopar% fn(i)  # parallel
#    names(NB2) = names(NB2.tmp)

# TC = union of neighbors of TRPs 
## TC2 = cluster growing from NB
growfn = function(i,NB) sort(unique(unlist(NB[NB[[i]]])))  ## grow from neighborhood NB
system.time(TC2 <- foreach(i= 1:length(NB)) %dopar% growfn(i, NB))
  names(TC2) = names(NB)

# TC3 = grow from TC2
system.time(TC3 <- foreach(i= 1:length(NB)) %dopar% growfn(i, TC2))
TC =  sapply(TC3, paste, collapse=' ') # pasted version, needed for tabulating
  names(TC) = names(NB)
  tcmap = TC  ## map of tx --> TC

#output: list of clusters and tx->cluster map 
clust = names(table(TC))  # clusters
TC = strsplit(clust,split=' ')
  names(TC) = clust

return(list(TC=TC, tcmap = tcmap))
}

## .......................................................
## .......................................................crpfun:  Cluster-response profile
## .......................................................
## input: TRP and TC=tcfun(TRP)$TC
## output: CRP
crpfun = function(TRP, TC)
{
## ALL transcripts
 TX = names(TRP)  
## TC = list of clusters
 clust.size = sapply(TC,length)  ##

### singleton clusters
 C1 = which(clust.size==1)   
 fn = function(x){a= matrix(1); colnames(a) = x ; rownames(a)='1'; return(a)}
 TC1 = lapply(names(C1),fn)
   names(TC1) = names(C1)

## clusters with >1 transcripts:
  TC2 = TC[clust.size>1]
   nc = length(TC2)

## .............upat =  unique occupancy pattern 
## assigning unique patterns based on a fixed-set of transcripts TXc typically from a cluster
## ADD ZERO cols to transcripts missing in trp mat. 
##
upat = function(trp, TXc){
  TX0 = setdiff(TXc, colnames(trp))  # missing tx in the trp
  n0= length(TX0)
  if (n0>0){
     trp0= matrix(0,nrow=nrow(trp),ncol=n0);
     colnames(trp0) = TX0
     trp = cbind(trp, trp0)
  }
  out =  trp[,TXc, drop=FALSE]  ## keep matrix property even for vectors
  out01= ifelse(out>0,1,0)      ## 0-1 pattern: 0=no expression. NOTE: Think about FALSE positive!!
  pat = apply(out01,1,paste,collapse='')
  rownames(out) = pat

  ## collapsing row to unique lines: be careful with 1 unique pattern!!
  nu = length(unique(pat))   ## number of unique pattern
  if (nu>1) out = apply(out,2,tapply,rownames(out),sum)
  if (nu==1) {
      a = matrix(apply(out,2,sum),ncol=ncol(out))
      colnames(a)= colnames(out)
      rownames(a) = unique(pat)
      out = a
  }
  return(out)
}

## for parallel
 incrpfun = function(i){
   TXc = TC2[[i]]  ## set of transcripts in cluster i
     tmp = TXc %in% TX   ## if TX is not complete, TXc may contain external tx
     TXc = TXc[tmp] 
     nt = length(TXc)
   TRPc= TRP[TXc]
   ordtrp = lapply(TRPc, upat, TXc)  ## apply upat to each trp 
   ## all unique patterns
   allpat = unique(unlist(lapply(ordtrp, rownames)))

   ## setup crp for a cluster
   crp = matrix(0, nrow=length(allpat), ncol=nt)
     rownames(crp) = allpat
     colnames(crp) = TXc
   for (nam in TXc){
     rname = rownames(ordtrp[[nam]])
     crp[rname,nam] = ordtrp[[nam]][,nam]
   }
  ## normalize to sum to 1:
  crp = t(t(crp)/colSums(crp))
  return(crp)
}
CRP = foreach(i= 1:length(TC2)) %dopar% incrpfun(i)
names(CRP) = names(TC2)

return(c(TC1,CRP)) ## combine the singletons and others 
}

## ...............................................
## ............................................... ccrpfun: combining paralogs (similar transcripts)
## ...............................................          within each crp
## usage: 
## set.seed(2016)  ## to fix the kmeans clustering results 
## CCRP <- lapply(CRP, ccrpfun, clim=50)
## clim= minimum condition number
##
ccrpfun = function(x, clim=50) {
 set.seed(2016)  ## to fix the kmeans clustering results
 newx = x
 repeat{
   sval = svd(newx)$d;     ## singular values
   #con = max(sval)/sval    ## condition numbers
   con = abs(max(sval)/sval)    ## condition numbers; note from Nghia: might use abs() to avoid the case of -Inf when dividing to zero. This happened randomly, hard to reproduce
   nsc = sum(con< clim)    ## number of clearly separated transcripts
   ## break if OK
   if (nsc== ncol(newx)) break
   ## but collapse if some tx are too similar, using kmeans clustering
   ##   to combine the similar tx's
   ## 24Jul2020/Nghia: increase the iterations and the number of randome sets to get converged results 
   subc = kmeans(t(x),nsc,iter.max = 100000,nstart=100) ## always start with original tx from the crp! 
   newx = t(subc$center)
   clust = subc$cluster
     clust.name = tapply(names(clust), clust, paste, collapse=' ')
     colnames(newx)= clust.name
 }
 return(newx)
}

## ....................................  LS-Pos function
## ....................................  estimation using positive-constrained LS
##  
#require(limSolve)
lsp = function(X,y){
 A = X
 B = y[rownames(X)] ## equalize order of y = rows of X
 p = ncol(X)
 G <- diag(p)
 H <- rep(0,p)
 lsp = lsei(A = A, B = B, G = G, H = H, type=2)
 bt= lsp$X
 err = y- c(X %*% bt)
 L2 = sum(err^2)
 return(list(beta=bt, error=err, L2 = L2))
}

##  
## ............. estimate using EM algorithm
## ............. beta0 = starting value 
##
EM = function(X,y,beta0=rep(sum(y)/ncol(X), ncol(X)), maxiter=50, maxerr=0.01, lim=0.01){
  csum = colSums(X)  ## should be =1, but here just in case
  for (i in 1:maxiter){
    X1 = t(t(X)*beta0)                 ## expected per cell under current est
      rsum = rowSums(X1)
      X2 = X1/ifelse(rsum>0, rsum,1)   ## normalize each row to get probabilities; avoid div by zero
    beta = c(t(X2) %*% y) /csum        ## update:add up info across rows, then normalize by colSum
    err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim)
    if (max(err) < maxerr) break
    beta0 = beta
  }
  return(beta)
}

##
## modify=TRUE: modify beta with digamma function
##
EMM = function(X,y,beta0=rep(sum(y)/ncol(X), ncol(X)), maxiter=100, 
               maxerr=0.01, lim=0.01, modify=TRUE){
   csum = colSums(X)  ## should be =1, but here just in case
                      ## sailfish X version, is row-normalized, so does not add-up to one
   logNorm = digamma(sum(y))
   for (i in 1:maxiter){
    gamma0= beta0
    if (modify) gamma0=exp(digamma(beta0+0.00000001) - logNorm)
    # beta0=ifelse(is.na(beta0),0,beta0)
    #print(c(i,beta0, sum(beta0)))
    X1 = t(t(X)*gamma0)                 ## expected per cell under current est
    rsum = rowSums(X1)
    X2 = X1/ifelse(rsum>0, rsum,1)    ## normalize each row to get probabilities; avoid div by zero
    beta = c(t(X2) %*% y) # /csum        ## update:add up info across rows, then normalize by colSum
    err = abs(beta-beta0)/ifelse(beta0>lim, beta0, lim)
    if (max(err) < maxerr) break
    beta0 = beta
    }
   return(beta)
 }

###crpcount function - for fast computation -but not suitable for 10x genomics
#crpcount = function(filename, weight_thres=0.9)
#{
#  options(stringsAsFactors=FALSE)
#  rawcount = read.table(filename,header=T)
#  sample.name1 = strsplit(filename,"/")[[1]]
#  sample.name = sample.name1[length(sample.name1)-1]
#
#  ## defining eqClasses in the CRPs
#  ## all CRPs are assumed matrices -- even 1x1 singletons
#  ##  *** Tx columns within CRP are sorted***
#  crp.eqfun = function(crp){
#    tx = colnames(crp)
#    txmat = outer(rep('',nrow(crp)), tx, paste0)
#    upat = matrix('', nrow=nrow(crp), ncol=ncol(crp))
#    upat[crp>0] = txmat[crp>0]
#    out = apply(upat,1,paste, collapse='')
#    return(out)
#  }
#  ## full list eqClasses from CRPs: 
#  CRP.eq <- lapply(CRP, crp.eqfun)# 15 Otc 2019/Nghia: use lapply instead of sapply, this avoids error when having only 1 CRP
#  CRP.names = rep(names(CRP.eq), lengths(CRP.eq))
#  CRP.eqID = paste(CRP.names, unlist(CRP.eq)) 
#  ## check: each eqClass is unique, no repeat
#  ## length(unique(unlist(CRP.eq))); length(unlist(CRP.eq))
  #
#  #map from crp to tx
#  c2t=sapply(names(CRP), strsplit, split=' ')
#  cLen=lengths(c2t)
#  cVec=rep(names(cLen),cLen)
#  tVec=unlist(c2t)  
#  names(cVec) = tVec
#  rawcount$crp = cVec[rawcount$Transcript]
  #
#  #do the weight filter here
#  matchID=match(rawcount$crp,names(cLen))
#  rawcount$eqlen=cLen[matchID]  
#  pick = (rawcount$eqlen>1)| 
#    (rawcount$eqlen==1 & rawcount$Weight>weight_thres )
#  # pick = rep(TRUE, nrow(rawcount))
  #
#  rawcount$eqCrp=paste(rawcount$crp,rawcount$eqClass)
#  eqClass = tapply(rawcount$Transcript[pick],rawcount$eqCrp[pick],
#                   function(x){paste(sort(x),collapse='')})
#  eqCount=tapply(rawcount$Count[pick],rawcount$eqCrp[pick],min)
#  eqCRP=tapply(rawcount$crp[pick],rawcount$eqCrp[pick],unique)
#  eqID = paste(eqCRP, eqClass)
#  eqCount2=tapply(eqCount,eqID,sum)
#  eqID2=names(eqCount2)
  #
#  ##  .. transfer data to CRP objects
#  ycount = rep(0, length(CRP.eqID))
#  names(ycount) = CRP.eqID
#  pick = eqID2 %in% CRP.eqID
#  ycount[eqID2[pick]]= eqCount2[pick]
  #
#  ## repack into CRPs
#  CRP.y = tapply(ycount, CRP.names, c)
#  out = Map(cbind, CRP, sample1=CRP.y[names(CRP)])
  #
#return(list(crpcount=out,samplename=sample.name))
#}

##crpcount function - orignal function from Wenjiang
crpcount = function(filename)
{
 options(stringsAsFactors=FALSE)
 rawcount = read.table(filename,header=T)
 
 sample.name1 = strsplit(filename,"/")[[1]]
 sample.name = sample.name1[length(sample.name1)-1]

 #remove zero eqc
 pick=which(rawcount$Count !=0)
 rawmat=rawcount[pick,]
 #reindex
 myeqID=unique(rawcount$eqClass)
 myeqID=sort(myeqID)
 matchID=match(rawcount$eqClass,myeqID)
 rawcount$eqClass=matchID
 ###

 tx1 = unique(rawcount$Transcript)
 nID=c(1:length(CRP))
 rapmap_crp=NULL
 extra_pattern = NULL
 for (j in nID)
 {
 #cat(j,'\n')
 crp = CRP[[j]]
 if(length(crp)==1)
 {
 tx = colnames(crp)
 sample1 = sum(rawcount[rawcount$Transcript==tx & rawcount$Weight>0.9,"Count"])
 crp1 = cbind(crp,sample1)
 crp.y = crp1
 }else{
   tx = colnames(crp)
   raw.c =NULL
   pat=NULL
   for(i in 1:length(tx))
      raw.c = rbind(raw.c,rawcount[rawcount$Transcript==tx[i],])
   eq = unique(raw.c$eqClass)
   if(length(eq)>0)
 {
 raw2 = matrix(0,length(eq),length(tx))
 colnames(raw2) = tx
 rownames(raw2) = eq
 for(z in 1:nrow(raw.c))
 raw2[as.character(raw.c[z,grep("eqClass", colnames(raw.c), ignore.case = T )]),raw.c[z,1]] = raw.c[z,grep("Count", colnames(raw.c), ignore.case = T )]##3 is counts, 6 is eqclass

 raw3 = ifelse(raw2>0,1,0)
 pat = apply(raw3,1,function(x)paste(x,collapse=""))
 raw3 = data.frame(pat=pat,raw3)
 raw2.1 = data.frame(pat=pat,raw2)
 a1 = apply(raw2,1,max)
 a2 = tapply(a1,pat,sum)
 pat=unique(pat)
 raw4 = data.frame(pat=names(a2),sample1=a2)
 crp1 = data.frame(pat=rownames(crp),crp)

 raw4 <- hamming_correction(raw4, crp1, hamming_dist = 1)
 
 m = merge(crp1,raw4,by="pat",all.x=T,sort=F)
 rownames(m) = as.character(m$pat)
 m[is.na(m)] = 0
 m = m[,grep("pat", colnames(m), ignore.case = T, invert = T)]
 crp.y = m[rownames(crp),]
 
 }else{
 crp.y = data.frame(crp,sample1=0)
}

}
crp.y=as.matrix(crp.y)
if(j ==1) rapmap_crp = list(crp.y)
if(j >1) rapmap_crp[[j]] = crp.y
}
names(rapmap_crp) = names(CRP)
return(list(crpcount=rapmap_crp,samplename=sample.name))
}

# 2020-11-16:
# Function added by: Lu, to fix the problems with no comparable equivalence class
hamming_correction <- function(raw_data, crp_data, hamming_dist = 1){   # by default, hamming correction of distance = 1
  raw_current=lapply(raw_data$pat, function(x) as.integer(unlist(strsplit(x,""))))
  crp_current=lapply(crp_data$pat, function(x) as.integer(unlist(strsplit(x,""))))
  
  #compute hamming distance
  res=lapply(raw_current, function(x){
    z=sapply(crp_current,function(y) sum(x!=y))
    return(z) 
  })
  res=do.call(rbind,res)
  dim(res)
  rownames(res)=raw_data$pat
  colnames(res)=crp_data$pat
  res=as.matrix(res)
  compare_hamming=res

#  raw_current <- data.frame(t(data.frame(strsplit(raw_data$pat, split = ""))))
#  for(i in 1:ncol(raw_current)){
#    raw_current[,i] <- as.numeric(as.character(raw_current[,i]))
#  }
#  row.names(raw_current) <- raw_data$pat
  #
#  crp_current <- data.frame(t(data.frame(strsplit(crp_data$pat, split = ""))))
#  for(i in 1:ncol(crp_current)){
#    crp_current[,i] <- as.numeric(as.character(crp_current[,i]))
#  }
#  row.names(crp_current) <- crp_data$pat
  #
#  compare_hamming <- NULL
#  for(p in 1:nrow(raw_current)){
#    temp <- NULL
#    for(q in 1:nrow(crp_current)){
#      current <- length(which(abs(crp_current[q,] - raw_current[p,]) > 0))
#      temp <- cbind(temp, current)
#    }
#    compare_hamming <- rbind(compare_hamming, temp)
#  }
  #
#  compare_hamming <- data.frame(compare_hamming)
#  colnames(compare_hamming) <- row.names(crp_current)
#  row.names(compare_hamming) <- row.names(raw_current)

  raw_corrected <- raw_data
  # current <- apply(compare_hamming,1,function(x){which(x == 1)})
  for(i in 1:nrow(compare_hamming)){

    #if this eqclass is found in crp
    if(sum(compare_hamming[i,] == 0) > 0){ 
      raw_corrected$pat[i] <- colnames(compare_hamming)[which(compare_hamming[i,] == 0)]
    }else{
      closeNum=sum(compare_hamming[i,] <= hamming_dist) #number of eqclasses with hamming_dist less than the threshold
      if( closeNum > 1){ #more than one close-eqclass
        # Lu noted here: If we choose only the one with highest value, then might be a bias towards high value in real data
        # might affect those with real cases of low values
        current <- crp_data[which(compare_hamming[i,] <= hamming_dist),]
        raw_corrected$pat[i] <- row.names(current)[which.max(apply(current,1,function(x){median(as.numeric(as.character(x[x>0])))}))]
      }

      if(closeNum == 1){ #only one-close eqclass
        raw_corrected$pat[i] <- colnames(compare_hamming)[which(compare_hamming[i,] <= hamming_dist)]
      }

      if(closeNum == 0){ #no close-eqclass
        current <- crp_data[which(compare_hamming[i,] <= min(compare_hamming[i,])),]
        raw_corrected$pat[i] <- row.names(current)[which.max(apply(current,1,function(x){median(as.numeric(as.character(x[x>0])))}))]
      }
    }
  }
  
#  if(sum(duplicated(raw_corrected$pat)) > 0){
  #
  #current <- raw_corrected[!(!(duplicated(raw_corrected$pat) | duplicated(raw_corrected$pat, fromLast = TRUE))),]
  #raw_corrected <- raw_corrected[(!(duplicated(raw_corrected$pat) | duplicated(raw_corrected$pat, fromLast = TRUE))),]
#  current <- split(current, current$pat)
  ## 2020-11-17: Modified by Lu, changed from select from the one that has max median to include all
  ## current <- lapply(current, function(x){x <- x[which.max(apply(data.frame(x[,grep("pat", colnames(x), ignore.case = T, invert = T)]),1,median)),]})
  #current <- lapply(current, function(x){x <- data.frame(pat = unique(x$pat), sample1 = sum(x$sample1))})
#  current <- do.call(rbind.data.frame, current)
#  raw_corrected <- data.frame(rbind(raw_corrected, current))
#  }
  
   a = tapply(raw_corrected$sample1,raw_corrected$pat,sum)
   raw_corrected = data.frame(pat=names(a),sample1=a)
 
  row.names(raw_corrected) <- raw_corrected$pat
  return(raw_corrected)
}

## 25 April 2017: 
## .. YP: estim.BETA and AEM allow modify=TRUE for digamma modification
##
## 18 April 2017
## ..  estim.BETA does not use big matrix Z  
##
## joint estimation of A and theta for sequgio model
## SEE sequgio-model and estimation.r
##

## to use LSP= positive-constrained least-squares from library(limSolve)
#require(limSolve)

##  ......................  Alternating EM algorithm applied to one CRP
##
##   X0   = starting value for design matrix X = CRP
##   Ymat = data matrix, ncol(Ymat)= sample size n; nrow(Ymat) = num of eqClasses
##   beta0 = starting value of beta; if missing, beta0 is set to equal expression
##
AEM = function(X0, Ymat, beta0, maxiter.X=100, maxerr=0.01, lim=0.01,
                         maxiter.bt=100, modify=TRUE){
 X.est0 = X0  # starting design matrix
 Ysum= colSums(Ymat)
 if (!missing(beta0)) beta=beta0
 #  default starting espression = average expression
 if (missing(beta0)) beta = rep(Ysum/ncol(X0), rep(ncol(X0), length(Ysum)))

 ## alternating estimation
 X.conv= beta.conv= FALSE   ## convergence codes
 for (i in 1:maxiter.X){
   BT = estim.BETA(X.est0, Ymat, beta0=beta, maxiter=maxiter.bt, modify=modify); 
   beta.conv = attr(BT,'conv')  ## 
   X.est = estim.X.em(X.est0,BT,Ymat); ## print(rbind(BT,X.est))
   ## check convergence
   err = abs(X.est-X.est0)/ifelse(abs(X.est0)>lim, X.est0, lim)  ## X.est0 could be zero
   if (max(err) < maxerr) {X.conv=TRUE; break}
   ## update for next iterate
   X.est0 = X.est
   beta=c(t(BT))   ## NOTE: BT is a matrix, saved byrow=T
 }
return(list(X = X.est, BETA=BT, beta.conv=beta.conv, X.conv=X.conv))
}


##
## EM without for-loop over sample size n, only over p
##
estim.X.em = function(X.em,BETA,Ymat){
  X0 = X.em  ## starting value
  Y = c(Ymat)
  q = nrow(Ymat)   ## = nrow(X.em)      ## number of eqClasses
  n = ncol(Ymat)   ## = nrow(BETA)     ## sample size
  p = ncol(BETA)  ## number of transcripts
 BTsum = rep(colSums(BETA), rep(q,p)) ## NOTE: BTsum might be zero!!
 x1 = NULL
 for (j in 1:p) x1=cbind(x1,kronecker(BETA[,j],X.em[,j]))
 sumx1 = rowSums(x1) 
   sumx1.pos= ifelse(sumx1>0, sumx1, 0.1)
 x2 = x1/sumx1.pos
 Ycell= x2*Y
 tmp = NULL 
 for (j in 1:p){
   tmp = c(tmp,colSums(matrix(Ycell[,j], nrow=n,byrow=TRUE)))
  }
 ## NOTE BTsum might be zero!!
 ## if csum=0, return to the original X0
 X.em = matrix(tmp/ifelse(BTsum>0, BTsum, 0.1), ncol=p)
   csum = colSums(X.em)
   X.em = t(t(X.em)/ifelse(csum>0,csum,0.1))  ## normalize to sum to 1
 for (j in 1:length(csum)){if (csum[j]< 0.00001) X.em[,j]=X0[,j]}
 return(X.em)
}



##
## .. modify=TRUE: EM modified by the digamma function
## .. only EM and WITHOUT big matrix Z
## .. given (crp X, data matrix Ymat): estimate Tx expression-matrix BETA
## .. Ymat size = qxn matrix: sample = columns
## .. beta0 vector of starting value for BETA: length = nxp
##
estim.BETA = function(X,Ymat,beta0, maxiter=100, maxerr=0.01, 
                         lim=0.01, modify=TRUE){  
  Y.list =  split(t(Ymat), seq(ncol(Ymat))) 
  q = nrow(X)    ## number of quasi-exons
  n = ncol(Ymat) ## number of samples
  csum = colSums(X)
  totCount= colSums(Ymat)
  logNorm = digamma(totCount+0.00000001)

  ## computing the old X2, allow modify=TRUE
  fn1 = function(bt, logNorm.i){
    gamma0= bt
    if (modify) gamma0=exp(digamma(bt+0.00000001) - logNorm.i)
    tmp = t(t(X)*gamma0)    ## NOTE: X is fixed, taken from the input
    rsum = rowSums(tmp)
    return(tmp/ifelse(rsum>0, rsum,1))
  }
  ## computing each (Xi'yi)
  fn2 = function(Xi,yi){
    bt = c(t(Xi) %*% yi) /csum  ## normalized by colSums(X)
    return(bt)
  }

beta.conv=FALSE
for (i in 1:maxiter){
  BETA = matrix(beta0, nrow=n, byrow=TRUE)
  BT.list = split(BETA, seq(nrow(BETA))) ## list version
  X2.list = mapply(fn1, BT.list, logNorm, SIMPLIFY=FALSE)
  beta= mapply(fn2, X2.list, Y.list, SIMPLIFY=TRUE)
  err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim); #print(max(err))
  if (max(err) < maxerr) {beta.conv=TRUE; break}
    beta0 = beta
}  ## endfor i

beta=matrix(beta, nrow=n, byrow=TRUE)  ## output
attr(beta,'conv') = beta.conv
return(beta)
}


## .. older version: EM with no modification 
## .. only EM and WITHOUT big matrix Z
## .. given (crp X, data matrix Ymat): estimate Tx expression-matrix BETA
## .. Ymat size = qxn matrix: sample = columns
## .. beta0 vector of starting value for BETA: length = nxp
##
estim.BETA.M0 = function(X,Ymat,beta0, maxiter=100, maxerr=0.01, lim=0.01){  
  Y.list =  split(t(Ymat), seq(ncol(Ymat))) 
  q = nrow(X)    ## number of quasi-exons
  n = ncol(Ymat) ## number of samples
  csum = colSums(X)

  ## computing the old X2 
  fn1 = function(bt){
    tmp = t(t(X)*bt)    ## NOTE: X is fixed, taken from the input
    rsum = rowSums(tmp)
    return(tmp/ifelse(rsum>0, rsum,1))
  }
  ## computing each (Xi'yi)
  fn2 = function(Xi,yi){
    bt = c(t(Xi) %*% yi) /csum  ## normalized by colSums(X)
    return(bt)
  }

for (i in 1:maxiter){
  BETA = matrix(beta0, nrow=n, byrow=TRUE)
  BT.list = split(BETA, seq(nrow(BETA))) ## list version
  X2.list = Map(fn1,BT.list)
  beta= mapply(fn2,X2.list,Y.list)
  err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim); #print(max(err))
  if (max(err) < maxerr) break
    beta0 = beta
}  ## endfor i

return(matrix(beta, nrow=n, byrow=TRUE))
}


##
## .. using potentially big matrix Z
## .. given (crp X, data matrix Ymat): estimate Tx expression-matrix BETA
## .. beta0 = starting value for BETA
##
estim.BETA.Z = function(X,Ymat,beta0, method='EM', maxiter=100, maxerr=0.01, lim=0.01){  
  Y = c(Ymat)
  q = nrow(X)   ## number of eqClasses
  n = ncol(Ymat)
  Z = diag(n) %x% X   ## kronecker product
 if (method=='OLS') {
    ZZ = t(Z)%*%Z
    bt.est = solve(ZZ, t(Z) %*%Y)
 }
 if (method=='LSP') {
   G = diag(ncol(Z))
   H = rep(0,ncol(Z))
   bt.est = lsei(A=Z,B=Y,G=G, H=H, verbose=FALSE)$X
 }
 if (method=='EM'){
    csum = colSums(Z)  ## should be =1, but here just in case
    for (i in 1:maxiter){
      X1 = t(t(Z)*beta0)                ## expected per cell under current est
       rsum = rowSums(X1)
      X2 = X1/ifelse(rsum>0, rsum,1)   ## normalize each row to get probabilities; avoid div by zero
     beta = c(t(X2) %*% Y) /csum        ## update:add up info across rows, then normalize by colSum
     err =  abs(beta-beta0)/ifelse(beta0>lim, beta0, lim); #print(max(err))
     if (max(err) < maxerr) break
     beta0 = beta
    }  ## end for i
  bt.est = beta
  } ## endif method=EM
 return(matrix(bt.est, nrow=n, byrow=TRUE))
}

## .. given (Tx expression, data Ymat, sample size n=ncol(Ymat)): estimate X
## .............. OLS: note Z is potentially a big matrix!!
##
estim.X = function(BETA,Ymat,method='LSP'){
 Y = c(Ymat)
 n = ncol(Ymat)    ## also = nrow(BETA)
 q = nrow(Ymat)    ## number of eqClasses per sample
 p = ncol(BETA)   ## number of transcripts
 sdy= sd(Y)
   Y1= Y/sdy       ## standardize to make lsp work better!!
 Z = kronecker(BETA, diag(q))
   Z = Z/sdy  ## normalize size
 if (method=='OLS'){
   ZZ = t(Z)%*%Z
   X = solve(ZZ, t(Z) %*%Y1)
 }
 if (method=='LSP'){
  G = diag(ncol(Z))
  H = rep(0,ncol(Z))
  X = lsei(A=Z,B=Y1,G=G, H=H, verbose=FALSE)$X
 }
 # normalize
 tmp = matrix(X, ncol=p)    
   X = t(t(tmp)/colSums(tmp))
 return(X)
}

## ................  EM based
## 
## calculate without using big matrix Z, but use for-loop
## the background pattern X.em is in matrix form
## .. given BETA= matrix of Tx expression, ** samples in rows ***
## Y = data vector, sample by sample
## n  = sample size n 
## ==> update X
#= X.em= matrix of starting values with OLS+
estim.X.em0 = function(X.em,BETA,Ymat){
  Y = c(Ymat)
  q = nrow(Ymat)   ## = nrow(X)      ## number of eqClasses
  n = ncol(Ymat)   ## = nrow(BETA)  ## sample size
  p = ncol(BETA)  ## number of transcripts
 Ymat = matrix(Y,ncol=n)
 Ycell = NULL
 for (j in 1:n){
  bt = BETA[j,]                    ## expression-vector of j'th sample
  y = Ymat[,j]                   ## data vector from j'th sample
  x1 = t(t(X.em)*bt)             ## expected per cell
    x2 = x1/rowSums(x1)          ## normalize each row to get probabilities
    ycell = c(x2* y) 
    Ycell = rbind(Ycell, ycell)
  }
  BTsum = rep(colSums(BETA), rep(q,p))
  X.em = matrix(colSums(Ycell)/BTsum, ncol=p)
    X.em = t(t(X.em)/colSums(X.em))  ## normalize to sum to 1
  # print(round(c(X.em),2))
  return(X.em)
}

##
## EM without for-loop over sample size n, only over p
##
estim.X.em = function(X.em,BETA,Ymat){
  X0 = X.em  ## starting value
  Y = c(Ymat)
  q = nrow(Ymat)   ## = nrow(X.em)      ## number of eqClasses
  n = ncol(Ymat)   ## = nrow(BETA)     ## sample size
  p = ncol(BETA)  ## number of transcripts
 BTsum = rep(colSums(BETA), rep(q,p)) ## NOTE: BTsum might be zero!!
 x1 = NULL
 for (j in 1:p) x1=cbind(x1,kronecker(BETA[,j],X.em[,j]))
 sumx1 = rowSums(x1) 
   sumx1.pos= ifelse(sumx1>0, sumx1, 0.1)
 x2 = x1/sumx1.pos
 Ycell= x2*Y
 tmp = NULL 
 for (j in 1:p){
   tmp = c(tmp,colSums(matrix(Ycell[,j], nrow=n,byrow=TRUE)))
  }
 ## NOTE BTsum might be zero!!
 ## if csum=0, return to the original X0
 X.em = matrix(tmp/ifelse(BTsum>0, BTsum, 0.1), ncol=p)
   csum = colSums(X.em)
   X.em = t(t(X.em)/ifelse(csum>0,csum,0.1))  ## normalize to sum to 1
 for (j in 1:length(csum)){if (csum[j]< 0.00001) X.em[,j]=X0[,j]}
 return(X.em)
}

getSE <-function(X,b,y){
  #general for a matrix of y
    yhat=X %*% t(b)
    err=y-yhat
    W = 1/(yhat+1)  ## weight matrix
    SE = NULL
    for (j in 1:ncol(W)){
      WX = (W[,j]*X)
      XWX = t(X) %*% WX
      XWXi = solve(XWX)
      XWeWX = t(WX) %*% (err[,j]^2 *WX)
      mat = XWXi %*% XWeWX %*% XWXi
      se = sqrt(diag(mat))
      SE = rbind(SE,se)
    }
    return(SE)
}

###### MAX2

buildMutantIso<-function(trueCount){
  #merge all mutant isoforms into a single mutant isoform
  pick=grep("mut",rownames(trueCount))
  wtIso=rownames(trueCount)[-pick]
  trueCountOut=trueCount[-pick,]
  mutCount=NULL
  for (tx in wtIso){
   pick1= grep(tx,rownames(trueCount))
   pick2=intersect(pick,pick1)
   countVal=colSums(trueCount[pick2,drop=FALSE,])
   mutCount=rbind(mutCount,countVal)
  }
  rownames(mutCount)=paste0(wtIso,"_mut_all")
  trueCountOut=rbind(trueCountOut,mutCount)
  trueCountOut=trueCountOut[order(rownames(trueCountOut)),]
  return(trueCountOut)
}

buildParalogIso<-function(countDat,paralogs){
  #merge isoforms to paralogs
  paraCount=NULL
  for (p in paralogs){
    tx=trimws(unlist(strsplit(p," ")))
    pick=rownames(countDat) %in% tx
    countVal=colSums(countDat[pick,drop=FALSE,])
    paraCount=rbind(paraCount,countVal)
  }
  rownames(paraCount)=paralogs
  return(paraCount)
}

mergeMut<-function(a,tt){

options(StringsAsFactors=F)
a_merge = NULL
for(i in unique(a$eqClass)){
  b =  NULL
   
  b = subset(a,eqClass==i)
  
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
  
rownames(a_merge) = c(1:nrow(a_merge))
x = a_merge
x$Transcript = gsub("_mut_[[:print:]]+$",'_mut_all',x$Transcript)#letters, numbers, punctuation, and whitespace.
#x$Transcript = gsub("_mut_[[:alnum:]]+$",'_mut_all',x$Transcript)
#x$Transcript = gsub("_mut_[0-9]{1,}",'_mut_all',x$Transcript)
#unique(x$Transcript)


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

x_merge = NULL
for(i in unique(x$eqClass)){ 
  x1 = NULL
  x1 = subset(x,eqClass==i)
  x1 = x1[order(x1$Transcript),]
  x1 = cbind(x1,Isoform=paste(x1$Transcript,collapse='_'))
  x_merge = rbind(x_merge,x1)
}



z_merge = NULL
allIsoforms=unique(x_merge$Isoform)
for(i in 1:length(allIsoforms)){ 
  b = NULL
  b = subset(x_merge,Isoform==allIsoforms[i])[,-7]
  b1 = rbind(tapply(b$Weight,b$Transcript,sum),tapply(b$Count,b$Transcript,sum),tapply(b$EffLength,b$Transcript,mean),tapply(b$RefLength,b$Transcript,mean))
  b1 = t(b1)
  b2 = data.frame(Transcript=rownames(b1),b1,eqClass=i)
  colnames(b2) = colnames(b)
  z_merge = rbind(z_merge,b2)
}
rownames(z_merge) = c(1:nrow(z_merge))

return(z_merge)
}



getCRP <- function(rawmat,H_thres=0.025,ncore=1){

  library(foreach)
  library(doParallel)
  registerDoParallel(ncore)

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

  ccrpfun = function(x, clim=50) {
   set.seed(2016)  ## to fix the kmeans clustering results
   newx = x
   repeat{
     sval = svd(newx)$d;     ## singular values
     #con = max(sval)/sval    ## condition numbers
     con = abs(max(sval)/sval)    ## condition numbers; note from Nghia: might use abs() to avoid the case of -Inf when dividing to zero. This happened randomly, hard to reproduce
     nsc = sum(con< clim)    ## number of clearly separated transcripts
     ## break if OK
     if (nsc== ncol(newx)) break
     ## but collapse if some tx are too similar, using kmeans clustering
     ##   to combine the similar tx's
     ## 24Jul2020/Nghia: increase the iterations and the number of randome sets to get converged results 
     subc = kmeans(t(x),nsc,iter.max = 100000,nstart=100) ## always start with original tx from the crp! 
     newx = t(subc$center)
     clust = subc$cluster
       clust.name = tapply(names(clust), clust, paste, collapse=' ')
       colnames(newx)= clust.name
   }
   return(newx)
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

  #export CRP to file
  #save(CRP,CCRP,txlength,file=fout)
  res=list(CRP=CRP,CCRP=CCRP,CRPCOUNT=CRPCOUNT,txlength=txlength)
  return(res)
}

