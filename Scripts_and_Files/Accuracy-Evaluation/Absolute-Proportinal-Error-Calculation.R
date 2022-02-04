#####
## R scripts to cacluate the Absolute-Proportinal-Error (APE) using Salmon and MAX estimations
## plot the estimates against true values from simulation
#####
source('APE-source.R')

load('True_counts_for_this_run.RData')
load('salmon_simulate-100-mutated.RData')

s1=salmon
s2=s1[,-1]
rownames(s2) = s1[,1]
s3=t(s2)
colnames(s3) = gsub('_m','_mut_all',colnames(s3))

## compare wt isoforms first
s3 = t(s3)[,c(1:10,12:100,11)]
result_est = s3[grep('mut',rownames(s3),invert=T),]
true1 = t3[1:15,]
Salmon.err = err(result_est,true1)
Salmon.err
## compare mutant isoforms##########################
result_est = s3[grep('mut',rownames(s3)),]
true=t3
true.1 = true[grep('mut',rownames(true)),]
true1 = data.frame(tx=rownames(true.1),true.1)
true1$tx = sapply(true1$tx,function(x)strsplit(x,'_m')[[1]][1])
true2=aggregate(true1[,-1], data.frame(true1$tx), sum)
rownames(true2) = true2[,1]
tr3 = true2[,-1]
tr3 = as.matrix(tr3)
Salmon.err = err(result_est,tr3)
Salmon.err

## draw the plot against true value

png('salmon.png',width=1000,height=1300,res=120)

true1=tr3
par(mfrow=c(4,4), mar=c(4,4,3,1))
for(i in 1:nrow(true1))
{
	plot(log(true1[i,]+0.1),log(result_est[i,]+0.1),xlab="True",ylab="Salmon",main=paste('Salmon ',rownames(true1)[i]),cex.main=1.2,cex.lab=1.2,pch=16)
	abline(0,1,col='brown',lty=1,lwd=1)
}
dev.off()



## MAX
source('APE-source.R')

load('True_counts_for_this_run.RData')
load('Beta_final_paralog.Rdata')
## compare wt isoforms first
est=beta.all[[1]][c(1:10,12:100,11),]
est1=t(est[,grep('mut',colnames(est),invert=T)])
result_est = est1
true = t3[1:15,]
true1 = expand(result_est,true)
MAX.err = err(result_est,true1)
MAX.err

par(mfrow=c(2,2), mar=c(4,4,3,1))
for(i in 1:nrow(true1))
{
	smoothScatter(log(true1[i,]+0.1),log(result_est[i,]+0.1),xlab="True",ylab="MAX",main=rownames(true1)[i],cex.main=0.6)
	abline(0,1,col='brown',lty=2,lwd=2)
}

#######################################################
## compare mutant isoforms##########################


est=beta.all[[1]][c(1:10,12:100,11),]
est1=t(est[,grep('mut',colnames(est))])
rownames(est1)=gsub('_mut_all','',rownames(est1))
result_est = est1

true = t3[grep('mut',rownames(t3)),]
true1 = data.frame(tx=rownames(true),true)
true1$tx = sapply(true1$tx,function(x)strsplit(x,'_m')[[1]][1])
true2=aggregate(true1[,-1], data.frame(true1$tx), sum)
rownames(true2) = true2[,1]
tr3 = true2[,-1]
tr3 = as.matrix(tr3)

true1 = expand(result_est,tr3)
MAX.err = err(result_est,true1)
MAX.err

par(mfrow=c(2,2), mar=c(4,4,3,1))
for(i in 1:nrow(true1))
{
	smoothScatter(log(true1[i,]+0.1),log(result_est[i,]+0.1),xlab="True",ylab="MAX",main=rownames(true1)[i],cex.main=0.6)
	abline(0,1,col='brown',lty=2,lwd=2)
}



