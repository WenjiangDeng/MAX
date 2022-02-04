
	
err = function(result_est,true){##sailfish2 and true are expanded

	seq_err = matrix(0,nrow(result_est),1)
	rownames(seq_err) = rownames(result_est)

	for(i in 1:nrow(result_est))
{
	seq_err2 = abs(result_est[i,] - true[i,])/(true[i,]+1)
	seq_err[i] = median(seq_err2,na.rm=TRUE)
}
	return(seq_err)
}


expand = function(result_est,Sailfish){##result_est from sequgio , Sailfish from different methods
sailfish1 = matrix(0,length(unique(unlist(strsplit(rownames(result_est)," ")))),100)
rownames(sailfish1) = unique(unlist(strsplit(rownames(result_est)," ")))
colnames(sailfish1) = colnames(result_est)
sailfish1[rownames(sailfish1),] = Sailfish[rownames(sailfish1),]

sailfish2 = matrix(0,nrow(result_est),100)
rownames(sailfish2) = rownames(result_est)
colnames(sailfish2) = colnames(result_est)
rowname = rownames(result_est)
tx = strsplit(rownames(result_est)," ")
ntx = sapply(tx,length)
rowname_sailfish = rownames(Sailfish)[ rownames(Sailfish) %in% rowname]
sailfish2[rowname_sailfish,] = Sailfish[rowname_sailfish,]
pick = ntx>1
tx2 = tx[pick]
for(i in 1:length(tx2)){
 tmp = sailfish1[tx2[[i]],];tmp2 = colSums(tmp);sailfish2[paste(tx2[[i]],collapse=" "),] = tmp2;
 }
## sailfish2 
return(sailfish2)
}
