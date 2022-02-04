## Salmon

salmon=/proj/snic2020-6-4/wenjiang/MSE/dbref/npm1/salmon-latest_linux_x86_64/bin/salmon

module load bioinfo-tools
module load R_packages/3.5.0
# Index the reference.fa Salmon
# $salmon index -t hg38_FLT3_mut_wt_clean.fa -i salmon_FLT3_index 

RNA_seq_fasta_dir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/fasta_file

workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/Salmon_FLT3

output=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/Salmon_FLT3

salmon_index_hg38_FLT3=/crex/proj/snic2020-6-4/wenjiang/MSE/flt3_2021/x-matirx07-27/salmon_FLT3_index

for fn in $(find  $RNA_seq_fasta_dir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir

$salmon quant -la -1 $fn -2 $fn2 -i $TxIndexer_idx_hg38 -o $outdir --useEM --hardFilter --numPreAuxModelSamples 0 -p 8

done

## R scripts to process Salmon results
R
options(stringsAsFactors=FALSE)
workdir="/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/Salmon_FLT3"
flist = list.files(paste(workdir,"",sep=""),pattern="quant.sf",recursive=TRUE,full.names = TRUE)
salmon = matrix(0,242,100)

for(i in 1:length(flist))
{
 cat(i,'\n')
 s = read.table(flist[i],header=TRUE) 
 salmon[,i] = s$NumReads
 }
 
rownames(salmon) = s$Name

tmp= s$Name
tmp=sapply(tmp,function(x)strsplit(x,'ut')[[1]][1])
salmon = data.frame(Name=tmp,salmon)
s1 = aggregate(salmon[,-1],by=list(salmon$Name),sum)
salmon = s1
save(salmon,file='salmon_simulate-100-mutated.RData')

######
salmon=/proj/snic2020-6-4/wenjiang/MSE/dbref/npm1/salmon-latest_linux_x86_64/bin/salmon

module load bioinfo-tools
module load R_packages/3.5.0

# $salmon index -t hg38_mut_final.fa -i salmonindex 

fqdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/fasta_file
workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/salmon
output=$workdir

TxIndexer_idx_hg38=/proj/snic2020-6-4/wenjiang/MSE/dbref/npm1/all_mut_2021/100-uniform-bioinfor-paper/salmon_TxIndexer_idx_npm1

for fn in $(find  $fqdir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir

$salmon quant -la -1 $fn -2 $fn2 -i $TxIndexer_idx_hg38 -o $outdir --useEM --hardFilter --numPreAuxModelSamples 0 -p 8

done


R
#workdir 
options(stringsAsFactors=FALSE)
workdir="/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/salmon"
flist = list.files(paste(workdir,"",sep=""),pattern="quant.sf",recursive=TRUE,full.names = TRUE)
salmon = matrix(0,92,100)

for(i in 1:length(flist))
{
 cat(i,'\n')
 s = read.table(flist[i],header=TRUE) 
 salmon[,i] = s$NumReads
 }
 
rownames(salmon) = s$Name

tmp= s$Name
tmp=sapply(tmp,function(x)strsplit(x,'ut')[[1]][1])
salmon = data.frame(Name=tmp,salmon)
s1 = aggregate(salmon[,-1],by=list(salmon$Name),sum)
salmon = s1
save(salmon,file='salmon_simulated-100-mutated.RData')

##############

#Salmon
salmon=/proj/snic2020-6-4/wenjiang/MSE/dbref/npm1/salmon-latest_linux_x86_64/bin/salmon

module load bioinfo-tools
module load R_packages/3.5.0

# $salmon index -t tp53_mut_final.fa2 -i salmonindex 

fqdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/fasta
workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/Salmon
output=$workdir

TxIndexer_idx_hg38=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/ref/salmonindex

for fn in $(find  $fqdir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir

$salmon quant -i $TxIndexer_idx_hg38  --seqBias --gcBias --posBias --biasSpeedSamp 5 -l IU -1 $fn  -2 $fn2 -p 8 -w 100000 --out $outdir

done

R
options(stringsAsFactors=FALSE)
workdir="/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/Salmon"
flist = list.files(paste(workdir,"",sep=""),pattern="quant.sf",recursive=TRUE,full.names = TRUE)
salmon = matrix(0,582,100)

for(i in 1:length(flist))
{
 cat(i,'\n')
 s = read.table(flist[i],header=TRUE) 
 salmon[,i] = s$NumReads
 }
 
rownames(salmon) = s$Name

tmp= s$Name
tmp=sapply(tmp,function(x)strsplit(x,'ut')[[1]][1])
salmon = data.frame(Name=tmp,salmon)
s1 = aggregate(salmon[,-1],by=list(salmon$Name),sum)
salmon = s1
save(salmon,file='salmon_simulate-100-mutated.RData')

