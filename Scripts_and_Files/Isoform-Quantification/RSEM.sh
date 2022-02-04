#!/bin/bash
#SBATCH -A snic2021-5-54
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1-12:00:00 
#SBATCH -J Simulated-RNA-seq-FLT3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wenjiang.deng@ki.se


module load bioinfo-tools
module load R_packages/3.5.0

RNA_seq_fasta_dir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML/fasta_file
workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML
output=$workdir

# RSEM 
# build index
# rsem-prepare-reference --bowtie2  hg38_FLT3_mut_wt_clean.fa hg38_FLT3 -p 16
# quantification 
rsem-calculate-expression -p 16 --paired-end --bowtie2 --estimate-rspd \
--no-qualities --bowtie2-k 10000 --no-bam-output $fn $fn2 \
$indexdir/bowtie2_index/hg38_FLT3 $ind 
done


## R scripts to process RSEM results
R
options(stringsAsFactors=FALSE)
workdir="/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/RSEM_FLT3"
flist = list.files(workdir,pattern=".fasta.isoforms.results",recursive=TRUE,full.names = F)
exp = NULL
sname = NULL
for(i in 1:length(flist))
{
 exp1 = read.table(flist[i],header=TRUE)
 sname = c(sname,gsub('_1.fasta.isoforms.results','',flist[i]))
 exp = cbind(exp,exp1[,5])
}
rownames(exp)=exp1[,1]
colnames(exp)=sname		
save(exp,file='RSEM_est_FLT3.RData')

######

#!/bin/bash
#SBATCH -A snic2021-5-54
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1-12:00:00
#SBATCH -J Simulated-RNA-seq-NPM1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wenjiang.deng@ki.se

module load bioinfo-tools
module load R_packages/3.5.0

RNA_seq_fasta_dir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/fasta_file
workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML
output=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML


for fn in $(find  $RNA_seq_fasta_dir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir

rsem-calculate-expression -p 16 --paired-end --bowtie2 --estimate-rspd --no-qualities --bowtie2-k 10000 --no-bam-output $fn $fn2 \
/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/rsem/NPM1_index/npm1_index $ind 

done

R
options(stringsAsFactors=FALSE)
workdir="/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/rsem"
flist = list.files(workdir,pattern=".fasta.isoforms.results",recursive=TRUE,full.names = F)
exp = NULL
sname = NULL
for(i in 1:length(flist))
{
 exp1 = read.table(flist[i],header=TRUE)
 sname = c(sname,gsub('_1.fasta.isoforms.results','',flist[i]))
 exp = cbind(exp,exp1[,5])
}
rownames(exp)=exp1[,1]
colnames(exp)=sname		
save(exp,file='RSEM_est_NPM1_simul.RData')

#####

#!/bin/bash
#SBATCH -A snic2021-5-54
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 6:00:00
#SBATCH -J Simulated-data-TP53
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wenjiang.deng@ki.se

export LD_LIBRARY_PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1/lib:$LD_LIBRARY_PATH
export PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1/bin:$PATH
MAX_home=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1

module load bioinfo-tools
module load R_packages/3.5.0

#TxIndexer -t hg38_final.fa -o TxIndexer_idx_hg38 
RNA_seq_fasta_dir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/fasta
workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/MAX
output=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/MAX

TxIndexer_idx_hg38=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/ref/TxIndexer_MAX

for fn in $(find  $RNA_seq_fasta_dir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir

# RSEM 
# build index
rsem-calculate-expression -p 16 --paired-end --bowtie2 --estimate-rspd --no-qualities --bowtie2-k 10000 --no-bam-output $fn $fn2 \
/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/ref/rsem_index/tp53_rsem_index $ind 
done


R
options(stringsAsFactors=FALSE)
workdir="/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/RSEM"
flist = list.files(workdir,pattern=".fasta.isoforms.results",recursive=TRUE,full.names = F)
exp = NULL
sname = NULL
for(i in 1:length(flist))
{
 exp1 = read.table(flist[i],header=TRUE)
 sname = c(sname,gsub('_1.fasta.isoforms.results','',flist[i]))
 exp = cbind(exp,exp1[,5])
}
rownames(exp)=exp1[,1]
colnames(exp)=sname		
save(exp,file='RSEM_est_tp53_simul.RData')
