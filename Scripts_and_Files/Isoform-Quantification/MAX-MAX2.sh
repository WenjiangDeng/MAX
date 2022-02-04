#!/bin/bash
#SBATCH -A snic2021-5-54
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1-12:00:00 
#SBATCH -J Simulated-RNA-seq-FLT3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wenjiang.deng@ki.se

export LD_LIBRARY_PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1/lib:$LD_LIBRARY_PATH
export PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1/bin:$PATH
MAX_home=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1

module load bioinfo-tools
module load R_packages/3.5.0

# Index the reference.fa in MAX
# TxIndexer -t hg38-flt3-wt-mut.fa -o TxIndexer_idx_hg38 --force

RNA_seq_fasta_dir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML/fasta_file
workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML
output=$workdir

TxIndexer_idx_hg38=/crex/proj/snic2020-6-4/wenjiang/MSE/flt3_2021/x-matirx07-27/TxIndexer_idx_wt_mut_flt3

for fn in $(find  $RNA_seq_fasta_dir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir

MAX -i $TxIndexer_idx_hg38 -l IU -1 $fn  -2 $fn2 -p 16 -o $outdir -w 100000 --strictIntersect TRUE
done

cd $workdir 
## X_matrix.RData should be referenced or saved in $workdir

Rscript $MAX_home/R/Create_count_matrix.R workdir=$workdir core=16
Rscript $MAX_home/R/AEM_update_X_beta.R workdir=$workdir design.matrix=$workdir/X_matrix.RData  core=16 merge.paralogs=FALSE remove.ycount=FALSE saveSubset=FALSE

## Since there are no singletons in FLT3, the MAX estimates are saved as beta.all[[1]] in Beta_final_paralog.Rdata



#!/bin/bash
#SBATCH -A snic2021-5-54
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1-12:00:00
#SBATCH -J Simulated-RNA-seq-NPM1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wenjiang.deng@ki.se

export LD_LIBRARY_PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1/lib:$LD_LIBRARY_PATH
export PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1/bin:$PATH
MAX_home=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/MAX-binary-0.1.1

module load bioinfo-tools
module load R_packages/3.5.0

#TxIndexer -t hg38_mut.fa -o TxIndexer_idx_hg38 --force

RNA_seq_fasta_dir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/fasta_file
workdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML
output=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML

TxIndexer_idx_hg38=/proj/snic2020-6-4/wenjiang/MSE/dbref/npm1/all_mut_2021/100-uniform-bioinfor-paper/TxIndexer_idx_only_npm1

for fn in $(find  $RNA_seq_fasta_dir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir

MAX -i $TxIndexer_idx_hg38 -l IU -1 $fn  -2 $fn2 -p 16 -o $outdir -w 100000 

done

cd $workdir 
Rscript $MAX_home/R/Create_count_matrix.R workdir=$workdir core=16
Rscript $MAX_home/R/AEM_update_X_beta.R workdir=$workdir design.matrix=$workdir/X_matrix.RData  core=16 merge.paralogs=FALSE remove.ycount=FALSE saveSubset=FALSE

###########

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

MAX -i $TxIndexer_idx_hg38 -l IU -1 $fn  -2 $fn2 -p 16 -o $outdir -w 100000 --strictIntersect

done

cd $output ## X_matrix.RData should be referenced or saved in $workdir
Rscript $MAX_home/R/Create_count_matrix.R workdir=$output core=16
Rscript $MAX_home/R/AEM_update_X_beta.R workdir=$output design.matrix=$output/X_matrix.RData  core=16 merge.paralogs=FALSE remove.ycount=FALSE saveSubset=FALSE
