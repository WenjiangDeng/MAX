## working directory on UPPMAX
## /crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/One_SNV_simul_WT_ref_XAEM/ASEP

#!/bin/bash
#SBATCH -A snic2021-5-54
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 5:00:00
#SBATCH -J BeatAML-data-npm1-all-mutation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wenjiang.deng@ki.se

module load bioinfo-tools 
module load  star/2.7.9a

module load samtools/1.14
module load VarScan/2.4.2

#FLT3 mapping
fd=/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/reference
fqdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/One_SNV_simul/fasta
# Nghia folder
# /proj/snic2020-6-4/nobackup/Nghia/Working/ASE_isoform/revision1/Data_Reference_MAX
wkdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/One_SNV_simul_WT_ref_XAEM/ASEP
output=$wkdir
cd $wkdir

for fn in $(find  $fqdir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output/$(basename $fn))
sample=$(basename $fn)
echo $sample
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir
mkdir -p $outdir
cd $outdir

STAR --genomeDir $fd/chr13_hg19_index \
--runThreadN 8 \
--readFilesIn $fn $fn2 \
--outFileNamePrefix $outdir/$sample"_" \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

samtools index  $outdir/$sample"_Aligned.sortedByCoord.out.bam"
 
samtools mpileup -r chr13:28577411-28674713 -f $fd/hg19_chr13.fa $outdir/$sample"_Aligned.sortedByCoord.out.bam" | java -jar  /crex/proj/snic2020-6-4/wenjiang/MSE/revision2/software/VarScan.v2.3.9.jar mpileup2snp --min-coverage 1  --min-avg-qual 15 --min-var-freq 0.0001 --p-value 0.95 >$outdir/$sample"_varscan_SNP"

samtools mpileup -r chr13:28577411-28674713 -f $fd/hg19_chr13.fa $outdir/$sample"_Aligned.sortedByCoord.out.bam" | java -jar  /crex/proj/snic2020-6-4/wenjiang/MSE/revision2/software/VarScan.v2.3.9.jar mpileup2indel --min-coverage 1  --min-avg-qual 15 --min-var-freq 0.0001 --p-value 0.95 >$outdir/$sample"_varscan_INDEL"

done


## Post process

wkdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/One_SNV_simul_WT_ref_XAEM/ASEP

cd $wkdir
R
workdir='/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/One_SNV_simul_WT_ref_XAEM/ASEP'
# flist1 = list.files(workdir,pattern="_varscan_SNP",recursive=TRUE,full.names = TRUE)
# flist2 = list.files(workdir,pattern="_varscan_INDEL",recursive=TRUE,full.names = TRUE)
options(stringsAsFactors = F)

load('Ycount.RData')
keep_sample=samplename1
wt.ratio=rep(1,length(keep_sample))
names(wt.ratio)=keep_sample
for(i in 1:length(keep_sample))
{
cat(i,'\n')
sname=keep_sample[i]
snp=read.table(paste0('/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/One_SNV_simul_WT_ref_XAEM/ASEP/',sname,'/',sname,'_varscan_SNP'),header=TRUE,sep=':')[,2:3]

snp.wt.ratio=indel.wt.ratio=NULL

if(nrow(snp)>0)snp.wt.ratio=mean(snp[,2]/snp[,1])

indel=read.table(paste0('/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/One_SNV_simul_WT_ref_XAEM/ASEP/',sname,'/',sname,'_varscan_INDEL'),header=TRUE,sep=':')[,2:3]

if(nrow(indel)>0)indel.wt.ratio=mean(indel[,2]/indel[,1])

wt.ratio[sname]=mean(c(snp.wt.ratio,indel.wt.ratio))
}

wt.ratio[is.na(wt.ratio)]=1

save(wt.ratio,file='keep_100_sample_WT_ratio.RData')


## combine XAEM results
load('keep_100_sample_WT_ratio.RData')
load('Beta_final_paralog.Rdata')
load('Ycount.RData')

est=beta.all[[1]]
rownames(est)=samplename1
est=as.data.frame(est)
est1=est
iso.mut=NULL
wt.ratio=wt.ratio[samplename1]

for(i in 1:length(colnames(est1)))
{
isoform=colnames(est1)[i]
iso.mut=c(iso.mut,paste0(isoform,'_mut_all'))
est=cbind(est,iso.mut=est1[,i]*(1-wt.ratio))
est[,i]=est1[,i]*wt.ratio

}

colnames(est)=c(colnames(est1),iso.mut)
save(est,file='Est_ASEP_Varscan.RData')

## pindel 

module load bioinfo-tools 
module load samtools/1.14
module load bwa/0.7.17
module load python/3.7.2
  

ref=/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/reference/hg19_chr13.fa
wkdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision2/FLT3/Pindel_BWA_alignment
fqdir=/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/fasta_file
cd $wkdir

for fn in $(find  $fqdir -type f -name "*_1.fasta")
do
echo "$fn"
outdir=$(echo $output$(basename $fn))
sample=$outdir
endpoint=$(expr $(echo "${#fn}") - 8)
fn2=$(echo $fn|cut -c 1-$endpoint)
fn2=$(echo $fn2"_2.fasta")
echo $outdir
mkdir -p $wkdir/$outdir
cd $wkdir/$outdir

read1=$fn
read2=$fn2
samplename=$sample

# BWA-MEM with -M
bwa mem /proj/snic2020-6-4/Nghia/referenceDB/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa $read1 $read2 -t 8 | samtools sort > $wkdir/$outdir/$samplename".bam"
# Index BAM
samtools index $wkdir/$outdir/$samplename".bam"

#make a config file
echo $wkdir/$outdir/$samplename".bam" 250 $samplename > $wkdir/$outdir/$samplename".config"

##run ScanITD
#/proj/snic2020-6-4/nobackup/Nghia/Working/ASE_isoform/revision2/Pindel/ScanITD/ScanITD.py -i $samplename".bam" -r /proj/snic2020-6-4/Nghia/referenceDB/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa -o $samplename"_scanITD"

ref=/proj/snic2020-6-4/Nghia/referenceDB/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
#run pindel
/proj/snic2020-6-4/nobackup/Nghia/Working/ASE_isoform/revision2/Pindel/pindel-master/pindel -f $ref -c chr13 -o $samplename -i $wkdir/$outdir/$samplename".config" --report_breakpoints FALSE --window_size 10 -T 8
#get vcf for TD only
/proj/snic2020-6-4/nobackup/Nghia/Working/ASE_isoform/revision2/Pindel/pindel-master/pindel2vcf -p $wkdir/$outdir/$samplename"_TD" -r $ref -R GRCh19 -d 20091123 -v $wkdir/$outdir/$samplename.vcf -mc 3 -he 0.0 -c chr13

cd $wkdir
done

##get vcf for all SV types
#/proj/snic2020-6-4/nobackup/Nghia/Working/ASE_isoform/revision2/Pindel/pindel-master/pindel2vcf -P $samplename -r $ref -R GRCh19 -d 20091123 -v $samplename.vcf -mc 1 -he 0.0

