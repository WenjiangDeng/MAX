### 2022-02-01 WJ Deng
### Simulate 100 mutated RNA-seq samples, in which both wt and mut allele are expressed
### Oncogenes in AML: FLT3, NPM1, TP53
### The simulation mimics the real mutation landscape from BeatAML patients

######
###### FLT3
######

load('FLT3_Polyester.RData') 

t1=revision1_simulation_matrix
full.202.mut.events=read.csv('FLT3-Mutation-Profile-BeatAML.csv')
set.seed(2021)
lab.id=sample(full.202.mut.events$labId,100)

for(i in 1:100)
{
lab.id1=lab.id[i]
tx_NM_004119=full.202.mut.events[full.202.mut.events$labId==lab.id1,'ID']*2-1

tx_NR_130706=full.202.mut.events[full.202.mut.events$labId==lab.id1,'ID']*2

t1[tx_NM_004119,i]=t1[241,i]
t1[tx_NR_130706,i]=t1[242,i]
}

true=t1
save(true,file='matrix.polyester_FLT3.RData')
save(true,file='True_read_counts.RData')

### Simulate RNA-seq reads using R package Polyester

cd /crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/fasta_file

module load bioinfo-tools
module load R_packages/3.6.0
R
library(polyester)
load('matrix.polyester_FLT3.RData')
fasta_file='FLT3_mut_wt_final.fa'
outdir='/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/flt3_ITD_simul_real_beatAML_correct/fasta_file'

simulate_experiment_countmat(fasta_file, readmat=true, outdir=outdir,error_rate=0, strand_specific=FALSE) 

######
###### NPM1
######

load('NPM1_Polyester.RData')
mut.profile=read.csv('NPM1-Mutation-Profile-BeatAML.csv')
set.seed(2021)
lab.id=sample(mut.profile$labId,100)

for(i in 1:100)
{	
lab.id1=lab.id[i]
tx_NM_001355006=mut.profile[mut.profile$labId==lab.id1,'id']*6-5
tx_NM_001355007=mut.profile[mut.profile$labId==lab.id1,'id']*6-4
tx_NM_001355010=mut.profile[mut.profile$labId==lab.id1,'id']*6-3
tx_NM_002520=mut.profile[mut.profile$labId==lab.id1,'id']*6-2
tx_NM_199185=mut.profile[mut.profile$labId==lab.id1,'id']*6-1
tx_NR_149149=mut.profile[mut.profile$labId==lab.id1,'id']*6

t3[tx_NM_001355006,i]=true.wt[2,i]
t3[tx_NM_001355007,i]=true.wt[3,i]
t3[tx_NM_001355010,i]=true.wt[5,i]
t3[tx_NM_002520,i]=true.wt[6,i]
t3[tx_NM_199185,i]=true.wt[7,i]
t3[tx_NR_149149,i]=true.wt[8,i]
}
true = t3

save(true,file='matrix.polyester_NPM1.RData')
save(true,file='True_read_counts.RData')

cd /crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/fasta_file

module load bioinfo-tools
module load R_packages/3.5.0
R
library(polyester)
load('matrix.polyester_NPM1.RData')
fasta_file='NPM1_mut_wt_final.fa'
outdir='/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/NPM1_simul_beatAML/fasta_file'
simulate_experiment_countmat(fasta_file, readmat=t3, outdir=outdir,error_rate=0, strand_specific=FALSE) 


####
#### TP53
####
load('TP53_Polyester.RData')
mut=tp53_isoform_list
true=tp53_wt_exp
t1=matrix(0,582,100)
colnames(t1)=colnames(true);rownames(t1)=mut[,1]
t1[1:15,]=true
t2=t1[-c(1:15),]
index=sapply(rownames(t2),function(x)strsplit(x,'_')[[1]][4])
index=as.numeric(index)
set.seed(2021)
mut.profile.tp53=read.csv('TP53-Mutation-Profile-BeatAML.csv')
sample.id=sample(mut.profile.tp53$mut_id,100,replace=T)
colnames(t2)=colnames(true)

for(id in 1:100)
{	
	mut.id=sample.id[id] # mut id
	cat(id,'  ',mut.id,'\n')
	t2.1=NULL
	t2.1=t2[index>(mut.id-1)*30  &   index<mut.id*30,]
	t2.1.1=sapply(rownames(t2.1),function(x)strsplit(x,'_mut')[[1]][1])
	t2[index>(mut.id-1)*30  &   index<mut.id*30,id]=true[t2.1.1,id]
	
	}
	
true=rbind(true,t2)

save(true,file='matrix.polyester_TP53.RData')
save(true,file='True_read_counts.RData')

module load bioinfo-tools
module load R_packages/3.5.0
R
library(polyester)
load('matrix.polyester_TP53.RData')
fasta_file='TP53_mut_wt_final.fa'
outdir='/crex/proj/snic2020-6-4/wenjiang/MSE/revision1/TP53_simul_beatAML/fasta'
simulate_experiment_countmat(fasta_file, readmat=true, outdir=outdir,error_rate=0, strand_specific=FALSE) 


