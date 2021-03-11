# MAX: quantification of mutant-allele expression in cancer from RNA sequencing data

## 1. What is MAX?
MAX is a novel method to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. 

## 2. The input data of MAX
MAX requires several files as input, such as:
- A list of mutations containing the information of chromosome, start position, end position, reference sequence, alternative sequence and gene names. Here is an example of the mutation file and the header: 

![image](https://user-images.githubusercontent.com/40486459/110524071-36484600-8113-11eb-9d86-6369007b391c.png)


#### (Please make sure if the mutations are detected using hg19 or hg38 assembly)

- The GTF annotation file, which can be dowlnloaded from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables). The hg19 or hg38 RefGene GTF can be downloaded by running:
```sh
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz # hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz #hg38
```
- The whole transcriptome reference. The clean version of hg19 or hg38 reference can be downloaded by running:
```sh
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz # hg19
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz #hg38
```
- The RNA-seq data from a group of samples.
#### Software requirements for XAEM:
- R version 3.6.0 or later with installed packages: GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38 (or BSgenome.Hsapiens.UCSC.hg19), foreach and doParallel
- C++11 compliant compiler (g++ >= 4.7)
## 2. Download and installation
#### If you use the binary verion of XAEM (recommended):

- Download the latest binary version from XAEM website:
```sh
wget https://github.com/WenjiangDeng/XAEM/raw/master/XAEM-binary-0.1.0.tar.gz
```
- Uncompress to folder
```sh
tar -xzvf XAEM-binary-0.1.0.tar.gz
```
- Move to the XAEM_home directory and do configuration for XAEM
```sh
cd XAEM-binary-0.1.0
bash configure.sh
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/XAEM-binary-0.1.0/lib:$LD_LIBRARY_PATH
export PATH=/path/to/XAEM-binary-0.1.0/bin:$PATH
```
#### If you want to build XAEM from sources:

- Download XAEM from XAEM website and move to XAEM_home directory
```sh
wget https://github.com/WenjiangDeng/XAEM/raw/master/XAEM-source-0.1.0.zip
unzip XAEM-source-0.1.0.zip
cd XAEM-source-0.1.0
bash configure.sh
```
XAEM requires information of flags from Sailfish including DFETCH_BOOST, DBOOST_ROOT, DTBB_INSTALL_DIR and DCMAKE_INSTALL_PREFIX. Please refer to the Sailfish website for more details of these flags.
- Do installation by the following command:
```sh
DBOOST_ROOT=/path/to/boostDir/ DTBB_INSTALL_DIR=/path/to/tbbDir/ DCMAKE_INSTALL_PREFIX=/path/to/expectedBuildDir bash install.sh
After the installation is finished, remember to add the paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
export LD_LIBRARY_PATH=/path/to/expectedBuildDir/lib:$LD_LIBRARY_PATH
export PATH=/path/to/expectedBuildDir/bin:$PATH
```
#### Do not forget to replace "/path/to/" by your local path.
## 3. Generate the wild-type + mutant reference, reference index and the design matrix X

This step will construct (1) the reference which contains both wild-type and mutant alleles; (2) the index for the reference and (3) the X matrix (design matrix). This step requires the following input files: a list of mutations, the GTF file, the wild-type transcriptome reference, the version of gene model ("hg19" or "hg38") and the working directory. When you have prepared these files, the command to start the analysis is:

```sh
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
#gunzip hg38.refGene.gtf.gz 
#wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
#gunzip refMrna.fa.gz

bash pipeline.sh -m mutation_list.txt -g hg38.refGene.gtf -r refMrna.fa -v hg38 -d /path/to/directory

```
The outputs will be **WT_Mut_tx_ref.final.fa**, **X_matrix.RData** and the **Index folder**.
## 4. Quantifcation of mutant-allele expression
Suppose we already created a working directory “MAX_project” (/path/to/MAX_project/) for the quantification.
### 4.1 Generate the equivalence class table and Y count matrix
- The command to generate equivalence class table for each sample is similar to [“salmon quant”](https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon). For example, we want to run MAX for sample1 and sample2 with 8 cpus:
```sh
MAX -i /path/to/Index -l IU -1 s1_read1.fasta -2 s1_read2.fasta -p 8 -o /path/to/MAX_project/sample1
MAX -i /path/to/Index -l IU -1 s2_read1.fasta -2 s2_read2.fasta -p 8 -o /path/to/MAX_project/sample2
```
- If the data is compressed in gz format, we can combine "gunzip" in the command:
```sh
XAEM -i /path/to/Index -l IU -1 <(gunzip -c s1_read1.gz) -2 <(gunzip -c s1_read2.gz) -p 8 -o /path/to/MAX_project/sample1
XAEM -i /path/to/Index -l IU -1 <(gunzip -c s2_read1.gz) -2 <(gunzip -c s2_read2.gz) -p 8 -o /path/to/MAX_project/sample2
```
- After running MAX there will be the output of the equivalence class table for multiple samples. We then create the Y count matrix. For example, if we want to run XAEM parallelly using 8 cores, the command is:

```sh
Rscript Create_count_matrix.R workdir=/path/to/XAEM_project core=8
```
### 4.2 Estimate the transcript expression using AEM algorithm
When the Y count matrix is constructed, we can use the AEM algorithm to quantify the mutant-allele expression. The command is as follows:

```sh
Rscript AEM_update_X_beta.R workdir=/path/to/XAEM_project core=8 design.matrix=/path/to/X_matrix.RData isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData merge.paralogs=FALSE isoform.method=average remove.ycount=TRUE
```
#### Parameter setting
- **workdir**: the path to working directory
- **core**: the number of cpu cores for parallel computing
- **design.matrix**: the path to the design matrix
- **remove.ycount** (default=TRUE): to clean all data of Ycount after use.
- xxxx
## 5. A complete run of MAX by copy and paste
This section shows the tutorial to run MAX pipeline. We can test MAX by just copy and paste of the example commands.

- Download the binary file of MAX
```sh
mkdir tmp_test
cd tmp_test
wget https://github.com/WenjiangDeng/XAEM/raw/master/XAEM-binary-0.1.0.tar.gz
tar -xzvf XAEM-binary-0.1.0.tar.gz
cd XAEM-binary-0.1.0
bash configure.sh
export LD_LIBRARY_PATH=/path/to/XAEM-binary-0.1.0/lib:$LD_LIBRARY_PATH
export PATH=/path/to/XAEM-binary-0.1.0/bin:$PATH
```
- Download fasta file and index it
```sh
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/transcripts.fa.gz
gunzip transcripts.fa.gz
TxIndexer -t transcripts.fa -o TxIndexer_idx
```
- Download the X matrix and RNA-seq data of sample1 and sample2
```sh
mkdir XAEM_project
cd XAEM_project
wget https://github.com/WenjiangDeng/XAEM/raw/master/X_matrix.RData
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample1_read1.fasta.gz
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample1_read2.fasta.gz
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample2_read1.fasta.gz
wget http://fafner.meb.ki.se/biostatwiki/2018_XAEM/sample2_read2.fasta.gz
cd ..
```
- Generate the eqclass table and Y count matrix
```sh
MAX -i TxIndexer_idx -l IU -1 <(gunzip -c XAEM_project/sample1_read1.fasta.gz) -2 <(gunzip -c XAEM_project/sample1_read2.fasta.gz) -p 4 -o XAEM_project/eqc_sample1
MAX -i TxIndexer_idx -l IU -1 <(gunzip -c XAEM_project/sample2_read1.fasta.gz) -2 <(gunzip -c XAEM_project/sample2_read2.fasta.gz) -p 4 -o XAEM_project/eqc_sample2
## R packages foreach and doParallel are required

Rscript Create_count_matrix.R workdir=$PWD/XAEM_project core=8
```
- Estimate isoform expression using AEM algorithm
```sh
Rscript AEM_update_X_beta.R workdir=$PWD/XAEM_project core=8
cd XAEM_project
```
Reference: tba
