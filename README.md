# MAX: quantification of mutant-allele expression in cancer from RNA sequencing data

## 1. What is MAX?
MAX is a novel method to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. 

## 2. The input data of MAX
MAX requires several files as input, such as:
- A list of mutations containing the information of chromosome, start position, end position, reference sequence and alternative sequence. Here is an example of the mutation file:

  ![image](https://user-images.githubusercontent.com/40486459/110201429-5520b100-7e63-11eb-9efd-e57f12793b66.png)
- The reference sequence of whole genome. In this study, we use the hg38 gene model downloaded [from UCSC repository](https://hgdownload.soe.ucsc.edu/downloads.html#human/).
- The GTF annotation file of the specific genes of interest. The hg38 GTF file can be retrieved from [UCSC Table Browser](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz). 
- The RNA-seq data from a group of samples.
#### Software requirements for XAEM:
- R version 3.3.0 or later with installed packages: foreach and doParallel
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
Do not forget to replace "/path/to/" by your local path.
## 5. A complete run of MAX by copy and paste
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
