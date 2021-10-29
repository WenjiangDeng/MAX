# MAX: quantification of mutant-allele expression in cancer from RNA sequencing data

## What is MAX?
MAX is a novel method to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. 
## 1. The input data of MAX
MAX requires several files as input, such as:
- A list of mutations containing the information of chromosome, start position, end position, reference sequence, alternative sequence and gene names. Here is an example of [the mutation file with the header](https://github.com/WenjiangDeng/MAX/raw/main/mutation_list.txt): 

![image](https://user-images.githubusercontent.com/40486459/110524071-36484600-8113-11eb-9d86-6369007b391c.png)


#### (Please make sure if the mutations are detected using hg19 or hg38 assembly)

- The GTF annotation file, which can be dowlnloaded from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables). The hg19 or hg38 RefGene GTF can be downloaded by running:
```sh
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz # hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz #hg38
```
- The whole transcriptome reference, or the isoforms of a subset of genes that you aim to focus on. The clean version of hg19 or hg38 reference can be downloaded by running:
```sh
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz # hg19
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz #hg38
```
- The RNA-seq data from multiple samples.
#### Software requirements for MAX:
- R version 3.6.0 or later with installed packages: GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38 (or BSgenome.Hsapiens.UCSC.hg19), polyester, Biostrings, foreach and doParallel
- C++11 compliant compiler (g++ >= 4.7)
## 2. Download and installation

#### If you use the binary verion of MAX (recommended):

- Download the latest binary version from the MAX release:
```sh
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.1.0/MAX-binary-0.1.0.tar.gz

```
- Uncompress to folder
```sh
tar -xzvf MAX-binary-0.1.0.tar.gz
```
- Move to the MAX_home directory and do configuration for MAX
```sh
cd MAX-binary-0.1.0
bash configure.sh
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/XAEM-binary-0.1.0/lib:$LD_LIBRARY_PATH
export PATH=/path/to/XAEM-binary-0.1.0/bin:$PATH
#Done
```
#### If you want to build MAX from sources:

- download MAX
```sh
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.1.0/MAX-source-0.1.0.tar.gz
tar -xzvf MAX-source-0.1.0.tar.gz
cd MAX-source-0.1.0

#config to run MAX
bash configure.sh
```
- install boost_1_55_0
```sh
wget http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.gz
tar -xvzf boost_1_58_0.tar.gz
cd boost_1_58_0

sudo apt-get update
sudo apt-get install build-essential g++ python-dev autotools-dev libicu-dev build-essential libbz2-dev libboost-all-dev
sudo apt-get install aptitude
aptitude search boost

./bootstrap.sh --prefix=boost_1_58_0_build
./b2
./b2 install
```
#The Boost C++ Libraries were successfully built!
#add the lib and folder to paths
```sh
export LD_LIBRARY_PATH=$PWD/boost_1_58_0_build/stage/lib:$LD_LIBRARY_PATH
export PATH=$PWD/boost_1_58_0_build:$PATH
```
- install tbb44_20160526oss
```sh
cd ..
wget https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb44_20160526oss_src_0.tgz
tar xvf tbb44_20160526oss_src_0.tgz
sudo apt-get install libtbb-dev

# install cmake for ubuntu: cmake 3.5.1
sudo apt install cmake
# install curl
sudo apt install curl
# install autoconf
sudo apt-get install autoconf
# install zlib
sudo apt install zlib1g-dev
sudo apt install zlib1g
# update all installations
sudo apt-get update

# install MAX
DBOOST_ROOT=$PWD/boost_1_58_0/boost_1_58_0_build/ DTBB_INSTALL_DIR=$PWD/tbb44_20160526oss/ DCMAKE_INSTALL_PREFIX=MAX-source-0.1.0 bash install.sh

#The MAX was successfully built!
###########

#add lib and bin folders to paths
export LD_LIBRARY_PATH=$PWD/Circall_0.1.0_build/lib:$LD_LIBRARY_PATH
export PATH=$PWD/Circall_0.1.0_build/bin:$PATH

#done
```
#### Do not forget to replace "/path/to/" by your local path.
## 3. Construct the wild-type + mutant reference, reference index and the X matrix
#### In Section 5 we show a test run of MAX just by copy and paste.
This step will produce (1) the reference which contains both wild-type and mutant alleles; (2) the index for the reference and (3) the X matrix (design matrix). This step requires the following input files: a list of mutations, the GTF file, the wild-type transcriptome reference, the version of gene model ("hg19" or "hg38") and the working directory. When you have prepared these files, the command to start the analysis is:

```sh
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz 
# wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz

bash MAX.sh -m mutation_list.txt -g hg38.refGene.gtf -r refMrna.fa -v hg38 -d /path/to/directory

```
The outputs will be: **WT_Mut_reference.fa**, **X_matrix.RData** and the **Index_reference** folder.
## 4. Quantifcation of mutant-allele expression
Suppose we already created a working directory “MAX_project” (/path/to/MAX_project/) for the quantification.
### 4.1 Generate the equivalence class table and Y count matrix
- The command to generate equivalence class table for each sample is similar with [“salmon quant”](https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon). The input parameters are the Index_reference folder and the RNA-seq data. The "-l" option defines the [library type of reads](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype). Below shows the example if we want to run MAX for sample1 and sample2 with 8 cpus:
```sh
MAX -i /path/to/Index_reference -l IU -1 s1_read1.fasta -2 s1_read2.fasta -p 8 -o /path/to/MAX_project/sample1 
MAX -i /path/to/Index_reference -l IU -1 s2_read1.fasta -2 s2_read2.fasta -p 8 -o /path/to/MAX_project/sample2 
```
- If the data is compressed in gz format, we can combine "gunzip" in the command:
```sh
MAX -i /path/to/Index_reference -l IU -1 <(gunzip -c s1_read1.gz) -2 <(gunzip -c s1_read2.gz) -p 8 -o /path/to/MAX_project/sample1 
MAX -i /path/to/Index_reference -l IU -1 <(gunzip -c s2_read1.gz) -2 <(gunzip -c s2_read2.gz) -p 8 -o /path/to/MAX_project/sample2 
```
- After running MAX there will be the output of the equivalence class table for multiple samples. We then create the Y count matrix. For example, if we want to run MAX parallelly using 8 cores, the command is:

```sh
Rscript Create_count_matrix.R design.matrix=/path/to/X_matrix.RData workdir=/path/to/MAX_project core=8
```
### 4.2 Estimate the transcript expression using AEM algorithm
When the Y count matrix is constructed, we can use the AEM algorithm to quantify the mutant-allele expression. The command is as follows:

```sh
Rscript AEM_update_X_beta.R workdir=/path/to/MAX_project design.matrix=/path/to/X_matrix.RData max.out=/path/to/mutant_expression.RData remove.ycount=TRUE core=8
```
#### Parameter setting
- **workdir**: the path to working directory
- **design.matrix**: the path to the design matrix
- **max.out**: the output file for isoform expression
- **remove.ycount** (default=TRUE): to clean all data of Ycount after use
- **core**: the number of cpu cores for parallel computing, default is 8.

The final results are in the mutant_expression.RData, which contains two objects: the **MAX_count** for the read counts value and **MAX_tpm** for the TPM (Transcripts Per Kilobase Million) value.
## 5. A trial run of MAX by copy and paste
This section shows a complete run for MAX pipeline. We can test MAX just by copy and paste of the commands. Here we focus on the mutations in the FLT3 gene, which is one of the most frequently mutated oncogenes in Acute Myeloid Leukemia (AML). 

- Download the binary file of MAX and configure the path
```sh
# Create a working folder
mkdir MAX_binary
cd MAX_binary
# Download the binary version of MAX
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.1.0/MAX-binary-0.1.0.tar.gz

# Configure the tool
tar -xzvf MAX-binary-0.1.0.tar.gz
cd MAX-binary-0.1.0
bash configure.sh
# Add the paths to system
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
export PATH=$PWD/bin:$PATH
cd ..



```
- Download the mutation file, GTF annotation and the sequences of isoforms from FLT3 gene
```sh
mkdir MAX_project
cd MAX_project

wget https://github.com/WenjiangDeng/MAX/releases/download/v0.1.0/mutation_list.txt
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.1.0/test_FLT3.gtf
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.1.0/isoform_ref_FLT3_gene.fa

```
- construct the wild-type+Mutant reference, reference index and the X matrix
```sh
# This step can take several minutes

bash ../MAX-binary-0.1.0/MAX.sh -m mutation_list.txt -g test_FLT3.gtf -r isoform_ref_FLT3_gene.fa -v hg19 -d $PWD



```
- Download the test RNA-seq data of 10 samples
```sh

wget https://github.com/WenjiangDeng/MAX/releases/download/v0.1.0/RNA-seq_FLT3.tar.gz
tar -xzvf RNA-seq_FLT3.tar.gz


```
- Generate the eqclass table and Y count matrix using 8 cores
```sh
for ind in $(seq -f %02.0f  10); do
MAX -i Index_reference -l IU -1 'fasta_flt3/sample_'$ind'_1.fasta' -2 'fasta_flt3/sample_'$ind'_2.fasta' -p 8 -o 'sample_'$ind 
done

Rscript ../MAX-binary-0.1.0/R/Create_count_matrix.R workdir=$PWD design.matrix=X_matrix.RData core=8



```
- Estimate mutant-allele expression using AEM algorithm with 8 cores
```sh
Rscript ../MAX-binary-0.1.0/R/AEM_update_X_beta.R workdir=$PWD design.matrix=X_matrix.RData max.out=mutant_expression.RData core=8


```
The final results are in the **mutant_expression.RData**, which contains the MAX_count and MAX_tpm objects. 
## MAX2
MAX2 is an extension of MAX for heterogenous RNA-seq data. In MAX2, the construction of X matrix and the quasi-mapping step are the same as in MAX. The only difference is that before isoform quantification, MAX2 will cluster heterogenous samples based on their mutation profile.

- Merge mutation to generate Mutated_eqClass.txt. 
```sh
# SampleMut has two columns, the first column is Sample name, the second column is the unique identifier of a mutation .SampleEq is the eqClass.txt file from previous quasi-mapping step.
Rscript mergeMutSingleSample.R sampleMut=$sampleMutFn sampleID=$sampleFn sampleEq=eqClass.txt
```
- Generate the ycount and quantify using the AEM algorithm
```sh
Rscript genCountSample.R xmatEq=$xmatEqFn sampleMut=$sampleMutFn sampleID=$sampleFn sampleEq=$outdir/Mutated_eqClass.txt YcountDir=$Ycount_outdir

Rscript estimateBeta.R workdir=$Ycount_outdir sampleMut=$sampleMutFn out=MAX_isoform_expression_AEM.RData

```
The estimation results will be saved as isoformCount in **MAX2_isoform_expression.RData**.
#### Reference: tba
