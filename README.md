# MAX: quantification of mutant-allele expression in cancer from RNA sequencing data

## 1. What is MAX?
MAX is a novel method to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. 

## 2. The input data of MAX
MAX requires several files as input, such as:
- A list of mutations containing the information of chromosome, start position, end position, reference sequence, alternative sequence and gene names. Here is an example of [the mutation file with the header](https://github.com/WenjiangDeng/MAX/blob/main/mutation_list.txt): 

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
## 3. Construct the wild-type + mutant reference, reference index and the design matrix X

This step will produce (1) the reference which contains both wild-type and mutant alleles; (2) the index for the reference and (3) the X matrix (design matrix). This step requires the following input files: a list of mutations, the GTF file, the wild-type transcriptome reference, the version of gene model ("hg19" or "hg38") and the working directory. When you have prepared these files, the command to start the analysis is:

```sh
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz 
# wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz

bash pipeline.sh -m mutation_list.txt -g hg38.refGene.gtf -r refMrna.fa -v hg38 -d /path/to/directory

```
The outputs will be: **WT_Mut_reference.fa**, **X_matrix.RData** and the **Index_reference** folder.
## 4. Quantifcation of mutant-allele expression
Suppose we already created a working directory “MAX_project” (/path/to/MAX_project/) for the quantification.
### 4.1 Generate the equivalence class table and Y count matrix
- The command to generate equivalence class table for each sample is similar with [“salmon quant”](https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon). The input parameters are the Index_reference folder and RNA-seq data. The "-l" option defines the [library type of reads](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype). Below shows the example if we want to run MAX for sample1 and sample2 with 8 cpus:
```sh
MAX -i /path/to/Index_reference -l IU -1 s1_read1.fasta -2 s1_read2.fasta -p 8 -o /path/to/MAX_project/sample1 -w 100000000
MAX -i /path/to/Index_reference -l IU -1 s2_read1.fasta -2 s2_read2.fasta -p 8 -o /path/to/MAX_project/sample2 -w 100000000
```
- If the data is compressed in gz format, we can combine "gunzip" in the command:
```sh
MAX -i /path/to/Index_reference -l IU -1 <(gunzip -c s1_read1.gz) -2 <(gunzip -c s1_read2.gz) -p 8 -o /path/to/MAX_project/sample1 -w 100000000
MAX -i /path/to/Index_reference -l IU -1 <(gunzip -c s2_read1.gz) -2 <(gunzip -c s2_read2.gz) -p 8 -o /path/to/MAX_project/sample2 -w 100000000
```
- After running MAX there will be the output of the equivalence class table for multiple samples. We then create the Y count matrix. For example, if we want to run XAEM parallelly using 8 cores, the command is:

```sh
Rscript Create_count_matrix.R design.matrix=/path/to/X_matrix.RData workdir=/path/to/MAX_project core=8
```
### 4.2 Estimate the transcript expression using AEM algorithm
When the Y count matrix is constructed, we can use the AEM algorithm to quantify the mutant-allele expression. The command is as follows:

```sh
Rscript AEM_update_X_beta.R workdir=/path/to/XAEM_project design.matrix=/path/to/X_matrix.RData max.out=/path/to/mutant_expression.RData remove.ycount=TRUE core=8
```
#### Parameter setting
- **workdir**: the path to working directory
- **design.matrix**: the path to the design matrix
- **max.out**: the output file for isoform expression
- **remove.ycount** (default=TRUE): to clean all data of Ycount after use
- **core**: the number of cpu cores for parallel computing, default is 8.

The final results are in the mutant_expression.RData, which contains two objects: the **MAX_count** for the read counts value and **MAX_tpm** for the TPM (Transcripts Per Kilobase Million) value.
## 5. An example run of MAX by copy and paste
This section shows a complete run for MAX pipeline. We can test MAX just by copy and paste of the commands. Here we focus on the mutations in the FLT3 gene, which is one of the most frequently mutated oncogenes in acute myeloid leukemia (AML). 

- Download the binary file of MAX and configure the path
```sh
# Create a working folder
mkdir MAX_binary
cd MAX_binary
# Download the binary version of XAEM
wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/XAEM_datasources/XAEM-binary-0.1.1.tar.gz

# Configure the tool
tar -xzvf XAEM-binary-0.1.1.tar.gz
cd XAEM-binary-0.1.1
bash configure.sh

# Add the paths to system
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
export PATH=$PWD/bin:$PATH
cd ..

```
- Download the mutation file, GTF annotation and the sequences of isoforms from FLT3 gene
```sh
mkdir XAEM_project
cd XAEM_project

wget https://github.com/WenjiangDeng/MAX/blob/main/mutation_list.txt
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
gunzip hg19.refGene.gtf.gz
wget https://github.com/WenjiangDeng/MAX/blob/main/isoform_ref_FLT3_gene.fa

```
- construct the wild-type+Mutant reference, reference index and the X matrix
```sh

bash pipeline.sh -m mutation_list.txt -g hg19.refGene.gtf -r isoform_ref_FLT3_gene.fa -v hg19 -d $PWD

```
- Download the test RNA-seq data from 10 samples
```sh

wget https://github.com/WenjiangDeng/MAX/blob/main/fasta_FLT3.tar.gz
tar -zxvf fasta_FLT3.tar.gz
ll -tr fasta_flt3/

```
- Generate the eqclass table and Y count matrix
```sh
for ind in $(seq -f %02.0f  10); do
MAX -i Index_reference -l IU -1 'fasta_flt3/sample_'$ind'_1.fasta' -2 'fasta_flt3/sample_'$ind'_2.fasta' -p 8 -o 'sample_'$ind -w 100000000
done

Rscript ../MAX_binary/XAEM-binary-0.1.1/R/Create_count_matrix.R workdir=$PWD design.matrix=X_matrix.RData core=8

```
- Estimate mutant-allele expression using AEM algorithm
```sh
Rscript ../MAX_binary/XAEM-binary-0.1.1/R/AEM_update_X_beta.R workdir=$PWD design.matrix=X_matrix.RData max.out=mutant_expression.RData core=8
```

#### Reference: tba
