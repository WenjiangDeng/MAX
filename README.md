# MAX: quantification of mutant-allele expression in cancer from RNA sequencing data

## What is MAX?
MAX is a novel method to quantify the Mutant-Allele eXpression (MAX) at isoform level from RNA-seq data. 
## 1. The input data of MAX
MAX requires several files as input, such as:
- A list of mutations containing the information of chromosome, start position, end position, reference sequence, alternative sequence and gene names. Here is an example of [the mutation file with the header](https://github.com/WenjiangDeng/MAX/blob/main/testData/test_mutation_list.txt): 

![image](https://user-images.githubusercontent.com/40486459/110524071-36484600-8113-11eb-9d86-6369007b391c.png)

**If using MAX2**, an extra column "Samples" consisting of IDs of samples carrying the mutations is required. The ID of a sample must be in the name of the fastq sequence of that sample (must be in gzip format). For example patient01_1.fastq.gz and patient01_2.fastq.gz are the sequencing data files of sample "patient01". Please see [an example of the mutation file for MAX2 here](https://github.com/WenjiangDeng/MAX/blob/main/testData/test_mutation_list_MAX2.txt): 

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
- R version 3.6.0 or later with installed packages: GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38 (or BSgenome.Hsapiens.UCSC.hg19), polyester, Biostrings, data.table, foreach and doParallel
- C++11 compliant compiler (g++ >= 4.7)
## 2. Download and installation

#### If you use the binary verion of MAX (recommended):

- Download the latest binary version from the MAX release:
```sh
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.2.0/MAX-binary-0.2.0.tar.gz

```
- Uncompress to folder
```sh
tar -xzvf MAX-binary-0.2.0.tar.gz
```
- Move to the MAX_home directory and do configuration for MAX
```sh
cd MAX-binary-0.2.0
bash configure.sh
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/MAX-binary-0.2.0/lib:$LD_LIBRARY_PATH
export PATH=/path/to/MAX-binary-0.2.0/bin:$PATH
#Done
```
#### If you want to build MAX from sources:

- download MAX
```sh
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.2.0/MAX-source-0.2.0.tar.gz
tar -xzvf MAX-source-0.2.0.tar.gz
cd MAX-source-0.2.0

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
```

- finally, install MAX
```sh
DBOOST_ROOT=$PWD/boost_1_58_0/boost_1_58_0_build/ DTBB_INSTALL_DIR=$PWD/tbb44_20160526oss/ DCMAKE_INSTALL_PREFIX=MAX_0.2.0_build bash install.sh

#The MAX was successfully built!
###########

#add lib and bin folders to paths
export LD_LIBRARY_PATH=$PWD/MAX_0.2.0_build/lib:$LD_LIBRARY_PATH
export PATH=$PWD/MAX_0.2.0_build/bin:$PATH

#done
```
#### Do not forget to replace "/path/to/" by your local path.


## 3. Parameter settings

The examples of parameter setting files can be downloaded here: [1) param file for running MAX](https://github.com/WenjiangDeng/MAX/blob/main/testData/test_params.sh) and [2) param file for running MAX2](https://github.com/WenjiangDeng/MAX/blob/main/testData/test_paramsMAX2.sh).
The param file requires the following parameter:

- **mutlist**: a file contains a list of mutations containing the information of chromosome, start position, end position, reference sequence, alternative sequence, and gene names. If you want to run MAX2, an extra column "Sample" is required. Examples of mutation lists are provided here: [example file for MAX](https://github.com/WenjiangDeng/MAX/blob/main/testData/test_mutation_list.txt) and [example file for MAX2](https://github.com/WenjiangDeng/MAX/blob/main/testData/test_mutation_list_MAX2.txt).
- **gtffile**: gtf of the gene carring mutations
- **fastafile**: sequences of the transcripts of the gene
- **hgversion**: annotation version of the human genome: hg19 or hg38
- **CPUNUM**: the number of threads
- **INPUT**: the path to the folder of the input fastq.gz files
- **OUTPUT**: the path to the output folder
- **useMAX2**: use MAX2 (1) or not (0)

## 4. Run MAX
Given a params.sh with the proper parameter setting for your data, there are two ways to run MAX: 1) using the installed version and 2) using docker.
- To run with the installed version:
```sh
# Download the binary version of MAX
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.2.0/MAX-binary-0.2.0.tar.gz

# Configure the tool
tar -xzvf MAX-binary-0.2.0.tar.gz
cd MAX-binary-0.2.0
bash configure.sh
# Add the paths to system
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
export PATH=$PWD/bin:$PATH
cd ..

# run MAX
runMAX.sh -param params.sh
```
- To run with docker,
```sh
# Pull the docker image of MAX:
docker pull nghiavtr/max:v0.2.0

#download the script to run MAX using docker
wget https://raw.githubusercontent.com/WenjiangDeng/MAX/main/runMAXdocker.sh

# run MAX via the docker:
runMAXdocker.sh -param params.sh
```

The final results are in the MAX_isoform_expression.RData or MAX2_isoform_expression.RData (if using MAX2), which contains two objects: the **isoform_count** for the read counts value and **isoform_tpm** for the TPM (Transcripts Per Kilobase Million) value.

## 5. A trial run of MAX by copy and paste
This section shows a complete run for MAX pipeline. We can test MAX just by copy and paste of the commands. Here we focus on the mutations in the FLT3 gene, which is one of the most frequently mutated oncogenes in Acute Myeloid Leukemia (AML). 


- Download the package containing a dataset of 10 samples, annotation and files for parameter settings to run MAX and MAX2

```sh

wget https://www.meb.ki.se/sites/biostatwiki/wp-content/uploads/sites/4/2022/02/testData.tar.gz
tar -xzvf testData.tar.gz

#move all files and folders of testData folder to the current folder
mv testData/* ./



# Uncompress the example dataset of 10 samples
tar -xzvf test_FLT3_data.tar.gz


ls

#List of the mutation file, GTF annotation and the sequences of isoforms from FLT3 gene are stored in the files below

# annotations for example gene: FLT3
#test_FLT3.gtf
#test_transcripts_FLT3.fa

# mutation list and parameters to run MAX
#test_mutation_list.txt
#test_params.sh

# mutation list and parameters to run MAX2
#test_mutation_list_MAX2.txt
#test_paramsMAX2.sh

```


- To run using installed MAX, need to download the binary file of MAX and configure the path
```sh
# Download the binary version of MAX
wget https://github.com/WenjiangDeng/MAX/releases/download/v0.2.0/MAX-binary-0.2.0.tar.gz

# Configure the tool
tar -xzvf MAX-binary-0.2.0.tar.gz
cd MAX-binary-0.2.0
bash configure.sh
# Add the paths to system
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
export PATH=$PWD/bin:$PATH
cd ..


### This step can take several minutes
# run MAX
bash MAX-binary-0.2.0/runMAX.sh -param test_params.sh

# run MAX2
bash MAX-binary-0.2.0/runMAX.sh -param test_paramsMAX2.sh

```


- To run using docker, we suppose docker was installed in your system
```sh
### This step can take several minutes

# Pull the docker image of MAX:
docker pull nghiavtr/max:v0.2.0

#download the script to run MAX using docker
wget https://raw.githubusercontent.com/WenjiangDeng/MAX/main/runMAXdocker.sh

# run MAX via the docker:
bash runMAXdocker.sh -param test_params.sh


# run MAX2 via the docker:
bash runMAXdocker.sh -param test_paramsMAX2.sh

```


#### Reference: tba
