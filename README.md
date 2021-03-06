# MAX: quantification of mutant-allele expression at isoform level

## 1. What is MAX?
MAX is a novel method to quantify the mutant-allele expression at isoform level from RNA-seq data. 

## 2. The input data of MAX.
MAX requires several files as input, such as:
- A list of mutations containing the information of chromosome, start position, end position, reference sequence and alternative sequence. Here is an example of the mutation file:

  ![image](https://user-images.githubusercontent.com/40486459/110201429-5520b100-7e63-11eb-9efd-e57f12793b66.png)
- The reference of genome or the specific chromosome. In this study, we use the [hg38 gene model downloaded from UCSC repository](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/).
- The reference of whole transcriptome or the isoforms from a specific gene. The transcritpome sequence can be retrieved from [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables). 
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

Download XAEM from XAEM website and move to XAEM_home directory
```sh
wget https://github.com/WenjiangDeng/XAEM/raw/master/XAEM-source-0.1.0.zip
unzip XAEM-source-0.1.0.zip
cd XAEM-source-0.1.0
bash configure.sh
```
XAEM requires information of flags from Sailfish including DFETCH_BOOST, DBOOST_ROOT, DTBB_INSTALL_DIR and DCMAKE_INSTALL_PREFIX. Please refer to the Sailfish website for more details of these flags.
Do installation by the following command:
```sh
DBOOST_ROOT=/path/to/boostDir/ DTBB_INSTALL_DIR=/path/to/tbbDir/ DCMAKE_INSTALL_PREFIX=/path/to/expectedBuildDir bash install.sh
After the installation is finished, remember to add the paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
export LD_LIBRARY_PATH=/path/to/expectedBuildDir/lib:$LD_LIBRARY_PATH
export PATH=/path/to/expectedBuildDir/bin:$PATH
```
Do not forget to replace "/path/to/" by your local path.
