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
