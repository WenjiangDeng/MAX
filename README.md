# MAX: quantification of mutant-allele expression at isoform level

## 1. What is MAX?
MAX is a novel method to quantify the mutant-allele expression at isoform level from RNA-seq data. 

## 2. The input data of MAX.
MAX requires several files as input. Such as:
- A list of mutations containing the information of chromosome, start position, end position, reference sequence and alternative sequence. Here is an example of the mutation file:
  chromosome	start_pos	end_pos	REF	ALT
13	28592620	28592620	T	C
13	28592622	28592622	G	T
13	28592623	28592623	T	A
13	28592629	28592629	T	C
13	28592634	28592637	CATG	C
13	28592640	28592640	A	C
13	28592640	28592640	A	T
13	28592641	28592641	T	A
![image](https://user-images.githubusercontent.com/40486459/110201401-286c9980-7e63-11eb-8957-211d678a2a56.png)

