#!/bin/bash
mut="";gtf="";ref="";workdir=$PWD;hg="hg19"
while getopts ":m:g:r:d:v:" opt
do
    case $opt in
		v)
        echo "hg version is:$OPTARG"
		hg=$OPTARG
        ;;
        m)
        echo "Mutation file is:$OPTARG"
		mut=$OPTARG
        ;;
        g)
        echo "GFT file is:$OPTARG"
		gtf=$OPTARG
        ;;
        r)
        echo "Wild-type reference is:$OPTARG"
		ref=$OPTARG
        ;;
		d)
        echo "Working directory is:$OPTARG"
		workdir=$OPTARG
        ;;
		:)                        
		echo "$varname" 
		echo "the option -$OPTARG requires an arguement"    
		exit 1
		;;
        ?)
        echo "Unkonwn parameter"
        exit 1;;
    esac
done

cd $workdir

Rscript /path/to/R/MAX.R mut=$mut gtf=$gtf ref=$ref workdir=$workdir hg=$hg


cat Mutant3.fa $ref > WT_Mut_tx_ref.1.fa


awk '{if($1~/^>/)printf("%s Gene_Num=G%s\n",$1,NR);if($1!~/^>/)print $0;}' WT_Mut_tx_ref.1.fa > WT_Mut_tx_ref.2.fa

awk '/^>/{f=!d[$1];d[$1]=1}f' WT_Mut_tx_ref.2.fa > WT_Mut_tx_ref.final.fa

rm WT_Mut_tx_ref.2.fa Mutant3.fa Mutant.fa Mutant.fa.1 Mutant2.fa WT_Mut_tx_ref.1.fa  tmp_gtf.sqlite
echo -e "\n WT_Mut.fa is generated. \n"

mv WT_Mut_tx_ref.final.fa WT_Mut_reference.fa

######################################
################### Generate X matrix
######################################
# export LD_LIBRARY_PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/XAEM-binary-0.1.1/lib:$LD_LIBRARY_PATH
# export PATH=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/XAEM-binary-0.1.1/bin:$PATH
# XAEM_home=/proj/snic2020-6-4/nobackup/private/wenjiang/ASE2020/XAEM-binary-0.1.1
# module load bioinfo-tools
# module load R_packages/4.0.0

chmod +x /path/to/bin/*

cd $workdir
# generate RNA-seq using polyester
Rscript /path/to/R/genPolyesterSimulation.R WT_Mut_reference.fa $workdir
# Index
/path/to/bin/TxIndexer -t WT_Mut_reference.fa -o Index_reference --force 

# generate eqClass table using GenTC
/path/to/bin/GenTC -i Index_reference -l IU -1 sample_01_1.fasta -2 sample_01_2.fasta -p 16 -o $workdir


Rscript /path/to/R/buildCRP.R in=Mutated_Combined_eqclass.txt out=$workdir/X_matrix.RData workdir=$workdir

#keep eqclass of X_matrix
mv eqClass.txt raw_Xmatrix.eq

rm tmp_*.fa
rm tmp*RData
rm *fasta
rm Mutated_Combined_eqclass.txt fragmentInfo.txt gene_tx_tmp.RData wild_type_subset_genes.fa
rm -rf LogDir
######################################
################### Generate X matrix
######################################
