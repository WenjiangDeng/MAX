#!/bin/bash

syntax="bash runMAX.sh -param [parameter input file]"

echo -e "
                                                                                                                                                                           
---------------------------------------------------------------------------
                     You are running MAX v0.2.0
---------------------------------------------------------------------------
"
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -param|--parameter)
     param=$2
     shift
     ;;

     *)
esac
shift
done

if [[ -n "$param" ]]; then
 echo ""
 echo "Loading parameters from file..."
 if [[ -f "$param" ]]; then
  source $param
 else
  echo "Your parameter file does not exist. Please check and try again!"
  exit
 fi
fi


#source params.sh
#mutlist="/source/referenceDB/mutation_list.txt"
#gtffile="/source/referenceDB/WT_transcript.gtf"
#fastafile="/source/referenceDB/WT_transcripts.fa"
##hgversion="hg19"
##CPUNUM=8
#INPUT="/source/input"
#OUTPUT="/source/output"

cd $OUTPUT

### internal input data
xmatEqFn="raw_Xmatrix.eq"
sampleMutFn="MAX_mut_list_keepMutID.txt"
if ((useMAX2 == 1)); then
 Ycount_MAX2_outdir="MAX2_Ycount"
 if [ ! -d $Ycount_MAX2_outdir ]; then
  mkdir $Ycount_MAX2_outdir
 fi
fi


if [[ ! -f X_matrix.RData || ! -f raw_Xmatrix.eq || ! -f MAX_mut_list_keepMutID.txt || ! -f WT_Mut_reference.fa || ! -d Index_reference ]]; then
#construct the wild-type+Mutant reference, reference index and the X matrix
bash /path/to/constructMutXmatrix.sh -m $mutlist -g $gtffile -r $fastafile -v $hgversion -d $OUTPUT
fi


#Generate the eqclass table and Y count matrix using 8 cores
for fn in $(find $INPUT -type f -name "*_1.fastq.gz")
do
 echo $fn
 endpoint=$(expr $(echo "${#fn}") - 11)
 fn2=$(echo $fn|cut -c 1-$endpoint)
 fn2=$(echo $fn2"_2.fastq.gz")
 echo $fn2
 sampleNameId=$(echo $(basename $fn))
 endpoint=$(expr $(echo "${#sampleNameId}") - 11)
 sampleNameId=$(echo $sampleNameId|cut -c 1-$endpoint)
 echo $sampleNameId

 if [[ ! -f $sampleNameId/eqClass.txt ]]; then
  /path/to/bin/MAX -i Index_reference -l IU -1 <(gunzip -c $fn) -2 <(gunzip -c $fn2) -p $CPUNUM -o $sampleNameId
 fi
 
 if ((useMAX2 == 1)); then
  ### some work for MAX2
  #merge mutation to generate Mutated_eqClass.txt
  Rscript /path/to/R/mergeMutSingleSample.R sampleMut=$sampleMutFn sampleID=$sampleNameId sampleEq=$sampleNameId/eqClass.txt
  #generate ycount
  Rscript /path/to/R/genCountSample.R xmatEq=$xmatEqFn sampleMut=$sampleMutFn sampleID=$sampleNameId sampleEq=$sampleNameId/Mutated_eqClass.txt YcountDir=$Ycount_MAX2_outdir
 fi
done

if ((useMAX2 == 0)); then
 ## MAX1 - Estimation
 # Get Ycount
 Rscript /path/to/R/Create_count_matrix.R workdir=$OUTPUT design.matrix=X_matrix.RData core=$CPUNUM
 # Estimate mutant-allele expression using AEM algorithm with 8 cores
 Rscript /path/to/R/AEM_update_X_beta.R workdir=$OUTPUT design.matrix=X_matrix.RData max.out=MAX_isoform_expression.RData core=$CPUNUM
fi

if ((useMAX2 == 1)); then
 ## MAX2 -  Estimation
 # use AEM - default
 Rscript /path/to/R/estimateBeta_MAX2.R workdir=$Ycount_MAX2_outdir sampleMut=$sampleMutFn out=MAX2_isoform_expression.RData
 if [ -d $Ycount_MAX2_outdir ]; then
  rm -r $Ycount_MAX2_outdir
 fi
fi

