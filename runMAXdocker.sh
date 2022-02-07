#!/bin/sh

syntax="bash runMAXdocker.sh -param [parameter input file]"

echo -e "
                                                                                                                                                                           
---------------------------------------------------------------------------
                     You are running MAX v0.2.0 using docker ....
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


printf "CPUNUM=${CPUNUM}\nhgversion=${hgversion}\nuseMAX2=${useMAX2}" > MAX_temp

path=`pwd`
params_path="${path}/MAX_temp"

cmd="sudo docker run -it"
cmd+=" -v ${mutlist}:/source/referenceDB/mutation_list.txt"
cmd+=" -v ${gtffile}:/source/referenceDB/WT_transcript.gtf"
cmd+=" -v ${fastafile}:/source/referenceDB/WT_transcripts.fa"
cmd+=" -v ${INPUT}:/source/input"
cmd+=" -v ${OUTPUT}:/source/output"
cmd+=" -v ${params_path}:/source/params.sh"  
cmd+=" nghiavtr/max:v0.2.0"

eval $cmd
rm MAX_temp
