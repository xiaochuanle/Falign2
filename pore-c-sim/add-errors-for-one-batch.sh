#!/bin/bash

if [ $# -ne 4 ]; then
    echo "USAGE:"
    echo "$0 ds-path FA-DIR FQ-DIR NREADS"
    exit 1
fi

DS=$1
WRK=$2
FQDIR=$3
NREADS=$4

echo "DeepSimulator: ${DS}"
echo "FA-dir: ${WRK}"
echo "FQ-dir: ${FQDIR}"
echo "NUM-reads: ${NREADS}"

for ((i=0;i<${NREADS};i++))
do
    IDLIST[$i]=$i 
done 

function process_one_read()
{
    id=$1

    fa=${WRK}/${id}.fasta
    wrk=${WRK}/w${id}
    rm -rf ${wrk}
    mkdir -p ${wrk}

    ${DS} -i ${fa} -o ${wrk} -B 2
    mv ${wrk}/pass.fastq ${FQDIR}/${id}.fastq

    rm -rf ${wrk}
}

export -f process_one_read
export DS
export WRK
export FQDIR

parallel -j48 process_one_read ::: ${IDLIST[@]}
