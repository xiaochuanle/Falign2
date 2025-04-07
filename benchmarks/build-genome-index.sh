#!/bin/bash

REF=(
/data3/cy/benchmarks/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.major.fa
/data3/cy/benchmarks/genome/dmel-all-chromosome-r6.43.major.fasta
/data3/cy/benchmarks/genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.major.fna
)

SIM_FQ_DIR=/data3/cy/benchmarks/sim-fqs/
REAL_FQ_DIR=/data3/cy/benchmarks/real-fqs/

NGMLR=/data3/cy/benchmarks/aligner/ngmlr-0.2.7/ngmlr
FALIGN2=/data3/cy/benchmarks/aligner/Falign2/Linux-amd64/bin/falign2
U4FALIGN2=/data3/cy/benchmarks/aligner/Falign2/Linux-amd64/bin/u4falign2
MINIMAP2=/data3/cy/benchmarks/aligner/minimap2/minimap2

SPECIES=(
arab
dmel
human
)

for ((i=0;i<3;i++))
do
	${NGMLR} -t 96 -r ${REF[$i]} -q ${SIM_FQ_DIR}/${SPECIES[$i]}.dpnii.fastq -o ${SPECIES[$i]}.ngmlr.sam -x ont
	if [ $? -ne 0 ]; then
		exit 1
	fi
	rm -f ${SPECIES[$i]}.ngmlr.sam
done

for ((i=0;i<3;i++))
do
	echo ${REF[$i]}
	${FALIGN2} index ${REF[$i]}
	if [ $? -ne 0 ]; then
		exit 1
	fi
done

for ((i=0;i<3;i++))
do
        echo ${REF[$i]}
        ${MINIMAP2} -x map-ont -t 96 -d ${REF[$i]}.mm2.idx ${REF[$i]}
	if [ $? -ne 0 ]; then
		exit 1
	fi
done

