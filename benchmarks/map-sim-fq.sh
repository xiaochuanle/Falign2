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

ENZYME_SEQ=(
"^GATC"
"CATG^"
)

ENZYME_NAME=(
dpnii
nlaiii
)

SIM_MAP_DIR=/data3/cy/benchmarks/sim-maps
mkdir -p ${SIM_MAP_DIR}

for ((i=0;i<3;i++))
do
        for ((j=0;j<2;j++))
        do
                echo ${REF[$i]} ${SPECIES[$i]} ${ENZYME_SEQ[$j]} ${ENZYME_NAME[$j]}
                /usr/bin/time -v  ${NGMLR} -t 96 -r ${REF[$i]} -q ${SIM_FQ_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.fastq -o ${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.ngmlr.sam -x ont
                if [ $? -ne 0 ]; then
                        exit 1
                fi
        done
done

for ((i=0;i<3;i++))
do
        for ((j=0;j<2;j++))
        do
                echo ${REF[$i]} ${SPECIES[$i]} ${ENZYME_SEQ[$j]} ${ENZYME_NAME[$j]}
                ${MINIMAP2} -x map-ont -a --eqx -t 96 -o ${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.mm2.sam ${REF[$i]}.mm2.idx ${SIM_FQ_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.fastq
        done
done

for ((i=0;i<3;i++))
do
        for ((j=0;j<2;j++))
        do
                echo ${REF[$i]} ${SPECIES[$i]} ${ENZYME_SEQ[$j]} ${ENZYME_NAME[$j]}
                ${FALIGN2} map -e ${ENZYME_SEQ[$j]} -u frag-bam -o ${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.falign.bam -t 96 ${REF[$i]} ${SIM_FQ_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.fastq
                if [ $? -ne 0 ]; then
                        exit 1
                fi
        done
done

