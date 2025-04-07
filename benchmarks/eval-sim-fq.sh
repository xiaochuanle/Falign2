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
		echo "*****************************************"
                echo ${REF[$i]} ${SPECIES[$i]} ${ENZYME_SEQ[$j]} ${ENZYME_NAME[$j]}
		SAM=${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.ngmlr.sam
		FSAM=${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.ngmlr.fixed.sam
		BAM=${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.ngmlr.bam
		FQ=${SIM_FQ_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.fastq

		${U4FALIGN2} fix-ngmlr-sam ${SAM} ${FSAM}
		if [ $? -ne 0 ]; then
			exit 1
		fi

                ${U4FALIGN2} bam-to-frag-bam ${REF[$i]} ${FQ} ${ENZYME_SEQ[$j]} ${FSAM} ${BAM}
                if [ $? -ne 0 ]; then
                        exit 1
                fi

                ${U4FALIGN2} eval-sim-frag-bam ${FQ} ${BAM}
                if [ $? -ne 0 ]; then
                        exit 1
                fi
		echo ${BAM}
        done
done

for ((i=0;i<3;i++))
do
        for ((j=0;j<2;j++))
        do
                echo "*****************************************"
                echo ${REF[$i]} ${SPECIES[$i]} ${ENZYME_SEQ[$j]} ${ENZYME_NAME[$j]}
                SAM=${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.mm2.sam
                BAM=${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.mm2.bam
                FQ=${SIM_FQ_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.fastq

                ${U4FALIGN2} bam-to-frag-bam ${REF[$i]} ${FQ} ${ENZYME_SEQ[$j]} ${SAM} ${BAM}
                if [ $? -ne 0 ]; then
                        exit 1
                fi

                ${U4FALIGN2} eval-sim-frag-bam ${FQ} ${BAM}
                if [ $? -ne 0 ]; then
                        exit 1
                fi
                echo ${BAM}
        done
done

for ((i=0;i<3;i++))
do
        for ((j=0;j<2;j++))
        do
                echo "*****************************************"
                echo ${REF[$i]} ${SPECIES[$i]} ${ENZYME_SEQ[$j]} ${ENZYME_NAME[$j]}
                BAM=${SIM_MAP_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.falign.bam
                FQ=${SIM_FQ_DIR}/${SPECIES[$i]}.${ENZYME_NAME[$j]}.fastq

                ${U4FALIGN2} eval-sim-frag-bam ${FQ} ${BAM}
                if [ $? -ne 0 ]; then
                        exit 1
                fi
                echo ${BAM}
        done
done
