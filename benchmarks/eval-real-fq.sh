#!/bin/bash

REAL_FQ_DIR=/data3/cy/benchmarks/real-fqs/

REAL_FQ=(
${REAL_FQ_DIR}/gm12878.dpnii.fq
${REAL_FQ_DIR}/gm24385.dpnii.fq
${REAL_FQ_DIR}/SRR11589392.fastq.gz
${REAL_FQ_DIR}/HG002_1_Dorado_v4_R1041_PoreC.fastq.gz
${REAL_FQ_DIR}/20220228-BNP2194-P7-PAK55062.pass.fastq.gz
${REAL_FQ_DIR}/20220427-NPL4197-P6-PAI99739.pass.fastq.gz
)

ENZYME_SEQ=(
"^GATC"
"^GATC"
"CATG^"
"CATG^"
"^GATC"
"^GATC"
)

ENZYME_NAME=(
dpnii
dpnii
nlaiii
nlaiii
dpnii
dpnii
)

DATA=(
gm12878
gm24385
gm12787
gm24385
arab
dmel
)

REF=(
/data3/cy/benchmarks/genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.major.fna
/data3/cy/benchmarks/genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.major.fna
/data3/cy/benchmarks/genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.major.fna
/data3/cy/benchmarks/genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.major.fna
/data3/cy/benchmarks/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.major.fa
/data3/cy/benchmarks/genome/dmel-all-chromosome-r6.43.major.fasta
)

NGMLR=/data3/cy/benchmarks/aligner/ngmlr-0.2.7/ngmlr
FALIGN2=/data3/cy/benchmarks/aligner/Falign2/Linux-amd64/bin/falign2
U4FALIGN2=/data3/cy/benchmarks/aligner/Falign2/Linux-amd64/bin/u4falign2
MINIMAP2=/data3/cy/benchmarks/aligner/minimap2/minimap2

REAL_MAP_DIR=/data3/cy/benchmarks/real-maps

for ((i=0;i<6;i++))
do
	SAM=${REAL_MAP_DIR}/${DATA[$i]}.${ENZYME_NAME[$i]}.ngmlr.sam
	FSAM=${REAL_MAP_DIR}/${DATA[$i]}.${ENZYME_NAME[$i]}.ngmlr.fixed.sam
	BAM=${REAL_MAP_DIR}/${DATA[$i]}.${ENZYME_NAME[$i]}.ngmlr.bam

	echo "*************************************"
	echo ${BAM}
	continue

	${U4FALIGN2} fix-ngmlr-sam ${SAM} ${FSAM}
	if [ $? -ne 0 ]; then
		exit 1
	fi

        ${U4FALIGN2} bam-to-frag-bam ${REF[$i]} ${REAL_FQ[$i]} ${ENZYME_SEQ[$i]} ${FSAM} ${BAM}
        if [ $? -ne 0 ]; then
                exit 1
        fi

        ${U4FALIGN2} frag-bam-stats ${ENZYME_SEQ[$i]} ${BAM}
        if [ $? -ne 0 ]; then
                exit 1
        fi
done

for ((i=0;i<6;i++))
do

        SAM=${REAL_MAP_DIR}/${DATA[$i]}.${ENZYME_NAME[$i]}.mm2.sam
        BAM=${REAL_MAP_DIR}/${DATA[$i]}.${ENZYME_NAME[$i]}.mm2.bam

        echo "*************************************"
	echo ${BAM}
	continue

        ${U4FALIGN2} bam-to-frag-bam ${REF[$i]} ${REAL_FQ[$i]} ${ENZYME_SEQ[$i]} ${SAM} ${BAM}
        if [ $? -ne 0 ]; then
                exit 1
        fi

        ${U4FALIGN2} frag-bam-stats ${ENZYME_SEQ[$i]} ${BAM}
        if [ $? -ne 0 ]; then
                exit 1
        fi
done

for ((i=0;i<6;i++))
do
	BAM=${REAL_MAP_DIR}/${DATA[$i]}.${ENZYME_NAME[$i]}.falign.bam
	echo "***********************************"
        echo ${BAM}

        ${U4FALIGN2} frag-bam-stats ${ENZYME_SEQ[$i]} ${BAM}
        if [ $? -ne 0 ]; then
                exit 1
        fi
        echo ${BAM}
done
