#!/bin/bash

U4FALIGN2=/data3/cy/benchmarks/aligner/Falign2/Linux-amd64/bin/u4falign2

GENOME_DIR=/data3/cy/benchmarks/genome
mkdir -p ${GENOME_DIR}

PATH1="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"

PATH2="http://ftp.flybase.net/releases/FB2021_06/dmel_r6.43/fasta/dmel-all-chromosome-r6.43.fasta.gz"

PATH3="http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"

wget -O ${GENOME_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz ${PATH1}
if [ $? -ne 0 ]; then
	exit 1
fi

wget -O ${GENOME_DIR}/dmel-all-chromosome-r6.43.fasta.gz ${PATH2}
if [ $? -ne 0 ]; then
	exit 1
fi

wget -O ${GENOME_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz ${PATH3}
if [ $? -ne 0 ]; then
	exit 1
fi

CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
echo ${CHR[@]}
${U4FALIGN2} extract-chr ${GENOME_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz ${CHR[@]} > ${GENOME_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.major.fna
if [ $? -ne 0 ]; then
	exit 1
fi

CHR=(1 2 3 4 5)
echo ${CHR[@]}
${U4FALIGN2} extract-chr ${GENOME_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz ${CHR[@]} > ${GENOME_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.major.fa
if [ $? -ne 0 ]; then
	exit 1
fi

CHR=(2L 2R 3L 3R 4 X Y)
echo ${CHR[@]}
${U4FALIGN2} extract-chr ${GENOME_DIR}/dmel-all-chromosome-r6.43.fasta.gz ${CHR[@]} > ${GENOME_DIR}/dmel-all-chromosome-r6.43.major.fasta
if [ $? -ne 0 ]; then
	exit 1
fi
