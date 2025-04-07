#!/bin/bash

# 1) full path of DeepSimulator
P_DSIM=/data2/chenying/pore-c-sim/DeepSimulator

# 2) full path of pore-c-sim
P_PSIM=/data1/chenying/dip3d-1/test3/Falign2/pore-c-sim

# 3) full path of reference genome
REF=/data2/chenying/pore-c-sim/genome/dmel-all-chromosome-r6.43.major.fasta

# 4) Restrict enzyme
ENZYME="^GATC"
#ENZYME="CATG^"

# 5) output of simulated pore-c-reads
READS=dmel.dpnii.fastq
#READS=dmel.nlaiii.fastq

#######################

WRK=wrk-pore-c-sim

if [ -f ${WRK} ]; then 
    echo "Remove directory ${WRK}"
    rm -rf ${WRK}
fi 
mkdir -p ${WRK}

u4falign2 sim ${REF} ${ENZYME} ${WRK}/perfect.fa
if [ $? -ne 0 ]; then 
    exit 1
fi 

python ${P_PSIM}/add-errors.py ${P_PSIM}/add-errors-for-one-batch.sh ${P_DSIM}/deep_simulator.PoreC.sh ${WRK}/perfect.fa ${WRK} ${READS}
if [ $? -ne 0 ]; then 
    exit 1
fi 