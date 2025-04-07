
# Introduction

This is a computational pipeline for generating simulated Pore-C reads. The simulated reads are used for evaluating the performance of Pore-C aligners.

<font size="16">**This whole pipeline works only with CPU threads. No GPU accelerator is required.**</font>

# Prerequisite

`pore-s-sim` requires [DeepSimulator](https://github.com/lykaust15/DeepSimulator.git) to add Oxford Nanopore-specific sequencing errors.

# Installation

```shell
git clone https://github.com/lykaust15/DeepSimulator.git
cp deep_simulator.PoreC.sh DeepSimulator
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.1.5_linux64.tar.gz
tar xzvf ont-guppy-cpu_3.1.5_linux64.tar.gz -C DeepSimulator/base_caller/guppy_3.1.5/
conda env create -f tensorflow_cdpm.yml
conda env create -f basecall.yml
```

The second to the last command creates a conda environment `tensorflow_cdpm`. It installs only the CPU version of `tensorflow`, don't panic. The last command creates a conda environment `basecall`.

# Usage

Edit the shell script `pore-c-sim.sh` and fill in the relative information:

```shell
#!/bin/bash

# 1) full path of DeepSimulator
P_DSIM=/data2/chenying/pore-c-sim/DeepSimulator

# 2) full path of pore-c-sim
P_PSIM=/data2/chenying/pore-c-sim/pore-c-sim

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

falign sim ${REF} ${ENZYME} ${WRK}/perfect.fa
if [ $? -ne 0 ]; then 
    exit 1
fi 

python ${P_PSIM}/add-errors.py ${P_PSIM}/add-errors-for-one-batch.sh ${P_DSIM}/deep_simulator.PoreC.sh ${WRK}/perfect.fa ${WRK} ${READS}
if [ $? -ne 0 ]; then 
    exit 1
fi 
```

And run this script in your terminal:
```shell
./pore-s-sim.sh
```
