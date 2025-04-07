#Introduction

In this benchmark, we compare performance of [`minimap2`(2.28-r1221-dirty)](https://github.com/lh3/minimap2), [`NGMLR`(0.2.7)](https://github.com/philres/ngmlr) and `Falign2` (1.0.0) on mapping simulated and experimental Pore-C sequencing reads.

# Genomes

We download three reference genomes [H. sapiens (GRCh38)](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz), [D. melanogaster (Release 6)](http://ftp.flybase.net/releases/FB2021_06/dmel_r6.43/fasta/dmel-all-chromosome-r6.43.fasta.gz) and [A. thaliana (TAIR10)](http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz).
For each reference, we remove the alternative and unknown sequences and retain the major chromosome sequences for benchmarking. All of the commands ares included in the shell script `build-genome.sh`. The genomes are save in the directory `/data3/cy/benchmarks/genome`.

# Genomic index

Before conducting the benchmarks, we build index for each reference genome for each aligner.
The shell script `build-genome-index.sh` builds genomic indices for the three aligners.

# Benchmarks on simulated Pore-C reads

### Datasets

We use three simulated Pore-C datasets created with the DpnII restriction enzyme:
1. [human.dpnii](https://zenodo.org/records/14912381/files/human.dpnii.fastq.gz)
2. [arab.dpnii](https://zenodo.org/records/14912381/files/arab.dpnii.fastq.gz)
3. [dmel.dpnii](https://zenodo.org/records/14912381/files/dmel.dpnii.fastq.gz)

And three simulated Pore-C datasets created with the NlaIII restriction enzyme:
1. [human.nlaiii](https://zenodo.org/records/14912381/files/human.nlaiii.fastq.gz)
2. [arab.nlaiii](https://zenodo.org/records/14912381/files/arab.nlaiii.fastq.gz)
3. [dmel.nlaiii](https://zenodo.org/records/14912381/files/dmel.nlaiii.fastq.gz)

The datasets are saved in the directory `/data3/cy/benchmarks/sim-fqs/`.

### Mapping

In the script `map-sim-fq.sh`, we map the six simulated datasets to their corresponding reference genome using three aligners. The alignment results are saved in directory `/data3/cy/benchmarks/sim-maps`.

### Evaluation

We use the script `eval-sim-fq.sh` to evaluate the performance of the three aligners on simulated datasets. We consider five metrics:
1. Fragment precision
2. Fragment recall
3. Fragment boundary precision
4. Contact precision
5. Contact recall

# Benchmarks on experimental Pore-C reads

### Datasets

We generate in-house four Pore-C datasets created with the DpnII restriction enzyme:
1. [GM12878.dpnii](https://download.cncb.ac.cn/gsa-human/HRA004983/HRR1274419/HRR1274419.fastq.gz)
2. [GM24385.dpnii](https://download.cncb.ac.cn/gsa-human/HRA004983/HRR1274420/HRR1274420.fastq.gz)
3. [arab.dpnii](https://ngdc.cncb.ac.cn/gsa/s/93mlUjl0)
4. [dmel.dpnii](https://ngdc.cncb.ac.cn/gsa/s/753wQzi6)

For `GM12878.dnpii` and `GM24385.dpnii`, we extract the first 10 Gbp bases of Pore-C reads for benchmarks using the following two commands:
``` shell
u4falign2 extract-fastq -s 10g HRR1274419.fastq.gz > gm12878.dpnii.fastq
u4falign2 extract-fastq -s 10g HRR1274420.fastq.gz > gm24385.dpnii.fastq
```

We also download two public released Pore-C reads created with the NlaIII restriction enzyme:
1. [GM12878.nlaiii](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&page_size=10&acc=SRR11589392&display=metadata)
2. [GM24385.nlaiii](https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_PoreC/Dorado_v4/HG002_1_Dorado_v4_R1041_PoreC_400bps_sup.fastq.gz)

The last dataset `GM24385.nlaiii` is sequenced by the latest Oxford Nanopore R10.4.1 chemistry.

The experimental datasets are saved in the directory `/data3/cy/benchmarks/real-fqs/`.

### Mapping

In shell script `map-real-fq.sh` we map each experimental Pore-C dataset to its corresponding reference genome using three aligners.

### Evaluation

We use the shell script `eval-real-fq.sh` to evaluation the performance on six experimental Pore-C datasets. We compare three metrics:
1. Number of Pore-C fragments called by the aligners
2. Proportion of Enzyme-aligned fragment boundaries
3. Number of pairwise contacts called by the aligners