# <a name="S-table-of-contents"></a> Table of Contents

* [Introduction](#S-introduction)
* [Tools](#S-tools)
* [Installation](#S-installation)
* [Usage](#S-usage)
* [Use `u4falign` for manipulating alignment results](#S-u4falign)
* [Difference between `read-SAM` and `frag-SAM`](#S-sam-diff)

# <a name="S-introduction"></a> Introduction

`Falign2` is a sequence alignment toolkit for fragmented and error-prone Pore-C reads.
`Falign2` derives from [`Falign`](https://github.com/xiaochuanle/Falign), which is developped in compliance with the [Dip3D](https://github.com/xiaochuanle/dip3d.git) pipeline. We develop `Falign2` here to make it specified for the mapping purpose.
`Falign2` is written in C and C++ programming language.
`Falign2` is built upon the [`HTSLIB`](https://github.com/samtools/htslib), the [`EDLIB`](https://github.com/Martinsos/edlib) and the [`KSW2`](https://github.com/lh3/ksw2) libraries.

# <a name="S-tools"></a> Tools

Two tools are released together in this toolkit.

* `falign2`. The alignment tool.
* `u4falign2`. A utility tool for benchmarking `Falign2` purpose.

# <a name="S-installation"></a> Installation

There are two different ways to build `Falign2`.

## Install from precompiled executable binary:

``` shell
wget https://github.com/xiaochuanle/Falign2/releases/download/v1.0.0/Falign2_1.0.0_Linux-amd64.tar.gz
tar xzvf Falign2_1.0.0_Linux-amd64.tar.gz
```

## Install from source codes

``` shell
git clone https://github.com/xiaochuanle/Falign2.git
cd Falign2/src/
wget --no-check-certificate -O third-party/htslib-1.19.1.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/1.19.1/htslib-1.19.1.tar.bz2/download
make
```

# <a name="S-usage"></a> Usage

## Download sample reads and reference

```shell
wget https://github.com/xiaochuanle/Falign2/releases/download/v1.0.0/arab-pore-c-sample-data.tar.gz
tar xzvf arab-pore-c-sample-data.tar.gz
```

By decomposing the gzip file `arab-pore-c-sample-data.tar.gz`, we have directory named `arab-pore-c-sample-data/`.

The directory `arab-pore-c-sample-data/` provides sample reference file `reference.fa.gz` and sample reads file `reads.fq.gz`. The reads file is created by the DpnII (`^GATC`) restriction enzyme.

## Step 1: Build index for the reference genome
```shell
./Falign2/Linux-amd64/bin/falign2 index arab-pore-c-sample-data/reference.fa.gz
```

## Step 2: Map reads to the reference
```shell
./Linux-amd64/bin/falign map -e ^GATC -t 4 arab-pore-c-sample-data/reference.fa.gz arab-pore-c-sample-data/reads.fq.gz > map.paf
```

The mapping results are output to the file `map.paf` in PAF format.

## <a name="output-format"></a> Output format

`Falign` supports different output formats (specified by the `-u` option):
* PAF
* SAM
* BAM
* FRAG-SAM
* FRAG-BAM

In each output result (note that every alignment has for offsets: read start, read end, reference start, reference end), `falign` adds the following additional fields:
* `qS:i:` the nearest restiction enzyme site to the read start position
* `qE:i:` the nearest restiction enzyme site to the read end position
* `vS:i:` the nearest restriction enzyme site to the reference start position
* `vE:i:` the nearest restriction enzyme site to the reference end position
* `pi:f:` percentage of identity of the alignment
* `SA:Z:` a homologous map of the fragment

# <a name="S-sam-diff"></a> Difference between `read-SAM` and `frag-SAM`

By convention, in SAM format, a read may contain multiple alignments. Among them, one and only one is desinated as primary and the others are categorized as either secondary (SAM flag 0x100 or 256) or supplementary (SAM flag 0x800 or 2048). Since a Pore-C read always captures multiway chromatin interactions, the aligner outputs multiple mapping results, such as
``` shell
0cd79600-51cf-4255-a6f7-0e9660721e85    2064    2    1794202    1 ...
0cd79600-51cf-4255-a6f7-0e9660721e85    2064    2    1791136    60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85    0    2    1733063    60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85    2064    4    12727850    60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85    2048    2    1713207    60 ...
```
From the results we know that:
1.	The name of the Pore-C read is 0cd79600-51cf-4255-a6f7-0e9660721e85.
2.	Falign2 decomposes the read into five fragments, each of which is assigned by one alignment line.
3.	The third alignment is the primary and the others are all supplementary. Falign2 always outputs no secondary alignments.
4.	In the first, second and the forth alignment, the fragment is mapped to the reference with its reversed strand (SAM flag 0x10 or 16).
For Pore-C reads, the alignments should be classified based on fragments, not on reads. Here we introduce the FRAG-SAM format here. In FRAG-SAM, the alignments above are represented as

``` shell
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:0000000000:0000000000:0001794201    16    2    1794202    1 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:0000000001:0000000000:0001791135    16    2    1791136    60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:0000000002:0000000000:0001733062    0    2    1733063    60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:0000000003:0000000001:0012727849    16    4    12727850    60 ...
0cd79600-51cf-4255-a6f7-0e9660721e85_0000000001:0000000004:0000000000:0001713206    0    2    1713207    60 ...
```


# <a name="S-maintainers"></a>Maintainers

* Chuan-Le Xiao, xiaochuanle@126.com
* Ying Chen, chenying2016@gmail.com

# <a name="S-citation"></a>Citation

# <a name="S-license"></a>License

GPLv3
