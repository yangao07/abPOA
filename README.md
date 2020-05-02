# abPOA: adaptive banded Partial Order Alignment
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/abPOA/total.svg?label=Download)](https://github.com/yangao07/abPOA/releases)
[![Latest Release](https://img.shields.io/github/release/yangao07/abPOA.svg?label=Release)](https://github.com/yangao07/abPOA/releases/latest)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/abpoa.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/abpoa)
[![PyPI](https://img.shields.io/pypi/v/pyabpoa.svg?style=flat)](https://pypi.python.org/pypi/pyabpoa)
[![Build Status](https://img.shields.io/travis/yangao07/abPOA/master.svg?label=Master)](https://travis-ci.org/yangao07/abPOA)
[![License](https://img.shields.io/badge/License-GPL-black.svg)](https://github.com/yangao07/abPOA/blob/master/LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/yangao07/abPOA.svg?label=Issues)](https://github.com/yangao07/abPOA/issues)
<!-- [![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-purple.svg)](https://doi.org/10.1093/bioinformatics/btz376) -->

## Getting started
Download the [latest release](https://github.com/yangao07/abPOA/releases):
```
wget https://github.com/yangao07/abPOA/releases/download/v1.0.1/abPOA-v1.0.1.tar.gz
tar -zxvf abPOA-v1.0.1.tar.gz && cd abPOA-v1.0.1
```
Install via conda and run with test data:
```
conda install -c bioconda abpoa
abpoa ./test_data/seq.fa > cons.fa
```
Or, make from source and run with test data:
```
make; ./bin/abpoa ./test_data/seq.fa > cons.fa
```
## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Installing abPOA](#conda)
  - [Building abPOA from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [General usage](#usage)
  - [Generate consensus sequence](#gen_cons)
  - [Generate row-column multiple sequence alignment](#gen_msa)
  - [Generate plot of alignment graph](#gen_plot)
- [Commands and options](#cmd)
- [Input](#input)
- [Output](#output)
  - [Consensus sequence](#cons)
  - [Row-column multiple sequence alignment](#msa)
  - [Plot of alignment graph](#plot)
- [For development](#dev)
- [Contact](#contact)

## <a name="introduction"></a>Introduction
abPOA is an extended version of [Partial Order Alignment (POA](10.1093/bioinformatics/18.3.452)) that performs adaptive banded dynamic programming (DP) on an SIMD implementation. 
abPOA can perform multiple sequence alignment (MSA) of a set of sequence and generate a consensus sequence applying the [heaviest bundling algorithm](10.1093/bioinformatics/btg109) on the final alignment graph.

The adaptive banded DP and the SIMD vecorization significantly accelerate the time-consuming POA procedure. 
abPOA supports flexible scoring schemes that allow linear, affine and convex gap penalties. 
It right now supports SSE2/SSE4.1/AVX2/AVX512F/AVX512BW vectorization and more advanced instructions will be supported in the future.


## <a name="install"></a>Installation

### <a name="conda"></a>Installing abPOA via conda
On Linux/Unix and Mac OS, abPOA can be installed via
```
conda install -c bioconda abpoa   # install abPOA program
```

### <a name="build"></a>Building abPOA from source files
You can also choose to build abPOA from source files.
It is recommended to download the latest release of abPOA 
from the [release page](https://github.com/yangao07/abPOA/releases).
```
wget https://github.com/yangao07/abPOA/releases/download/v1.0.1/abPOA-v1.0.1.tar.gz
tar -zxvf abPOA-v1.0.1.tar.gz
cd abPOA-v1.0.1; make
```
Or, you can use `git clone` command to download the source code.
This gives you the latest version of abPOA, which might be still under development.
```
git clone https://github.com/yangao07/abPOA.git
cd abPOA; make
```

### <a name="binary"></a>Pre-built binary executable file for Linux/Unix 
If you meet any compiling issue, please try the pre-built binary file:
```
wget https://github.com/yangao07/abPOA/releases/download/v1.0.1/abPOA-v1.0.1_x64-linux.tar.gz
tar -zxvf abPOA-v1.0.1_x64-linux.tar.gz
```

## <a name="usage"></a>General usage
### <a name="gen_cons"></a>Generate consensus sequence

```
abpoa seq.fa > cons.fa
```

### <a name="gen_cons"></a>Generate row-column multiple sequence alignment in PIR format

```
abpoa seq.fa -r2 > cons.out
```

### <a name="gen_plot"></a>Generate plot of alignment graph

```
abpoa seq.fa -g poa.png > cons.fa
```

## <a name="cmd"></a>Commands and options
```
abpoa: adaptive banded Partial Order Alignment

Usage: abpoa [option] <in.fa/fq> > cons.fa/msa.out

Options:
  Alignment:
    -m --aln-mode INT       alignment mode [0]
                              0: global, 1: local, 2: extension
    -M --match    INT       match score [2]
    -X --mismatch INT       mismatch penalty [4]
    -O --gap-open INT(,INT) gap opening penalty (O1,O2) [4,24]
    -E --gap-ext  INT(,INT) gap extension penalty (E1,E2) [2,1]
                            abPOA provides 3 gap penalty modes, penalty of a g-long gap:
                            - convex (default): min{O1+g*E1, O2+g*E2}
                            - affine (set O2 as 0): O1+g*E1
                            - linear (set O1 as 0): g*E1
  Adaptive banded DP:
    -b --extra-b  INT       first part of extra band [10]
                            set b as < 0 to disable adaptive banded DP
    -f --extra-f  FLOAT     second part of extra band: f * L, L is the length of input sequence [0.01]
                            width of extra band is b + f * L
  Input/Output:
    -l --in-list            input file is a list of sequence file [False]
                            each line is one sequence file
    -o --output   FILE      ouput to FILE [stdout]
    -r --result   INT       output result mode [0]
                            - 0: consensus (FASTA format)
                            - 1: MSA (PIR format)
                            - 2: both 0 & 1
    -g --out-pog  FILE      dump final alignment graph to FILE (.pdf/.png) [Null]

    -h --help               print this help usage information
    -v --version            show version number

```

## <a name="input"></a>Input
abPOA works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats. The input file is expected to contains multiple reads which will be processed as a whole set. 

abPOA also can take a list of file names as the input file with option `-l`, where each line is the path to one file containing multiple sequences.

## <a name="output"></a>Output
### <a name="cons"></a>Consensus sequence 
abPOA outputs consensus sequence in FASTA format with the name field as "Consensus_sequence".
For example:
```
>Consensus_sequence
ACGTGTACACGTTGAC
```
### <a name="msa"></a>Row-column multiple sequence alignment
abPOA outputs the row-column multiple sequence alignment of input sequences in PIR format with a FASTA header. For example:
```
>Multiple_sequence_alignment
ACGTGTACA-GTTGAC
A-G-GTACACGTT-AC
A-GTGT-CACGTTGAC
ACGTGTACA--TTGAC
```
The `-` in the sequence stands for alignment gap. 

### <a name="plot"></a>Plot of alignment graph

abPOA can generate a plot of the final partial order alignment graph with the help of `DOT` programs. For example:

![pog](https://github.com/yangao07/abPOA/blob/master/pog.png)

The numbers inside the nodes are the node IDs. The numbers on the edges are the edge weights.
`S` and `E` are virtual start and end nodes that have no sequence bases.

## <a name="dev"></a>For development
abPOA is not only a standalone tool for MSA and consensus calling, it can also work as a programming library. [example.c](example.c) shows how to use the C APIs of abPOA to perform MSA and generate a consensus from a set of sequences. Basically, the library file `libabpoa.a` and two header files [abpoa.h](include/abpoa.h) and [simd_instruction.h](include/simd_instruction.h) are needed to make abPOA work in your program.

abPOA also provides Python bindings of all the C APIs. Refer to [python/README.md](python/README.md) for more details.

## <a name="contact"></a>Contact
Yan Gao yangao07@hit.edu.cn

Yi Xing XINGYI@email.chop.edu

Yadong Wang ydwang@hit.edu.cn

[github issues](https://github.com/yangao07/abPOA/issues)
