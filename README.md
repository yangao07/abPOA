# abPOA: adaptive banded Partial Order Alignment
[![Latest Release](https://img.shields.io/github/release/yangao07/abPOA.svg?label=Release)](https://github.com/yangao07/abPOA/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/abPOA/total.svg?label=Download)](https://github.com/yangao07/abPOA/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/abpoa.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/abpoa)
[![PyPI](https://img.shields.io/pypi/dm/pyabpoa.svg?label=pip%20install)](https://pypi.python.org/pypi/pyabpoa)
[![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963)
[![GitHub Issues](https://img.shields.io/github/issues/yangao07/abPOA.svg?label=Issues)](https://github.com/yangao07/abPOA/issues)
[![Build Status](https://img.shields.io/travis/yangao07/abPOA/master.svg?label=Master)](https://travis-ci.org/yangao07/abPOA)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/abPOA/blob/master/LICENSE)
<!-- [![PyPI](https://img.shields.io/pypi/v/pyabpoa.svg?style=flat)](https://pypi.python.org/pypi/pyabpoa) -->
## Updates (v1.2.1)

- Add minimizer-based seeding, reduce memory usage for long input sequences
- Remove redundant topological sorting

## Getting started
Download the [latest release](https://github.com/yangao07/abPOA/releases):
```
wget https://github.com/yangao07/abPOA/releases/download/v1.2.1/abPOA-v1.2.1.tar.gz
tar -zxvf abPOA-v1.2.1.tar.gz && cd abPOA-v1.2.1
```
Make from source and run with test data:
```
make; ./bin/abpoa ./test_data/seq.fa > cons.fa
```
Or, install via conda and run with test data:
```
conda install -c bioconda abpoa
abpoa ./test_data/seq.fa > cons.fa
```
## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Installing abPOA via conda](#conda)
  - [Building abPOA from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [General usage](#usage)
  - [To generate consensus sequence](#gen_cons)
  - [To generate row-column multiple sequence alignment](#gen_msa)
  - [To generate graph information in GFA format](#gen_gfa)
  - [To align sequence to an existing graph in GFA/MSA format](#aln_to_gfa)
  - [To generate a plot of the alignment graph](#gen_plot)
- [Commands and options](#cmd)
- [Input](#input)
- [Output](#output)
  - [Consensus sequence](#cons)
  - [Row-column multiple sequence alignment](#msa)
  - [Full graph information](#gfa)
  - [Plot of alignment graph](#plot)
- [For development](#dev)
- [Evaluation datasets](#eval)
- [Contact](#contact)

## <a name="introduction"></a>Introduction
abPOA is an extended version of [Partial Order Alignment (POA](10.1093/bioinformatics/18.3.452)) 
that performs adaptive banded dynamic programming (DP) with an SIMD implementation. 
abPOA can perform multiple sequence alignment (MSA) on a set of input sequences and 
generate a consensus sequence by applying the [heaviest bundling algorithm](10.1093/bioinformatics/btg109) 
to the final alignment graph.

abPOA can generate high-quality consensus sequences from error-prone long reads and offer 
significant speed improvement over existing tools.

abPOA supports three alignment modes (global, local, extension) and flexible scoring schemes that allow linear, affine and convex gap penalties. 
It right now supports SSE2/SSE4.1/AVX2/AVX512F/AVX512BW vectorization and more advanced instructions 
will be supported in the future.

For more information, please refer to our [preprint paper](https://doi.org/10.1101/2020.05.07.083196).

## <a name="install"></a>Installation

### <a name="conda"></a>Installing abPOA via conda
On Linux/Unix and Mac OS, abPOA can be installed via
```
conda install -c bioconda abpoa   # install abPOA program
```

### <a name="build"></a>Building abPOA from source files
You can also build abPOA from source files. 
Make sure you have gcc (>=6.4.0) and zlib installed before compiling.
It is recommended to download the [latest release](https://github.com/yangao07/abPOA/releases).
```
wget https://github.com/yangao07/abPOA/releases/download/v1.2.1/abPOA-v1.2.1.tar.gz
tar -zxvf abPOA-v1.2.1.tar.gz
cd abPOA-v1.2.1; make
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
wget https://github.com/yangao07/abPOA/releases/download/v1.2.1/abPOA-v1.2.1_x64-linux.tar.gz
tar -zxvf abPOA-v1.2.1_x64-linux.tar.gz
```

## <a name="usage"></a>General usage
### <a name="gen_cons"></a>To generate consensus sequence

```
abpoa seq.fa > cons.fa
```

### <a name="gen_cons"></a>To generate row-column multiple sequence alignment in PIR format

```
abpoa seq.fa -r2 > cons.out
```

### <a name="gen_gfa"></a>To generate graph information in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format

```
abpoa seq.fa -r3 > out.gfa
```
To include the generated consensus sequence as a path in the GFA file:
```
abpoa seq.fa -r4 > out.gfa
```

### <a name="aln_to_gfa"></a>To align sequence to an existing graph in GFA/MSA format
```
abpoa -i in.gfa seq.fa -r3 > out.gfa
abpoa -i in.msa seq.fa -Ar1 > out.msa
```
For GFA input file, `S` and `P` lines are required and are used to reconstruct the alignment graph.
For MSA input file, which is generally a FASTA format file, `-` in the sequence indicates the alignment gap.
If you want to use abPOA to generate a MSA output file and then perform the incremental graph alignment, please do not forget `-A` to include the FASTA header of each sequence:

```
abpoa seq1.fa -Ar1 > seq1.msa
abpoa -i seq1.msa seq2.fa > cons.fa
```

### <a name="gen_plot"></a>To generate a plot of the alignment graph

```
abpoa seq.fa -g poa.png > cons.fa
```
See [Plot of alignment graph](#plot) for more details about the plot file.

## <a name="cmd"></a>Commands and options
```
Usage: abpoa [options] <in.fa/fq> > cons.fa/msa.out/abpoa.gfa

Options:
  Alignment:
    -m --aln-mode INT       alignment mode [0]
                              0: global, 1: local, 2: extension
    -M --match    INT       match score [2]
    -X --mismatch INT       mismatch penalty [4]
    -O --gap-open INT(,INT) gap opening penalty (O1,O2) [4,24]
    -E --gap-ext  INT(,INT) gap extension penalty (E1,E2) [2,1]
                            abPOA provides three gap penalty modes, cost of a g-long gap:
                            - convex (default): min{O1+g*E1, O2+g*E2}
                            - affine (set O2 as 0): O1+g*E1
                            - linear (set O1 as 0): g*E1
    -s --amb-strand         ambiguous strand mode [False]
                            for each input sequence, try the reverse complement if the current
                            alignment score is too low, and pick the strand with a higher score
  Adaptive banded DP:
    -b --extra-b  INT       first adaptive banding parameter [10]
                            set b as < 0 to disable adaptive banded DP
    -f --extra-f  FLOAT     second adaptive banding parameter [0.01]
                            the number of extra bases added on both sites of the band is
                            b+f*L, where L is the length of the aligned sequence
  Minimizer-based seeding and partition (only effective in global alignment mode):
    -N --no-seeding         disable seeding [False]
    -k --k-mer       INT    minimizer k-mer size [19]
    -w --window      INT    minimizer window size [10]
    -n --min-poa-win INT    min. size of window to perform POA [50]
    -p --progressive        build guide tree and perform progressive partial order alignment [False]
  Input/Output:
    -l --in-list            input file is a list of sequence file names [False]
                            each line is one sequence file containing a set of sequences
                            which will be aligned by abPOA to generate a consensus sequence
    -i --incrmnt  FILE      incrementally align sequences to an existing graph/MSA [Null]
                            graph could be in GFA or MSA format generated by abPOA
    -o --output   FILE      ouput to FILE [stdout]
    -r --result   INT       output result mode [0]
                            - 0: consensus (FASTA format)
                            - 1: MSA (PIR format)
                            - 2: both 0 & 1
                            - 3: graph (GFA format)
                            - 4: graph with consensus path (GFA format)
    -A --msa-header         add read ID as header of each sequence in MSA output [False]
    -g --out-pog  FILE      dump final alignment graph to FILE (.pdf/.png) [Null]

    -h --help               print this help usage information
    -v --version            show version number
```

## <a name="input"></a>Input
abPOA works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats. The input file is 
expected to contains multiple sequences which will be processed sequentially to perform the iterative 
sequence-to-graph (partial order) alignment.

abPOA can also take a list of file names as input with option `-l`, where each line is the path to one 
file containing multiple sequences. Each sequence file is then individually aligned by abPOA to generate a
consensus sequence.

## <a name="output"></a>Output
### <a name="cons"></a>Consensus sequence 
By default, abPOA only outputs the consensus sequence generated from the final alignment agraph.
It is in FASTA format with the name field set as "Consensus_sequence".
For example:
```
>Consensus_sequence
ACGTGTACACGTTGAC
```
### <a name="msa"></a>Row-column multiple sequence alignment
abPOA can also output the row-column multiple sequence alignment (RC-MSA) of all the aligned sequences in PIR format
with an additional FASTA header `>Multiple_sequence_alignment`. For example:
```
>Multiple_sequence_alignment
ACGTGTACA-GTTGAC
A-G-GTACACGTT-AC
A-GTGT-CACGTTGAC
ACGTGTACA--TTGAC
```
The `-` in the sequence stands for alignment gap. 

### <a name="gfa"></a>Full graph information
abPOA can output the final alignment graph in GFA format.
Each segment line (`S` line) represents one node and each link line (`L` line) represents one edge between two nodes.
The original input sequences and the generated consensus sequence are described as paths in `P` lines.

abPOA outputs two graph-related numbers in the header line (`H` line):
`NS` and `NL`, which denote the total number of nodes and edges in the GFA file, respectively.

Please refer to the [GFA specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) for more details of the GFA format.

### <a name="plot"></a>Plot of alignment graph

abPOA can generate a plot of the final partial order alignment graph with the help of `graphviz dot`. 
For example:

![pog](https://github.com/yangao07/abPOA/blob/master/pog.png)

The numbers inside the nodes are the node IDs. The numbers on the edges are the edge weights.
`S` and `E` are the auxiliary start and end nodes that have no sequence bases.

Make sure you have `dot` installed beforing using abPOA to generate the plot.
For Linux/Unix systems: `sudo apt-get install graphviz`.

## <a name="dev"></a>For development
abPOA is not only a stand-alone tool for MSA and consensus calling, it can also work as a programming library. [example.c](example.c) shows how to use the C APIs of abPOA to take a set of sequences as input and perform MSA and consensus calling. Basically, the library file `libabpoa.a` and two header files [abpoa.h](include/abpoa.h) and [simd_instruction.h](include/simd_instruction.h) are needed to make the abPOA library work in your program.

abPOA also provides Python bindings to all the primary C APIs. Refer to [python/README.md](python/README.md) for more details.

## <a name="eval"></a>Evaluation datasets
The evaluation datasets and scripts used in [abPOA paper](https://dx.doi.org/10.1093/bioinformatics/btaa963) can be found in [abPOA-v1.0.5](https://github.com/yangao07/abPOA/releases/tag/v1.0.5).

## <a name="contact"></a>Contact
Yan Gao gaoy286@mail.sysu.edu.cn

Yi Xing xingyi@email.chop.edu

Yadong Wang ydwang@hit.edu.cn

[github issues](https://github.com/yangao07/abPOA/issues)
