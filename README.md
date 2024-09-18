# abPOA: adaptive banded Partial Order Alignment
[![Latest Release](https://img.shields.io/github/release/yangao07/abPOA.svg?label=Release)](https://github.com/yangao07/abPOA/releases/latest)
[![Github All Releases](https://img.shields.io/github/downloads/yangao07/abPOA/total.svg?label=Download)](https://github.com/yangao07/abPOA/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/abpoa.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/abpoa)
[![PyPI](https://img.shields.io/pypi/dm/pyabpoa.svg?label=pip%20install)](https://pypi.python.org/pypi/pyabpoa)
[![Published in Bioinformatics](https://img.shields.io/badge/Published%20in-Bioinformatics-blue.svg)](https://dx.doi.org/10.1093/bioinformatics/btaa963)
[![GitHub Issues](https://img.shields.io/github/issues/yangao07/abPOA.svg?label=Issues)](https://github.com/yangao07/abPOA/issues)
[![C/C++ CI](https://github.com/yangao07/abPOA/actions/workflows/linux-CI.yml/badge.svg)](https://github.com/yangao07/abPOA/actions/workflows/linux-CI.yml)
[![C/C++ CI](https://github.com/yangao07/abPOA/actions/workflows/macos-CI.yml/badge.svg)](https://github.com/yangao07/abPOA/actions/workflows/macos-CI.yml)
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/abPOA/blob/main/LICENSE)
## Updates (v1.5.3)

- Fix a score matrix bug in pyabpoa
- Add consensus mode for pyabpoa: most frequent base at each pos

## Getting started
Download the [latest release](https://github.com/yangao07/abPOA/releases):
```
wget https://github.com/yangao07/abPOA/releases/download/v1.5.3/abPOA-v1.5.3.tar.gz
tar -zxvf abPOA-v1.5.3.tar.gz && cd abPOA-v1.5.3
```
Make from source and run with test data:
```
# Linux or MacOS with armv8
make; ./bin/abpoa ./test_data/seq.fa > cons.fa

# MacOS with armv7
make armv7=1; ./bin/abpoa ./test_data/seq.fa > cons.fa
```
Or, install via conda and run with test data:
```
conda install -c bioconda abpoa
abpoa ./test_data/seq.fa > cons.fa
```
## Table of Contents

- [abPOA: adaptive banded Partial Order Alignment](#abpoa-adaptive-banded-partial-order-alignment)
  - [Updates (v1.5.3)](#updates-v152)
  - [Getting started](#getting-started)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Installation](#installation)
    - [Installing abPOA via conda](#installing-abpoa-via-conda)
    - [Building abPOA from source files](#building-abpoa-from-source-files)
    - [Pre-built binary executable file for Linux/Unix or MacOS](#pre-built-binary-executable-file-for-linuxunix-or-macos)
  - [General usage](#general-usage)
    - [Generate a consensus sequence](#generate-a-consensus-sequence)
    - [Generate multiple consensus sequences](#generate-multiple-consensus-sequences)
    - [Generate row-column multiple sequence alignment in FASTA format](#generate-row-column-multiple-sequence-alignment-in-fasta-format)
    - [Generate graph information in GFA format](#generate-graph-information-in-gfa-format)
    - [Align sequence to an existing graph in GFA/MSA format](#align-sequence-to-an-existing-graph-in-gfamsa-format)
    - [Generate a consensus sequence for amino acid sequences](#generate-a-consensus-sequence-for-amino-acid-sequences)
    - [Generate a plot of the alignment graph](#generate-a-plot-of-the-alignment-graph)
  - [Input](#input)
  - [Output](#output)
    - [Consensus sequence](#consensus-sequence)
    - [Row-column multiple sequence alignment](#row-column-multiple-sequence-alignment)
    - [Full graph information](#full-graph-information)
    - [Plot of alignment graph](#plot-of-alignment-graph)
  - [Algorithm description](#algorithm-description)
    - [Adaptive banding](#adaptive-banding)
    - [Minimizer-based seeding mode](#minimizer-based-seeding-mode)
    - [Minimizer-based progressive tree](#minimizer-based-progressive-tree)
    - [Multiple consensus sequences](#multiple-consensus-sequences)
  - [For development](#for-development)
  - [Evaluation datasets](#evaluation-datasets)
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
It right now supports SSE2/SSE4.1/AVX2 vectorization.

For more information, please refer to our [paper](https://dx.doi.org/10.1093/bioinformatics/btaa963) published in Bioinformatics.

## Installation

### Installing abPOA via conda
On Linux/Unix and Mac OS, abPOA can be installed via
```
conda install -c bioconda abpoa   # install abPOA program
```

### Building abPOA from source files
You can also build abPOA from source files. 
Make sure you have gcc (>=6.4.0) and zlib installed before compiling.
It is recommended to download the [latest release](https://github.com/yangao07/abPOA/releases).
```
wget https://github.com/yangao07/abPOA/releases/download/v1.5.3/abPOA-v1.5.3.tar.gz
tar -zxvf abPOA-v1.5.3.tar.gz
cd abPOA-v1.5.3; make
```
Or, you can use `git clone` command to download the source code.
This gives you the latest version of abPOA, which might be still under development.
```
git clone --recursive https://github.com/yangao07/abPOA.git
cd abPOA; make
```

### Pre-built binary executable file for Linux/Unix or MacOS
If you meet any compiling issue, please try the pre-built binary file for linux:
```
wget https://github.com/yangao07/abPOA/releases/download/v1.5.3/abPOA-v1.5.3_x64-linux.tar.gz
tar -zxvf abPOA-v1.5.3_x64-linux.tar.gz
```
or for macos:
```
wget https://github.com/yangao07/abPOA/releases/download/v1.5.3/abPOA-v1.5.3_arm64-macos.tar.gz
tar -zxvf abPOA-v1.5.3_arm64-macos.tar.gz
```

## General usage
### Generate a consensus sequence

```
abpoa seq.fa > cons.fa
```

abPOA provides two conensus calling methods: 
* heaviest bundlding (default): the path with largest weight along the partial order graph
* most frequent bases: pick the most common base at each alignment position

Sometimes these two methods will generate different consensus sequences [#67](https://github.com/yangao07/abPOA/issues/67)


To use `most frequent bases` method:
```
abpoa seq.fa -a1 > cons.fa
```

### Generate multiple consensus sequences

```
abpoa heter.fa -d2 > 2cons.fa
```

### Generate row-column multiple sequence alignment in FASTA format

```
abpoa seq.fa -r1 > out.msa
abpoa seq.fa -r2 > out_cons.msa
```

### Generate graph information in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format

```
abpoa seq.fa -r3 > out.gfa
```
To include the generated consensus sequence as a path in the GFA file:
```
abpoa seq.fa -r4 > out.gfa
```

### Align sequence to an existing graph in GFA/MSA format
```
abpoa -i in.gfa seq.fa -r3 > out.gfa
abpoa -i in.msa seq.fa -r1 > out.msa
```
For GFA input file, `S` and `P` lines are required and are used to reconstruct the alignment graph.
For MSA input file, which is generally a FASTA format file, `-` in the sequence indicates the alignment gap.
```
abpoa seq1.fa -r1 > seq1.msa
abpoa -i seq1.msa seq2.fa > cons.fa
```

### Generate a consensus sequence for amino acid sequences
```
abpoa -c -t BLOSUM62.mtx input_aa.fa > output_aa_cons.fa
```
abPOA provides two score matrix files for amino acid sequences: `BLOSUM62.mtx`, `HOXD70.mtx`.

You can also use any score matrix, as long as it has the same format as the above two.

### Generate a plot of the alignment graph

```
abpoa seq.fa -g poa.png > cons.fa
```
See [Plot of alignment graph](#plot-of-alignment-graph) for more details about the plot file.

## Input
abPOA works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats. The input file is 
expected to contains multiple sequences which will be processed sequentially to perform the iterative 
sequence-to-graph (partial order) alignment.

abPOA can also take a list of filenames as input with option `-l`, where each line is the path to one 
file containing multiple sequences. Each sequence file is then individually aligned by abPOA to generate a
consensus sequence.

## Output
### Consensus sequence 
By default, abPOA only outputs the consensus sequence generated from the final alignment graph.
It is in FASTA format with the name field set as "Consensus_sequence".
For example:
```
>Consensus_sequence
ACGTGTACACGTTGAC
```

For diploid input sequences, you may want to generate two or more consensus sequences, simply set `-d/--max-num-cons` as a desired value:
```
abpoa heter.fa -d2
```
and this gives you two consensus sequences:
```
>Consensus_sequence_1
CCATTCCCACCATCCTTACCATCAACATCACCATCCCCACCATCCCCAACACCATTCCCACCATCCCTACCATCACCATCACCATCCCCACCAACATCCCCACCACCATCCTCACTACCATCCCCACCACCATTTCCACCATTCCCACCACAGTCACCATCACCCCCACCATCCCCATCATCATCCGCACCATCCCCACCATCCCCACCACCATCTCCATTACCATCCCCACCACCATCTCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCTCGTTACCATCCCCACCACCTTTTCCACCATTCCCACCATCTCCAACACCTCCCCCACCATCATCCCCACCATCCCCACCACCTTCTCCACCATCATTCTCACCATCCCCACCACCATCTCCACCACCATTCTCACCATCTCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCATCATCCCCACCATCC
>Consensus_sequence_2
CCATTCCCACCATCCTTACCATCAACATCACCATCCCCACCATCCCCAACACCATTCCCACCATCCCTACCATCACCATCACCATCCCCACCAACATCCCCACCACCATCCTCACTACCATCCCCACCACCATTTCCACCATTCCCACCACAGTCACCATCACCCCCACCATCCCCATCATCATCCGCACCATCCCCACCATCCCCACCACCATCTCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCTCGTTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCCCATTACCATCCCCACCACCATCCCCATTACCATCCCCACCACCATTTCCACCATTCCCACCATCATCCCCACCACCATCCTCGTTACCATCCCCACCACCTTTTCCACCATTCCCACCATCATCCCCACCGCCATCCTCGTTACCATCCCCACCACCTTTTCCACCATTCCCACCATCTCCAACACCTCCCCCACCATCATCCCCACCATCCCCACCACCTTCTCCACCATCATTCTCACCATCCCCACCACCATCTCCACCACCATTCTCACCATCTCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCAACATCCCCACCATCCCCACCCCCATGCCCACCATCATCCCCACCATCC
```
### Row-column multiple sequence alignment
abPOA can also output the row-column multiple sequence alignment (RC-MSA) of all the aligned sequences in FASTA format.
For example:
```
>1
ACGTGTACA-GTTGAC
>2
A-G-GTACACGTT-AC
>3
A-GTGT-CACGTTGAC
>4
ACGTGTACA--TTGAC
```
The `-` in the sequence stands for alignment gap. 

### Full graph information
abPOA can output the final alignment graph in GFA format.
Each segment line (`S` line) represents one node and each link line (`L` line) represents one edge between two nodes.
The original input sequences and the generated consensus sequence are described as paths in `P` lines.

abPOA outputs two graph-related numbers in the header line (`H` line):
`NS` and `NL`, which denote the total number of nodes and edges in the GFA file, respectively.

Please refer to the [GFA specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) for more details of the GFA format.

### Plot of alignment graph

abPOA can generate a plot of the final partial order alignment graph with the help of `graphviz dot`. 
For example:

![pog](https://github.com/yangao07/abPOA/blob/main/pog.png)

The numbers inside the nodes are the node IDs. The numbers on the edges are the edge weights.
`S` and `E` are the auxiliary start and end nodes that have no sequence bases.

Make sure you have `dot` installed before using abPOA to generate the plot.
For Linux/Unix systems: `sudo apt-get install graphviz`.

## Algorithm description
### Adaptive banding
To understand how the adaptive banding working, please refer to our [Bioinformatics paper](https://dx.doi.org/10.1093/bioinformatics/btaa963).

### Minimizer-based seeding mode
As abPOA always allocates quadratic size of memory, for very long input sequences (>10 kb), memory usage will be a challenge.

To solve this issue, we develop a minimizer-based seeding and partition method to split the sequence and graph with a small window.
The full POA DP matrix can be split into several smaller ones and adaptive banded POA can be performed within each small window separately.

In more detail, abPOA extracts all the minimizers from all the input sequences, then all the minimizer hits between each pair of two sequences can be found.
For each pair of sequences, the minimizer hits are first chained together using relatively stringent criteria to make sure that no big gap exists in the chain.
This usually leads to several separated local chains of minimizer hits.
A second round of chaining is then performed on all the local minimizer chains to generate a global chain going through the entire sequence.
With this global chain, abPOA selects a series of minimizer hits as partition anchors which has at least a distance of 500 bp (by default, -n/--min-poa-win).
Within each partitioned window, abPOA performs banded partial order alignment separately and combines all the alignment results at the end.

### Minimizer-based progressive tree
Instead of aligning all the sequences in the original order, abPOA can alternatively build a progressive tree to guide the alignment order.
The generation of the progressive tree is also based on minimizers.
For each pair of sequences, abPOA calculates their similarity score which is the Jaccard similarity of the minimizers, i.e. the number of minimizer hits divided by the total number of all minimizers from the two sequences.
With all the similarity scores (minimizer-based Jaccard similarity), abPOA builds the progressive tree in the following way:

1. Pick the first two sequences that have the highest scores. The progressive tree set is initialized as these first two sequences.
2. For each remaining sequence, sum the scores between the remaining sequence and all the sequences from the current progressive tree set. Pick the one with the highest sum score, and push it to the progressive tree set.
3. Repeat step 2, until no sequence remains.

Then, abPOA performs partial order alignment following the order of sequences in this progressive tree set.

### Multiple consensus sequences
abPOA supports generating multiple consensus sequences from the final alignment graph (set -d/--max-num-cons as >1).

The general underlying idea is to group input sequences into multiple clusters based on the heterozygous bases in the graph,
Then, one consensus sequence is separately generated for each cluster of input sequences.
The minimum allele frequency for each heterozygous base is 0.25 (by default, -q/--min-freq). 

## For development
abPOA is not only a stand-alone tool for MSA and consensus calling, it can also work as a programming library. [example.c](example.c) shows how to use the C APIs of abPOA to take a set of sequences as input and perform MSA and consensus calling. Basically, the library file `libabpoa.a` and two header files [abpoa.h](include/abpoa.h) and [simd_instruction.h](include/simd_instruction.h) are needed to make the abPOA library work in your program.

abPOA also provides Python bindings to all the primary C APIs. Refer to [python/README.md](python/README.md) for more details.

## Evaluation datasets
The evaluation datasets and scripts used in [abPOA paper](https://dx.doi.org/10.1093/bioinformatics/btaa963) can be found in [abPOA-v1.0.5](https://github.com/yangao07/abPOA/releases/tag/v1.0.5).

## Contact
Yan Gao yangao@ds.dfci.harvard.edu

Yi Xing xingyi@chop.edu

Yadong Wang ydwang@hit.edu.cn

[github issues](https://github.com/yangao07/abPOA/issues)
