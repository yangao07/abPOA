This Repository
==================

This repository was created because the [Google Code repository](https://code.google.com/archive/p/pbsim/downloads) provided by the original authors is not being maintained, and Google Code is now defunct. 


About PBSIM
==============

PacBio sequencers produce two types of characteristic reads as below.

  Continuous Long Read (CLR)    : long and high error rate.
  Circular consensus Read (CCS) : short and low error rate.

We have developed a PacBio reads simulater (called PBSIM) in which
sampling-based and model-based simulations are implemented.


Building PBSIM
================

To build PBSIM run;

```bash

 autoreconf -i
 ./configure
 make
 ```
 
A new executable pbsim will be available in the src/ directory

Run PBSIM with sample data
=============================

To run model-based simulation:

  pbsim --data-type CLR
        --depth 20
        --model_qc data/model_qc_clr
        sample/sample.fasta

In the example above, simulated read sequences are randomly sampled from
a reference sequence ("sample/sample.fasta") and differences (errors) of
the sampled reads are introduced.
Data type is CLR, and coverage depth is 20.
If the reference sequence is multi-FASTA file, the simulated data is created
for each FASTA. Three output files are created for each FASTA.
"sd_0001.ref" is a single-FASTA file which is copied from the reference
sequence.
"sd_0001.fastq" is a simulated read dataset in the FASTQ format.
"sd_0001.maf" is a list of alignments between reference sequence and
simulated reads in the MAF format.
The length and accuracy of reads are simulated based on our model of PacBio
read.

To run sampling-based simulation:

  pbsim --data-type CLR 
        --depth 20
        --sample-fastq sample/sample.fastq
        sample/sample.fasta

In the sampling-based simulation, read length and quality score are
the same as those of a read taken randomly in the sample PacBio dataset
("sample/sample.fastq").

If you want to create several simulated data with different coverage depths
using the same PacBio sample, you would be better to use --sample-profile-id
option as below. You can save time to parse "sample/sample.fastq".

  (1) storing profile

    pbsim --data-type CLR 
          --depth 20
          --prefix depth20
          --sample-fastq sample/sample.fastq
          --sample-profile-id pf1
          sample/sample.fasta

  (2) reusing profile

    pbsim --data-type CLR 
          --depth 30
          --prefix depth30
          --sample-profile-id pf1
          sample/sample.fasta

    pbsim --data-type CLR 
          --depth 40
          --prefix depth40
          --sample-profile-id pf1
          sample/sample.fasta


Model-based simulation
=========================

For each read, the length is randomly drawn from the log-normal distribution
with given mean and standard deviation.

How to simulate the accuracy of each read is different between CLR and CCS
read. For CLR reads, the accuracy is randomly drawn from the normal
distribution with given mean and standard deviation. For CCS reads,
an exponential function which is fit to the the real distribution is
utilized to simulate with fixed mean and standard deviation.

Errors from single molecule sequencing which generates PacBio reads are
considered to be stochastical, therefore quality scores are randomly chosen
from a frequency table of quality scores (named "quality profile") for each
accuracy of a read. For accuracies of 0-59% and 86-100% of CLR readsi and
0-84% of CCS reads, uniform distributions are used because real PacBio
datasets are not sufficiently large.
"data/model_qc_clr" is quality profile for CLR and 
"data/model_qc_ccs" is for CCS.

Simulated read sequences are randomly sampled from a reference sequence.
The percentage of both directions of reads is same. Differences (errors)
of the sampled reads are introduced as follows.
The substitutions and insertions are introduced according to the quality
scores. Their probabilities are computed for each positions of a simulated
read from the error probability of the position (computed from the quality
score of the position) and the ratios of differences given by the user.
Patterns of substitutions are randomly sampled.  
We observed that inserted nucleotides are often the same as their following
nucleotides. According to the observed bias, half of inserted nucleotides
are chosen to be the same as their following nucleotides, and the other
half are randomly chosen.
The deletion probability is uniform for all positions of all simulated
reads, which is computed from the mean error probability of the read set
and the ratios of differeces. 

By setting minimum and maximum of the length, the range of length chosen
from the distribution model can be restricted. Note that mean and standard
deviation of the chosen length are influenced by this restriction.
The accuracy can be restricted in the same way, however unlike the length,
the restriction of accuracy is not strict, and can be used in only case
of CLR reads.


Sampling-based simulation
============================

The lengths and quality scores of reads  are simulated by randomly
sampling them in a real library of PacBio reads provided by the user.
randomly in a real PacBio dataset given by user. Subsequently, their
nucleotide sequences are simulated by the same method with the model-
based simulation. The restriction of length and accuracy are also
the same as model-based simulation.


Input files
==============

PBSIM requires reference sequences in the single- or multi-FASTA Format. 

A real PacBio read data is required for sampling-based simulation,
specified with the --sample-fastq option.
FASTQ format must be Sanger standard (fastq-sanger).


Output files
===============

If a reference sequence file is multi-FASTA format, simulated datasets
are generated for each reference sequence numbered sequentially.
Three output files are created for each reference sequence.

"sd_<num>.ref" is a single-FASTA file which is copied from the reference
sequence.
"sd_<num>.fastq" is a simulated read dataset in the FASTQ format.
"sd_<num>.maf" is a list of alignments between reference sequence and
simulated reads in the MAF format.

"sd" is prefix which can be specified with the --prefix option.


Quality profile
==================

Quality profiles are derived from frequencies of real quality scores 
for each accuracy of a read. 
"data/model_qc_clr" is quality profile for CLR, "data/model_qc_ccs" is for CCS.
In "data/model_qc_clr", 1st column is accuracies of a read, and 2nd-23th
columns are proportions of phred quality scores (0-21).
In "data/model_qc_ccs", 1st column is accuracies of a read, and 2nd-95th
columns are proportions of phred quality scores (0-93).


Runtime and memory
=====================

When a coverage depth is 100x and a length of reference sequence is about 10M,
PBSIM generates simulated dataset in several minutes.
The runtime is roughly proportional to the coverage depth and the length of
reference sequence. 

PBSIM requires memory of the length of reference sequence plus several mega
bytes. 


Contributors
=====================
@kiwiroy - autotools and warning corrections

@jumpinsky - fixed memory leaks
