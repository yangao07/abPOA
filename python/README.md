# pyabpoa: abPOA Python interface
## Introduction
pyabpoa provides an easy-to-use interface to [abPOA](https://github.com/yangao07/abPOA).

## Installation

### Install pyabpoa with pip

pyabpoa can be installed with pip:

```
pip install pyabpoa
```

### Install pyabpoa from source
Alternatively, you can install pyabpoa from source:
```
git clone https://github.com/yangao07/abPOA.git
cd abPOA
make install_py
```

## Examples
The following code illustrates how to use pyabpoa.
```
import pyabpoa as pa
a = pa.msa_aligner()
seqs=[
'CCGAAGA',
'CCGAACTCGA',
'CCCGGAAGA',
'CCGAAGA'
]
res=a.msa(seqs, out_cons=True, out_msa=True, out_pog='pog.png') # perform multiple sequence alignment 
                                                                # generate a figure of alignment graph to pog.png

for seq in res.cons_seq:
    print(seq)  # print consensus sequence

res.print_msa() # print row-column multiple sequence alignment in PIR format
```
You can also try the example script:
```
python ./python/example.py
```


## APIs

### Class pyabpoa.msa_aligner
```
pyabpoa.msa_aligner(aln_mode='g', ...)
```
This constructs a multiple sequence alignment handler of pyabpoa, it accepts the following arguments:
* **aln_mode**: alignment mode. 'g': global, 'l': local, 'e': extension; default: 'g'
* **match**: match score; default: 2
* **gap_open1**: first gap opening penalty; default: 4
* **gap_ext1**: first gap extension penalty; default: 2
* **gap_open2**: second gap opening penalty; default: 24
* **gap_ext2**: second gap extension penalty; default: 1
* **extra_b**: first part of extra band width; default: 10
* **extra_f**: second part of extra band width; Total extra band width: b+f\*L, L is the sequence lengthl default : 0.01
* **is_diploid**: set as 1 if input is diploid datal default: 0
* **min_freq**: minimum frequency of each consensus to output for diploid datal default: 0.3

```
pyabpoa.msa_aligner.msa(seqs, out_cons, out_msa, out_pog=None)
```
This method performs mutliple sequence alignment and generates
* **consensus sequence** if `out_cons` is set as `True`
* **row-column multiple sequence alignment in PIR format** if `out_msa` is set as `True`
* **plot of alignment graph** if `out_pog` is set as a file name with suffix as `.png` or `.pdf`

### Class pyabpoa.msa_result
```
pyabpoa.msa_result(seq_n, cons_n, cons_len, ...)
```
This class describes the information of the generated consensus sequence and row-column multiple sequence alignment. The returned result of `pyabpoa.msa_aligner.msa()` is an object of this class and it has the following properties:
* **seq_n**: number of input sequences
* **cons_n**: number of generated consensus sequences
* **cons_len**: an array of consensus sequence length
* **cons_seq**: an array of consensus sequence
* **msa_len**: size of each row in the row-column multiple sequence alignment
* **msa_seq**: an array containing `seq_n` rows of the row-column multiple sequence alignment

`pyabpoa.msa_result()` also has a function of `print_msa`. It prints the row-column multiple sequence alignment.
