import pyabpoa as pa

#@parameters of msa_aligner:
#   aln_mode='g' # g: global, l: local, e: extension
#   match=2
#   mismatch=4
#   gap_open1=4
#   gap_open2=24
#   gap_ext1=2
#   gap_ext2=1
#   extra_b = 10 # 1st part of extra band, -1 to disable banded DP
#   extra_f = 0.01  # 2nd part of eatra band. w = extra_b+extra_f*L (L is sequence length)
#   end_bonus=-1 # for extension alignment
#   zdrop=-1     # for extension alignment
#   cons_agrm=0  # 0: heaviest bundling, 1: heaviest column
#   is_diploid=0 # 1: diploid data
#   min_freq=0.3 # minimum frequence of each consensus to output for diploid data

# construct msa aligner
a = pa.msa_aligner()

seqs=[
'CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA',
'CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC',
'CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC',
'CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC',
'CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC',
'CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC',
'CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC',
'CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC',
'CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC',
'CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT'
]

#@parameters of msa
#seqs: multiple sequences
out_cons=True # generate consensus sequence, set as False to disable
out_msa=True  # generate row-column multiple sequence alignment, set as False to disable
out_pog="example.png" # generate plot of alignment graph, set None to disable

# multiple sequence alignment for 'seqs'
res=a.msa(seqs, out_cons=out_cons, out_msa=out_msa, out_pog=out_pog)

# output result
if out_cons:
    print(">Consensus_sequence")
    for seq in res.cons_seq:
        print(seq)

if out_msa:
    print(">Multiple_sequence_alignment")
    res.print_msa()
