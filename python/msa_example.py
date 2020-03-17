import pyabPOA as pa

# aln_mode='g'
# gap_mode='c'
# match=2
# mismatch=4
# gap_open1=4
# gap_open2=24
# gap_ext1=2
# gap_ext2=1
band_width=-1 # band width during DP, -1 to disable banded DP
# end_bonus=-1
# zdrop=-1
ystop_iter_n=-1 # ystop heuristic, -1 to disable
# ystop_min_processed_n=5
# ystop_min_frac=0.5,
# cons_agrm=0
# multip=1
# min_freq=0.3

# construct msa aligner
a = pa.msa_aligner(band_width = band_width, ystop_iter_n = ystop_iter_n)

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
out_cons=True
out_msa=True

# multiple sequence alignment for 'seqs'
res=a.msa(seqs, out_cons, out_msa)

# output result
if out_cons:
    print(">Consensus_sequence")
    for seq in res.cons_seq:
        print(seq)

if out_msa:
    print(">Multiple_sequence_alignment")
    res.print_msa()