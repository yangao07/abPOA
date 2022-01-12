import re, sys, os
from libc.stdlib cimport malloc, free
from libc.stdint cimport uint8_t
from collections import defaultdict as dd
cimport cython
from cabpoa cimport *


cdef class msa_result:
    cdef int seq_n
    cdef int cons_n
    cdef cons_len, cons_seq # _cons_len:[int], _cons_seq:['']
    cdef int msa_len
    cdef msa_seq # _msa_seq:['']

    def __cinit__(self, seq_n, cons_n, cons_len, cons_seq, msa_len, msa_seq):
        self.seq_n = seq_n
        self.cons_n = cons_n
        self.cons_len = cons_len
        self.cons_seq = cons_seq
        self.msa_len = msa_len
        self.msa_seq = msa_seq

    @property
    def seq_n(self): return self.seq_n

    @property
    def cons_n(self): return self.cons_n

    @property
    def cons_len(self): return self.cons_len

    @property
    def cons_seq(self): return self.cons_seq

    @property
    def msa_len(self): return self.msa_len

    @property
    def msa_seq(self): return self.msa_seq

    def print_msa(self): 
        if not self.msa_seq: return
        for s in self.msa_seq:
            print(s)
        return 


def set_seq_int_dict(m):
    if m == 5: # ACGTN ==> 01234, U ==> 4
        seqs = 'ACGUTN'
        ints = [0, 1, 2, 3, 3, 4]
    elif m == 27: # ACGTN    ==> 01234, BDEFH... ==> 56789...
        seqs = 'ACGTNBDEFHIJKLMOPQRSUVWXYZ*'
        ints = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
    else:
        raise Exception('Unexpected m: {}'.format(m))

    seq2int_dict = dd(lambda: m-1)
    int2seq_dict = dd(lambda: '-')
    for s, i in zip(seqs, ints):
        seq2int_dict[s] = i
        seq2int_dict[s.lower()] = i
        int2seq_dict[i] = s
    return seq2int_dict, int2seq_dict

cdef class msa_aligner:
    cdef abpoa_t *ab
    cdef abpoa_para_t abpt
    cdef seq2int_dict, int2seq_dict

    def __cinit__(self, aln_mode='g', is_aa=False, match=2, mismatch=4, score_matrix=b'', gap_open1=4, gap_open2=24, gap_ext1=2, gap_ext2=1,
            extra_b=10, extra_f=0.01, end_bonus=-1, zdrop=-1, cons_agrm=ABPOA_HB, is_diploid=0, min_freq=0.3):
        self.ab = abpoa_init()

        if aln_mode == 'g':
            self.abpt.align_mode = ABPOA_GLOBAL_MODE
        elif aln_mode == 'l':
            self.abpt.align_mode = ABPOA_LOCAL_MODE
        elif aln_mode == 'e':
            self.abpt.align_mode = ABPOA_EXTEND_MODE
        else:
            raise Exception('Unknown align mode: {}'.format(aln_mode))
        if is_aa:
            self.abpt.m = 27
            self.abpt.mat = <int*>malloc(27 * 27 * cython.sizeof(int))
        else:
            self.abpt.m = 5
            self.abpt.mat = <int*>malloc(25 * cython.sizeof(int))
        self.abpt.match = match
        self.abpt.mismatch = mismatch

        if score_matrix:
            if isinstance(score_matrix, str): 
                score_matrix = bytes(score_matrix, 'utf-8')
            if os.path.exists(score_matrix.decode('utf-8')):
                abpoa_set_mat_from_file(&self.abpt, score_matrix)
            else:
                raise Exception('Matrix file not exist: {}'.format(score_matrix.decode('utf-8')))

        self.abpt.gap_open1 = gap_open1
        self.abpt.gap_open2 = gap_open2
        self.abpt.gap_ext1 = gap_ext1
        self.abpt.gap_ext2 = gap_ext2
        self.abpt.ret_cigar = 1

        self.abpt.wb = extra_b
        self.abpt.wf = extra_f 
        self.abpt.end_bonus = end_bonus
        self.abpt.zdrop = zdrop

        self.abpt.cons_agrm = cons_agrm
        self.abpt.is_diploid = is_diploid 
        self.abpt.min_freq = min_freq

        self.seq2int_dict, self.int2seq_dict = set_seq_int_dict(self.abpt.m)

    def __dealloc__(self):
        free(self.abpt.mat)
        abpoa_free(self.ab)

    def __bool__(self):
        return self.ab != NULL


    def msa(self, seqs, out_cons, out_msa, out_pog=b'', incr_fn=b''):
        cdef int seq_n = len(seqs)
        cdef int exist_n = 0
        cdef int tot_n = seq_n
        cdef uint8_t *bseq
        cdef abpoa_res_t res
        cdef uint8_t **cons_seq
        cdef uint8_t **msa_seq
        cdef int *cons_len
        cdef int cons_n=0
        cdef int msa_l=0

        _cons_len, _cons_seq, _msa_seq = [], [], []

        if out_cons: self.abpt.out_cons = 1
        else: self.abpt.out_cons = 0
        if out_msa: self.abpt.out_msa = 1
        else: self.abpt.out_msa = 0
        if out_pog: 
            if isinstance(out_pog, str): out_pog = bytes(out_pog, 'utf-8')
            self.abpt.out_pog = out_pog
        else: self.abpt.out_pog = NULL

        abpoa_post_set_para(&self.abpt)
        abpoa_reset_graph(self.ab, &self.abpt, len(seqs[0]))
        if incr_fn:
            if isinstance(incr_fn, str):
                incr_fn = bytes(incr_fn, 'utf-8')
            self.abpt.incr_fn = incr_fn
            abpoa_restore_graph(self.ab, &self.abpt)
            exist_n = self.ab[0].abs[0].n_seq
            tot_n += exist_n
        else:
            self.abpt.incr_fn = NULL

        self.ab[0].abs[0].n_seq += seq_n

        if tot_n < 2: return msa_result(tot_n, cons_n, _cons_len, _cons_seq, msa_l, _msa_seq)

        for read_i, seq in enumerate(seqs):
            seq_l = len(seq)
            bseq = <uint8_t*>malloc(seq_l * cython.sizeof(uint8_t))
            for i in range(seq_l):
                bseq[i] = self.seq2int_dict[seq[i]]
            res.n_cigar = 0
            abpoa_align_sequence_to_graph(self.ab, &self.abpt, bseq, seq_l, &res)

            abpoa_add_graph_alignment(self.ab, &self.abpt, bseq, seq_l, NULL, res, exist_n+read_i, tot_n, 1)
            free(bseq)
            if res.n_cigar: free(res.graph_cigar)

        if self.abpt.out_cons:
            abpoa_generate_consensus(self.ab, &self.abpt, NULL, &cons_seq, NULL, &cons_len, &cons_n)
            for i in range(cons_n):
                _cons_len.append(cons_len[i])
                cons_seq1 = ''
                for c in cons_seq[i][:cons_len[i]]:
                    if isinstance(c, bytes): c = ord(c)
                    cons_seq1 += self.int2seq_dict[c]
                _cons_seq.append(cons_seq1)
            if cons_n > 0:
                for i in range(cons_n):
                    free(cons_seq[i])
                free(cons_seq)
                free(cons_len)
        if self.abpt.out_msa:
            abpoa_generate_rc_msa(self.ab, &self.abpt, NULL, &msa_seq, &msa_l)
            for i in range(tot_n):
                msa_seq1 = ''
                for c in msa_seq[i][:msa_l]:
                    if isinstance(c, bytes): c = ord(c)
                    msa_seq1 += self.int2seq_dict[c]
                _msa_seq.append(msa_seq1)
            if msa_l > 0:
                for i in range(tot_n):
                    free(msa_seq[i])
                free(msa_seq)
        if self.abpt.out_pog:
            abpoa_dump_pog(self.ab, &self.abpt)
        return msa_result(tot_n, int(cons_n), _cons_len, _cons_seq, int(msa_l), _msa_seq)
