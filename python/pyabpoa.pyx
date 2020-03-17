import re, sys
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


cdef class msa_aligner:
    cdef abpoa_t *ab
    cdef abpoa_para_t abpt
    cdef seq_nt4_dict, nt4_seq_dict

    def __cinit__(self, aln_mode='g', match=2, mismatch=4, gap_open1=4, gap_open2=24, gap_ext1=2, gap_ext2=1,
            band_width=-1, end_bonus=-1,zdrop=-1, ystop_iter_n=-1, ystop_min_processed_n=YSTOP_MIN_NUM, ystop_min_frac=YSTOP_MIN_FRAC,
            cons_agrm=ABPOA_HB, multip=1, min_freq=0.3):
        self.ab = abpoa_init()

        if aln_mode == 'g':
            self.abpt.align_mode = ABPOA_GLOBAL_MODE
        elif aln_mode == 'e':
            self.abpt.align_mode = ABPOA_EXTEND_MODE
        else:
            print('Unknown align mode: {}'.format(aln_mode))
            sys.exit(1)
        self.abpt.match = match
        self.abpt.mismatch = mismatch
        self.abpt.gap_open1 = gap_open1
        self.abpt.gap_open2 = gap_open2
        self.abpt.gap_ext1 = gap_ext1
        self.abpt.gap_ext2 = gap_ext2
        self.abpt.ret_cigar = 1

        self.abpt.m = 5
        self.abpt.mat = <int*>malloc(25 * cython.sizeof(int))

        self.abpt.bw = band_width
        self.abpt.end_bonus = end_bonus
        self.abpt.zdrop = zdrop
        self.abpt.ystop_iter_n = ystop_iter_n
        self.abpt.ystop_min_processed_n = ystop_min_processed_n
        self.abpt.ystop_min_frac = ystop_min_frac

        self.abpt.cons_agrm = cons_agrm
        self.abpt.multip = multip
        self.abpt.min_freq = min_freq

        self.abpt.simd_flag = simd_check()


        self.seq_nt4_dict = dd(lambda:4)
        self.seq_nt4_dict['A'] = 0
        self.seq_nt4_dict['a'] = 0
        self.seq_nt4_dict['C'] = 1
        self.seq_nt4_dict['c'] = 1
        self.seq_nt4_dict['G'] = 2
        self.seq_nt4_dict['g'] = 2
        self.seq_nt4_dict['T'] = 3
        self.seq_nt4_dict['t'] = 3
        self.nt4_seq_dict = dd(lambda:'-')
        self.nt4_seq_dict[0] = 'A'
        self.nt4_seq_dict[1] = 'C'
        self.nt4_seq_dict[2] = 'G'
        self.nt4_seq_dict[3] = 'T'

    def __dealloc__(self):
        free(self.abpt.mat)
        abpoa_free(self.ab, &self.abpt)

    def __bool__(self):
        return self.ab != NULL


    def msa(self, seqs, out_cons, out_msa, out_pog=None):
        cdef int seq_n = len(seqs)
        cdef int tot_n = seq_n 
        cdef uint8_t *bseq
        cdef abpoa_res_t res
        cdef int check_ystop
        cdef int ystop
        cdef int ystop_n = 0
        cdef uint8_t **cons_seq
        cdef uint8_t **msa_seq
        cdef int *cons_len
        cdef int cons_n=0
        cdef int msa_l=0

        _cons_len, _cons_seq, _msa_seq = [], [], []

        if seq_n < 2: return msa_result(seq_n, cons_n, _cons_len, _cons_seq, msa_l, _msa_seq)

        if out_cons: self.abpt.out_cons = 1
        if out_msa: self.abpt.out_msa = 1
        abpoa_post_set_para(&self.abpt)

        # abpoa_reset_graph(self.ab, seq_n, &self.abpt)
        for read_i, seq in enumerate(seqs):
            seq_l = len(seq)
            bseq = <uint8_t*>malloc(seq_l * cython.sizeof(uint8_t))
            for i in range(seq_l):
                bseq[i] = self.seq_nt4_dict[seq[i]]
            res.n_cigar = 0
            abpoa_align_sequence_to_graph(self.ab, &self.abpt, bseq, seq_l, &res)
            if self.abpt.ystop_iter_n > 0 and self.abpt.out_cons and not self.abpt.out_msa and read_i>= self.abpt.ystop_min_processed_n:
                self.abpt.ystop_min_wei = int((read_i + 1) * self.abpt.ystop_min_frac)
                check_ystop = 1
            else: check_ystop = 0

            ystop = abpoa_add_graph_alignment(self.ab, &self.abpt, bseq, seq_l, res.n_cigar, res.graph_cigar, check_ystop, read_i, seq_n)
            free(bseq)

            if res.n_cigar: free(res.graph_cigar)
            if check_ystop:
                # fprintf(stderr, "%c", "NS"[ystop]); \
                if ystop: ystop_n += 1
                else: ystop_n = 0
                if ystop_n >= self.abpt.ystop_iter_n:
                    tot_n = read_i+1
                    break
        if self.abpt.out_cons:
            abpoa_generate_consensus(self.ab, self.abpt.cons_agrm, self.abpt.multip, self.abpt.min_freq, tot_n, NULL, &cons_seq, &cons_len, &cons_n)
            for i in range(cons_n):
                _cons_len.append(cons_len[i])
                cons_seq1 = ''
                for c in cons_seq[i][:cons_len[i]]:
                    if isinstance(c, bytes): c = ord(c)
                    cons_seq1 += self.nt4_seq_dict[c]
                _cons_seq.append(cons_seq1)
            if cons_n > 0:
                for i in range(cons_n):
                    free(cons_seq[i])
                free(cons_seq)
                free(cons_len)
        if self.abpt.out_msa:
            abpoa_generate_rc_msa(self.ab, seq_n, NULL, &msa_seq, &msa_l)
            for i in range(seq_n):
                msa_seq1 = ''
                for c in msa_seq[i][:msa_l]:
                    if isinstance(c, bytes): c = ord(c)
                    msa_seq1 += self.nt4_seq_dict[c]
                _msa_seq.append(msa_seq1)
            if msa_l > 0:
                for i in range(seq_n):
                    free(msa_seq[i])
                free(msa_seq)
        return msa_result(seq_n, int(cons_n), _cons_len, _cons_seq, int(msa_l), _msa_seq)
        

cdef class s2g_result:
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

cdef class s2g_aligner:
    cdef abpoa_para_t abpt
    cdef abpoa_t *ab
    cdef seq_nt4_dict, nt4_seq_dict

    def __cinit__(self, gfa_fn, aln_mode='g', match=2, mismatch=4, gap_open1=4, gap_open2=24, gap_ext1=2, gap_ext2=1,
            band_width=-1, end_bonus=-1, zdrop=-1):
        if gfa_fn:
            self.ab = abpoa_init()
            print(gfa_fn) # build graph from GFA
        if aln_mode == 'g':
            self.abpt.align_mode = ABPOA_GLOBAL_MODE
        elif aln_mode == 'e':
            self.abpt.align_mode = ABPOA_EXTEND_MODE
        else:
            print('Unknown align mode: {}'.format(aln_mode))
            sys.exit(1)
        self.abpt.match = match
        self.abpt.mismatch = mismatch
        self.abpt.gap_open1 = gap_open1
        self.abpt.gap_open2 = gap_open2
        self.abpt.gap_ext1 = gap_ext1
        self.abpt.gap_ext2 = gap_ext2
        self.abpt.ret_cigar = 1

        self.abpt.m = 5
        self.abpt.mat = <int*>malloc(25 * cython.sizeof(int))

        self.abpt.bw = band_width
        self.abpt.end_bonus = end_bonus
        self.abpt.zdrop = zdrop
        self.abpt.ystop_iter_n = -1
        # self.abpt.ystop_min_processed_n = ystop_min_processed_n
        # self.abpt.ystop_min_frac = ystop_min_frac

        # self.abpt.cons_agrm = cons_agrm
        self.abpt.multip = 1 #multip
        # self.abpt.min_freq = min_freq

        self.abpt.simd_flag = simd_check()

        self.seq_nt4_dict = dd(lambda:4)
        self.seq_nt4_dict['A'] = 0
        self.seq_nt4_dict['a'] = 0
        self.seq_nt4_dict['C'] = 1
        self.seq_nt4_dict['c'] = 1
        self.seq_nt4_dict['G'] = 2
        self.seq_nt4_dict['g'] = 2
        self.seq_nt4_dict['T'] = 3
        self.seq_nt4_dict['t'] = 3
        self.nt4_seq_dict = dd(lambda:'-')
        self.nt4_seq_dict[0] = 'A'
        self.nt4_seq_dict[1] = 'C'
        self.nt4_seq_dict[2] = 'G'
        self.nt4_seq_dict[3] = 'T'

    def __dealloc__(self):
        free(self.abpt.mat)
        abpoa_free(self.ab, &self.abpt)

    def __bool__(self):
        return self.ab != NULL

    def print_graph(self):
        return

    def dot_plot(self):
        return

    def s2g(self, seq):
        return

