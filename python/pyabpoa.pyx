import re, sys, os
from libc.stdlib cimport malloc, free
from libc.stdint cimport uint8_t
from collections import defaultdict as dd
cimport cython
from cabpoa cimport *
from libc.string cimport strcpy


cdef class msa_result:
    cdef int n_seq
    cdef int n_cons
    cdef clu_n_seq, clu_read_ids, cons_len, cons_seq, cons_cov, cons_qv  # _cons_len:[int], _cons_seq:['']
    cdef int msa_len
    cdef msa_seq  # _msa_seq:['']

    def __cinit__(self, n_seq, n_cons, clu_n_seq, clu_read_ids, cons_len, cons_seq, cons_cov, cons_qv, msa_len, msa_seq):
        self.n_seq = n_seq
        self.n_cons = n_cons
        self.clu_n_seq = clu_n_seq
        self.clu_read_ids = clu_read_ids
        self.cons_len = cons_len
        self.cons_seq = cons_seq
        self.cons_cov = cons_cov
        self.cons_qv = cons_qv
        self.msa_len = msa_len
        self.msa_seq = msa_seq

    @property
    def n_seq(self): return self.n_seq

    @property
    def n_cons(self): return self.n_cons

    @property
    def clu_n_seq(self): return self.clu_n_seq

    @property
    def clu_read_ids(self): return self.clu_read_ids

    @property
    def cons_len(self): return self.cons_len

    @property
    def cons_seq(self): return self.cons_seq

    @property
    def cons_cov(self): return self.cons_cov

    @property
    def cons_qv(self): return self.cons_qv

    @property
    def msa_len(self): return self.msa_len

    @property
    def msa_seq(self): return self.msa_seq

    def print_msa(self): 
        if not self.msa_seq: return
        for i, s in enumerate(self.msa_seq):
            if i < self.n_seq:
                print(f'>Seq_{i+1}')
            else:
                if self.n_cons > 1:
                    cons_id = f'_{i-self.n_seq+1} {",".join(list(map(str, self.clu_read_ids[i-self.n_seq])))}'
                else:
                    cons_id = ''
                print(f'>Consensus_sequence{cons_id}')
            print(s)
        return


def set_seq_int_dict(m):
    if m == 5:  # ACGTN ==> 01234, U ==> 4
        seqs = 'ACGUTN'
        ints = [0, 1, 2, 3, 3, 4]
    elif m == 27:  # ACGTN    ==> 01234, BDEFH... ==> 56789...
        seqs = 'ACGTNBDEFHIJKLMOPQRSUVWXYZ*'
        ints = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
    else:
        raise Exception(f'Unexpected m: {m}')

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
    cdef char* _score_mat_fn
    cdef bytes _score_mat_fn_b

    def __cinit__(self, aln_mode='g', is_aa=False,
                  match=2, mismatch=4, score_matrix=b'', gap_open1=4, gap_open2=24, gap_ext1=2, gap_ext2=1,
                  extra_b=10, extra_f=0.01,
                  cons_algrm='HB'):
        self.ab = abpoa_init()

        if aln_mode == 'g':
            self.abpt.align_mode = ABPOA_GLOBAL_MODE
        elif aln_mode == 'l':
            self.abpt.align_mode = ABPOA_LOCAL_MODE
        elif aln_mode == 'e':
            self.abpt.align_mode = ABPOA_EXTEND_MODE
        else:
            raise Exception(f'Unknown align mode: {aln_mode}')
        if is_aa:
            self.abpt.m = 27
            self.abpt.mat = <int*>malloc(27 * 27 * cython.sizeof(int))
        else:
            self.abpt.m = 5
            self.abpt.mat = <int*>malloc(25 * cython.sizeof(int))
        self.abpt.match = match
        self.abpt.mismatch = mismatch

        if score_matrix:
            self.abpt.use_score_matrix = 1
            _score_mat_fn_b = score_matrix.encode('utf-8')
            _score_mat_fn = <char*> malloc(len(_score_mat_fn_b) + 1)
            strcpy(_score_mat_fn, _score_mat_fn_b)
            self.abpt.mat_fn = _score_mat_fn

        self.abpt.gap_open1 = gap_open1
        self.abpt.gap_open2 = gap_open2
        self.abpt.gap_ext1 = gap_ext1
        self.abpt.gap_ext2 = gap_ext2
        self.abpt.ret_cigar = 1

        self.abpt.wb = extra_b
        self.abpt.wf = extra_f 
        self.abpt.use_qv = 0
        self.abpt.end_bonus = -1 # disable end_bonus/zdrop
        self.abpt.zdrop = -1
        self.abpt.disable_seeding = 1
        self.abpt.progressive_poa = 0
        if cons_algrm.upper() == 'MF':
            self.abpt.cons_algrm = ABPOA_MF
        elif cons_algrm.upper() == 'HB':
            self.abpt.cons_algrm = ABPOA_HB
        else:
            raise Exception(f'Unknown conseneus calling mode: {cons_algrm}')

        self.seq2int_dict, self.int2seq_dict = set_seq_int_dict(self.abpt.m)

    def __dealloc__(self):
        free(self.abpt.mat)
        abpoa_free(self.ab)

    def __bool__(self):
        return self.ab != NULL

    cdef int* _prepare_weights(self, seq, qscore):
        cdef Py_ssize_t i, seq_l = len(seq)
        cdef int *weights = NULL

        weights = <int*>malloc(seq_l * cython.sizeof(int))
        if weights == NULL:
            raise MemoryError()

        try:
            if len(qscore) != seq_l:
                raise ValueError('Each qscore array must have the same length as its sequence.')
            for i in range(seq_l):
                weights[i] = int(qscore[i])
                if weights[i] < 0:
                    raise ValueError('Qscores must be non-negative integers.')
        except:
            free(weights)
            raise
        return weights

    cdef void _add_sequences(self, seqs, qscores, int exist_n, int tot_n):
        cdef uint8_t *bseq
        cdef int *weights = NULL
        cdef int seq_l
        cdef Py_ssize_t i
        cdef abpoa_res_t res

        if qscores is not None and len(qscores) != len(seqs):
            raise ValueError('qscores must contain one entry per input sequence.')

        for read_i, seq in enumerate(seqs):
            seq_l = len(seq)
            bseq = <uint8_t*>malloc(seq_l * cython.sizeof(uint8_t))
            if bseq == NULL:
                raise MemoryError()
            weights = NULL
            res.n_cigar = 0
            res.graph_cigar = NULL
            try:
                for i in range(seq_l):
                    bseq[i] = self.seq2int_dict[seq[i]]
                if qscores is not None:
                    weights = self._prepare_weights(seq, qscores[read_i])
                abpoa_align_sequence_to_graph(self.ab, &self.abpt, bseq, seq_l, &res)
                abpoa_add_graph_alignment(self.ab, &self.abpt, bseq, weights, seq_l, NULL, res, exist_n+read_i, tot_n, 1)
            finally:
                free(bseq)
                if weights != NULL:
                    free(weights)
                if res.n_cigar:
                    free(res.graph_cigar)

    def msa(self, seqs, out_cons, out_msa, max_n_cons=1, min_freq=0.25, out_pog=b'', incr_fn=b'', qscores=None):
        cdef int seq_n = len(seqs)
        cdef int exist_n = 0
        cdef int tot_n = seq_n
        cdef abpoa_cons_t abc

        if out_cons: self.abpt.out_cons = 1
        else: self.abpt.out_cons = 0
        if out_msa: self.abpt.out_msa = 1
        else: self.abpt.out_msa = 0
        if max_n_cons < 1 or max_n_cons > 2:
            raise Exception('Error: max number of consensus sequences should be 1 or 2.')
        self.abpt.max_n_cons = max_n_cons
        self.abpt.min_freq = min_freq
        self.abpt.use_qv = 1 if qscores is not None else 0
        if out_pog: 
            if isinstance(out_pog, str): out_pog = bytes(out_pog, 'utf-8')
            self.abpt.out_pog = out_pog
        else: self.abpt.out_pog = NULL

        abpoa_post_set_para(&self.abpt)
        abpoa_reset(self.ab, &self.abpt, len(seqs[0]))
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
        self._add_sequences(seqs, qscores, exist_n, tot_n)

        if self.abpt.out_msa:
            abpoa_generate_rc_msa(self.ab, &self.abpt)
        elif self.abpt.out_cons:
            abpoa_generate_consensus(self.ab, &self.abpt)
        abc = self.ab[0].abc[0]

        n_cons, clu_n_seq, clu_read_ids, cons_len, cons_seq, cons_cov, cons_qv, msa_len, msa_seq = 0, [], [], [], [], [], [], 0, []
        n_cons = abc.n_cons
        for i in range(n_cons):
            clu_n_seq.append(abc.clu_n_seq[i])
            cons_len.append(abc.cons_len[i])
            clu_read_ids1, cons_seq1, cons_cov1, cons_qv1 = [], '', [], ''
            for j in range(abc.clu_n_seq[i]):
                clu_read_ids1.append(abc.clu_read_ids[i][j])
            clu_read_ids.append(clu_read_ids1)
            for j in range(abc.cons_len[i]):
                c = abc.cons_base[i][j]
                if isinstance(c, bytes): c = ord(c)
                cons_seq1 += self.int2seq_dict[c]
                cons_cov1.append(abc.cons_cov[i][j])
                if abc.cons_phred_score != NULL:
                    cons_qv1 += chr(abc.cons_phred_score[i][j])
            cons_seq.append(cons_seq1)
            cons_cov.append(cons_cov1)
            cons_qv.append(cons_qv1)

        msa_len = abc.msa_len
        if msa_len > 0:
            for i in range(abc.n_seq + n_cons):
                msa_seq1 = ''
                for c in abc.msa_base[i][:msa_len]:
                    if isinstance(c, bytes): c = ord(c)
                    msa_seq1 += self.int2seq_dict[c]
                msa_seq.append(msa_seq1)

        if self.abpt.out_pog:
            abpoa_dump_pog(self.ab, &self.abpt)
        return msa_result(tot_n, n_cons, clu_n_seq, clu_read_ids, cons_len, cons_seq, cons_cov, cons_qv, msa_len, msa_seq)

    def msa_align(self, seqs, out_cons, out_msa, max_n_cons=1, min_freq=0.25, incr_fn=b'', qscores=None):
        cdef int seq_n = len(seqs)
        cdef int exist_n = 0
        cdef int tot_n = seq_n

        if out_cons: self.abpt.out_cons = 1
        else: self.abpt.out_cons = 0
        if out_msa: self.abpt.out_msa = 1
        else: self.abpt.out_msa = 0
        if max_n_cons < 1 or max_n_cons > 2:
            raise Exception('Error: max number of consensus sequences should be 1 or 2.')
        self.abpt.max_n_cons = max_n_cons
        self.abpt.min_freq = min_freq
        self.abpt.use_qv = 1 if qscores is not None else 0

        abpoa_post_set_para(&self.abpt)
        abpoa_reset(self.ab, &self.abpt, len(seqs[0]))
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
        self._add_sequences(seqs, qscores, exist_n, tot_n)
        return self

    def msa_add(self, new_seqs, qscores=None):
        if isinstance(new_seqs, str):
            raise TypeError('Expected a list of strings. If you want to add a single sequence, pass it as a list: ["ACGT..."]')

        cdef int exist_n = 0
        exist_n = self.ab[0].abs[0].n_seq
        if (exist_n == 0):
            raise Exception('Error: no existing sequences in the graph. Please run msa() or msa_align() first.')
        cdef int seq_n = len(new_seqs)
        cdef int tot_n = seq_n

        tot_n += exist_n
        self.ab[0].abs[0].n_seq += seq_n
        if qscores is not None:
            self.abpt.use_qv = 1
        self._add_sequences(new_seqs, qscores, exist_n, tot_n)
        return self

    def msa_output(self):
        cdef abpoa_cons_t abc
        if self.abpt.out_msa:
            abpoa_generate_rc_msa(self.ab, &self.abpt)
        elif self.abpt.out_cons:
            abpoa_generate_consensus(self.ab, &self.abpt)
        seq_n = self.ab[0].abs[0].n_seq
        abc = self.ab[0].abc[0]

        n_cons, clu_n_seq, clu_read_ids, cons_len, cons_seq, cons_cov, cons_qv, msa_len, msa_seq = 0, [], [], [], [], [], [], 0, []
        n_cons = abc.n_cons
        for i in range(n_cons):
            clu_n_seq.append(abc.clu_n_seq[i])
            cons_len.append(abc.cons_len[i])
            clu_read_ids1, cons_seq1, cons_cov1, cons_qv1 = [], '', [], ''
            for j in range(abc.clu_n_seq[i]):
                clu_read_ids1.append(abc.clu_read_ids[i][j])
            clu_read_ids.append(clu_read_ids1)
            for j in range(abc.cons_len[i]):
                c = abc.cons_base[i][j]
                if isinstance(c, bytes): c = ord(c)
                cons_seq1 += self.int2seq_dict[c]
                cons_cov1.append(abc.cons_cov[i][j])
                if abc.cons_phred_score != NULL:
                    cons_qv1 += chr(abc.cons_phred_score[i][j])
            cons_seq.append(cons_seq1)
            cons_cov.append(cons_cov1)
            cons_qv.append(cons_qv1)

        msa_len = abc.msa_len
        if msa_len > 0:
            for i in range(abc.n_seq + n_cons):
                msa_seq1 = ''
                for c in abc.msa_base[i][:msa_len]:
                    if isinstance(c, bytes): c = ord(c)
                    msa_seq1 += self.int2seq_dict[c]
                msa_seq.append(msa_seq1)
        return msa_result(seq_n, n_cons, clu_n_seq, clu_read_ids, cons_len, cons_seq, cons_cov, cons_qv, msa_len, msa_seq)
