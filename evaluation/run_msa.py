import os, re, sys
import Levenshtein as lv
import mappy as mp

sim_read_cnt = 100  # 2000 * 5 * 5 = 50000 in total
sim_read_len = [300, 1000, 2000, 3000]
sim_read_depth = [3, 10, 20]
sim_names = ['nanosim', 'pbsim']
msa_out_dir = '/home/yan/data/abpoa_simu/low_depth'

msa_abPOA = '/home/yan/program/abPOA/msa_abPOA'
msa_spoa = '/home/yan/software/spoa/msa_spoa'

# spoa-like 5 4 8/10 6/4
match = 2
mismatch = 4
gap_open1 = 4
gap_ext1 = 2
gap_oe1 = 6
gap_open2 = 24
gap_ext2 = 1
gap_oe2 = 25

# 2. generate mulutiple sequence
eval_out = msa_out_dir + '/new.eval.out'
cons_abPOA = msa_out_dir + '/abPOA.cons'
abpoa_cons_algr = 0
cons_spoa = msa_out_dir + '/spoa.cons'
# val_out = msa_out_dir + '/val.out'

def cal_err_rate(cons_fa, cons_res, length, out_fn):
    with open(cons_fa) as cons_fp, open(cons_res) as res_fp:
        err = []
        res_len = []
        for s1, s2 in zip(cons_fp, res_fp):
            if s1.startswith('>'): continue
            a = mp.Aligner(seq = s1[:-1])
            hit = 0
            for h in a.map(s2[:-1]):
                if h.is_primary:
                    # print(h.NM)
                    err.append(h.NM)
                    res_len.append(len(s2[:-1]))
                    hit = 1
                    break
            if hit == 0:
                print('>S1\n{}>S2\n{}'.format(s1, s2))
            else:
                if h.NM / length > 0.2:
                    err.pop()
                    print(h.NM, s1)
        err_cnt = sum(err)/len(err)
        ave_len = sum(res_len)/len(err)
        # err_rate = 100 * sum(err)/len(err)/ length
        err_rate = 100 * sum(err)/len(err)/ ave_len
        cnt = len(err)
        # print(err_cnt, err_rate, ave_len)
        with open(eval_out, 'a') as out_fp:
            out_fp.write('{:.2f}\t'.format(err_rate))
        # os.system('echo \"\terr: {:.2f} ({:.2f}%), ave_len: {:.2f}, cnt: {}\" >> {}'.format(err_cnt, err_rate, ave_len, cnt, eval_out))
        return sum(err)/len(err), len(err)

def cal_raw_err_rate(cons_fa, seq_fa, n):
    with open(cons_fa) as cons_fp, open(seq_fa) as seq_fp:
        err = []
        res_len = []
        for cons in cons_fp:
            if cons.startswith('>'): continue
            a = mp.Aligner(seq = cons[:-1])
            i = 0

            for seq in seq_fp:
                if seq.startswith('>'):continue
                for h in a.map(seq[:-1]):
                    if h.is_primary:
                        err.append(h.NM)
                        res_len.append(len(seq[:-1]))
                        break
                i += 1
                if i == n: break
        err_rate = 100 * sum(err)/len(err)/ (sum(res_len) / len(res_len))
        return err_rate



if __name__ == '__main__':
    os.system('rm {} -f'.format(eval_out))
    for length in sim_read_len:
        for depth in sim_read_depth:
            cons = msa_out_dir + '/len{}_depth{}_cnt{}.cons.fa'.format(length, depth, sim_read_cnt)
            # run evaluation
            band_width = [10, -1] 
            for name in sim_names:
                seq = msa_out_dir + '/len{}_depth{}_cnt{}.{}.fa'.format(length, depth, sim_read_cnt, name)
                raw_rate = cal_raw_err_rate(cons, seq, depth)

                with open(eval_out, 'a') as out_fp:
                    out_fp.write('{}\t{}\t{}\t{:.2f}\t'.format(name, length, depth, raw_rate))
                
                os.system('rm -f {}'.format(cons_spoa))
                print('{} -l1 -n {} -s {} -m {} -x -{} -o -{} -e -{} -q -{} -c -{} > {} 2>> {}'.format(msa_spoa, depth, seq, match, mismatch, gap_oe1, gap_ext1, gap_oe2, gap_ext2, cons_spoa, eval_out))
                os.system('{} -l1 -n {} -s {} -m {} -x -{} -o -{} -e -{} -q -{} -c -{} > {} 2>> {}'.format(msa_spoa, depth, seq, match, mismatch, gap_oe1, gap_ext1, gap_oe2, gap_ext2, cons_spoa, eval_out))
                cal_err_rate(cons,  cons_spoa, length, eval_out) 

                for w in band_width:
                    os.system('rm -f {}'.format(cons_abPOA))
                    print('{} -n {} -s {} -C {} -w {} -m {} -x {} -o{},{} -e{},{} > {} 2>> {}'.format(msa_abPOA, depth, seq, abpoa_cons_algr, w, match, mismatch, gap_open1, gap_open2, gap_ext1, gap_ext2, cons_abPOA, eval_out))
                    os.system('{} -n {} -s {} -C {} -w {} -m {} -x {} -o{},{} -e{},{} > {} 2>> {}'.format(msa_abPOA, depth, seq, abpoa_cons_algr, w, match, mismatch, gap_open1, gap_open2, gap_ext1, gap_ext2, cons_abPOA, eval_out))
                    cal_err_rate(cons,  cons_abPOA, length, eval_out) 
                with open(eval_out, 'a') as out_fp:
                    out_fp.write('\n')
