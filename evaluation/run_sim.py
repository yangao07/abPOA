# input ref.fa
# output num_len_err.fa, num_len_err.cons
import sys, os
import random
from datetime import datetime
import utils as ut

threads = 4
fxtools = 'fxtools'
pbsim = 'pbsim.sh'
pbsim_model = 'CLR'
pbsim_min_acc = 0.85
pbsim_max_acc = 1.0
pbsim_mean_acc = 0.88
pbsim_seed = 0

simlord = 'simlord'
max_pass_pb = 5
nanosim = 'simulator.py'
nano_basecaller = 'guppy'
train_model = '/home/yan/software/NanoSim/models/human_NA12878_DNA_FAB49712_guppy/training'


sim_read_cnt = 100  # 2000 * 5 * 5 = 50000 in total
sim_read_len = [300, 1000, 2000, 3000]
sim_read_depth = [3, 10, 20]
nanosim_depth = {3:3, 5:6, 10:11, 15:17, 20:23, 30:34}

# ref_fa = '/home/gaoy1/data/genome/hg38/chr1.fa' 
ref_fa = '/home/yan/data/chr16.fa'
#ref_amb = '/home/gaoy1/data/genome/hg38/chr1.fa.amb' 
ref_amb = '/home/yan/data/chr16.fa.amb'
# chrom = 'chr1' 
chrom = 'chr16'
# chrom_len = 248956422
chrom_len = 90354753
out_dir = '/home/yan/data/abpoa_simu/pbsim'


# amb: amb_start:amb_len
def get_amb(in_amb, len_max):
    first = True
    amb = dict()
    with open(in_amb) as amb_fp:
        for line in amb_fp:
            if first:
                first = False
                continue
            amb_start = int(line.rsplit()[0]) + 1
            amb_len = int(line.rsplit()[1])
            amb[amb_start] = amb_len
            if amb_start >= len_max:
                return amb
    return amb


def ovlp_len(s1, l1, s2, l2):
    e1 = s1 + l1 - 1
    e2 = s2 + l2 - 1
    if e1 < s2 or e2 < s1:
        return 0
    else:
        return 1


# randomly select one position in ref
def get_random_start(chrom_len, len, amb):
    start = -1
    while start < 0:
        start = random.randint(1, chrom_len - len + 1)
        for s in amb:
            if ovlp_len(s, amb[s], start, len):
                start = -1
                break
    return start


def get_seq(fxtools, chrom, start, len, out_fa):
    ut.exec_cmd(sys.stderr, 'get_seq',
                '{} sd {} {} {} {} > {}'.format(fxtools, ref_fa, chrom, start, start + len - 1, out_fa))


def run_simlord(ref_fa, depth, length, tmp_out, sim_read_fa):
    ut.exec_cmd(sys.stderr, 'SimLoRD', '{} -rr {} -n {} -mp {} -fl {} --no-sam {}'.format(simlord, ref_fa, depth, max_pass_pb, length, tmp_out))
    ut.exec_cmd(sys.stderr, 'merge reads', '{} qa {}.fastq >> {}'.format(fxtools, tmp_out, sim_read_fa))
    ut.exec_cmd(sys.stderr, 'rm tmp', 'rm {}.fastq'.format(tmp_out))

def run_pbsim(ref_fa, depth, length, tmp_out, sim_read_fa):
    ut.exec_cmd(sys.stderr, 'pbsim', '{} {} {} {} {} {} {} {} {} {} {}'.format(pbsim, ref_fa, tmp_out, pbsim_model, depth, length, length, pbsim_min_acc, pbsim_max_acc, pbsim_mean_acc, pbsim_seed))
    ut.exec_cmd(sys.stderr, 'merge reads', 'cat {}*.fastq | {} qa - | head -n {} >> {}'.format(tmp_out, fxtools, depth*2, sim_read_fa))
    ut.exec_cmd(sys.stderr, 'rm tmp', 'rm {}*.ref {}*.maf {}*.fastq'.format(tmp_out, tmp_out, tmp_out))
    return

def run_nanosim(ref_fa, depth, length, tmp_out, sim_read_fa):
    depth = nanosim_depth[depth]
    ut.exec_cmd(sys.stderr, 'NanoSim', '{} genome -rg {} -c {} -t {} -n {} -max {} -min {} -b {} -s 0 -o {}'.format(nanosim, ref_fa, train_model, threads, depth, length, int(length*0.95), nano_basecaller, tmp_out))

    ut.exec_cmd(sys.stderr, 'merge reads', 'cat {}_aligned_reads.fasta >> {}'.format(tmp_out, sim_read_fa))
    ut.exec_cmd(sys.stderr, 'rm tmp', 'rm {}_aligned_reads.fasta {}_error*fasta'.format(tmp_out, tmp_out))


if __name__ == '__main__':
    # if len(sys.argv) != 3:
    #     print "Usage:"
    #     print "{} in_fa per_cnt ".format(sys.argv[0])
    #     sys.exit(1)

    random.seed(datetime.now())

    amb = get_amb(ref_amb, chrom_len)

    # simulate reads
    for length in sim_read_len:
        for depth in sim_read_depth:
            out_pre = out_dir + '/len{}_depth{}_cnt{}'.format(length, depth, sim_read_cnt)

            cons_fa = out_pre + '.cons.fa'
            simlord_read_fa = out_pre + '.simlord.fa'
            pbsim_read_fa = out_pre + '.pbsim.fa'
            nanosim_read_fa = out_pre + '.nanosim.fa'
            tmp_ref = out_pre + '.tmp.ref.fa'
            tmp_out = out_pre + '.tmp'
            os.system('rm -rf {} {} {}'.format(pbsim_read_fa, nanosim_read_fa, cons_fa))
            for i in range(sim_read_cnt):
                start = get_random_start(chrom_len, length, amb)
                get_seq(fxtools, chrom, start, length, tmp_ref)

                # run_simlord(tmp_ref, depth, length, tmp_out, simlord_read_fa)
                run_pbsim(tmp_ref, depth, length, tmp_out, pbsim_read_fa)
                # run_nanosim(tmp_ref, depth, length, tmp_out, nanosim_read_fa)

                ut.exec_cmd(sys.stderr, 'merge cons', 'cat {} | {} qa - >> {}; rm {}'.format(tmp_ref, fxtools, cons_fa, tmp_ref))
