#!/usr/bin/env python3

from qtlseq.params import Params


pm = Params('qtlseq')
args = pm.set_options()


import os
import sys
import glob
import subprocess as sbp
from qtlseq.refindex import RefIndex
from qtlseq.trim import Trim
from qtlseq.alignment import Alignment
from qtlseq.bamfilt import BamFilt
from qtlseq.mpileup import Mpileup
from qtlseq.vcf2index import Vcf2Index
from qtlseq.utils import time_stamp, clean_cmd, call_log


class QTLseq(object):

    def __init__(self, args):
        args = pm.check_max_threads(args)
        self.N_fastq = pm.check_args(args)
        self.args = args
        self.out = args.out
        self.parent_fastq, self.parent_bam = self.get_files(args.parent)
        self.bulk1_fastq, self.bulk1_bam = self.get_files(args.bulk1)
        self.bulk2_fastq, self.bulk2_bam = self.get_files(args.bulk2)

        self.mkdir()
        self.link_ref()
        self.link_bam('parent', self.parent_bam)
        self.link_bam('bulk1', self.bulk1_bam)
        self.link_bam('bulk2', self.bulk2_bam)

    def get_files(self, input_names):
        fastq = []
        bam = []
        for input_name in input_names:
            if ',' in input_name:
                fastq.append(input_name)
            else:
                bam.append(input_name)
        return fastq, bam

    def mkdir(self):
        os.mkdir('{}'.format(self.out))
        os.mkdir('{}/log'.format(self.out))
        if self.N_fastq > 0 and self.args.trim:
            os.mkdir('{}/00_fastq'.format(self.out))
        os.mkdir('{}/10_ref'.format(self.out))
        os.mkdir('{}/20_bam'.format(self.out))

    def link_ref(self):
        path_to_ref = os.path.abspath(self.args.ref)
        ref = os.path.basename(self.args.ref)
        sym_ref = '{}/10_ref/{}'.format(self.out, ref)
        os.symlink(path_to_ref, sym_ref)
        self.args.ref = sym_ref

    def link_bam(self, label, input_files):
        for i, input_name in enumerate(input_files):
            path_to_bam = os.path.abspath(input_name)
            index = '{}.1{:0>3}'.format(label, i)
            os.symlink(path_to_bam,
                       '{}/20_bam/{}.bam'.format(self.out, index))

    def refindex(self):
        ri = RefIndex(args)
        ri.run()

    def trimming(self):
        tm = Trim(args)
        for i, fastq in enumerate(self.parent_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'parent.0{:0>3}'.format(i)
            tm.run(fastq1, fastq2, index)

        for i, fastq in enumerate(self.bulk1_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'bulk1.0{:0>3}'.format(i)
            tm.run(fastq1, fastq2, index)

        for i, fastq in enumerate(self.bulk2_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'bulk2.0{:0>3}'.format(i)
            tm.run(fastq1, fastq2, index)

    def alignment(self):
        aln = Alignment(args)
        for i, fastq in enumerate(self.parent_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'parent.0{:0>3}'.format(i)
            aln.run(fastq1, fastq2, index)

        for i, fastq in enumerate(self.bulk1_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'bulk1.0{:0>3}'.format(i)
            aln.run(fastq1, fastq2, index)

        for i, fastq in enumerate(self.bulk2_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'bulk2.0{:0>3}'.format(i)
            aln.run(fastq1, fastq2, index)

    def bamfilt(self):
        bt = BamFilt(self.args)
        bt.run()

    def mpileup(self):
        os.mkdir('{}/30_vcf'.format(self.out))
        mp = Mpileup(self.args)
        mp.run()

    def qtlplot(self):
        cmd = 'qtlplot -v {0}/30_vcf/qtlseq.vcf.gz \
                           -F {1} \
                           -n1 {2} \
                           -n2 {3} \
                           -t {4} \
                           -w {5} \
                           -s {6} \
                           -N {7} \
                           -D {8} \
                           -d {9} \
                           -o {0}/40_qtlseq'.format(self.out,
                                                    self.args.filial,
                                                    self.args.N_bulk1,
                                                    self.args.N_bulk2,
                                                    self.args.threads,
                                                    self.args.window,
                                                    self.args.step,
                                                    self.args.N_rep,
                                                    self.args.max_depth,
                                                    self.args.min_depth)
        if self.args.snpEff is not None:
            cmd = cmd + ' -e {}'.format(self.args.snpEff)
        
        if self.args.species is not None:
            cmd = cmd + ' --species {}'.format(self.args.species)

        cmd = clean_cmd(cmd)
        p = sbp.Popen(cmd,
                      stdout=sbp.PIPE,
                      stderr=sbp.STDOUT,
                      shell=True)

        for line in iter(p.stdout.readline, b''):
            print(line.rstrip().decode('utf-8'), flush=True)

    def run(self):
        self.refindex()
        if self.args.trim:
            self.trimming()
        else:
            self.alignment()
        self.bamfilt()
        self.mpileup()
        self.qtlplot()

def main():
    print(time_stamp(), 'start to run QTL-seq.', flush=True)
    QTLseq(args).run()
    print(time_stamp(), 'QTL-seq successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
