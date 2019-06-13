
import os
import subprocess as sbp
from qtlseq.utils import time_stamp
from qtlseq.utils import clean_cmd
from qtlseq.alignment import Alignment


class Trim(object):

    def __init__(self, args, config):
        self.out = args.out
        self.args = args
        self.config = config

    def run(self, fastq1, fastq2, index):
        print(time_stamp(),
        'start trimming for {} and {}.'.format(fastq1, fastq2),
        flush=True)

        trim1 = '{}/00_fastq/{}.1.trim.fastq.gz'.format(self.out,
                                                        index)
        trim2 = '{}/00_fastq/{}.2.trim.fastq.gz'.format(self.out,
                                                        index)
        unpaired1 = '{}/00_fastq/{}.1.unpaired..fastq.gz'.format(self.out,
                                                                 index)
        unpaired2 = '{}/00_fastq/{}.2.unpaired..fastq.gz'.format(self.out,
                                                                 index)

        cmd = 'trimmomatic PE -threads {} \
                              -phred{} {} {} {} {} {} {} \
                              ILLUMINACLIP:{} \
                              LEADING:{} \
                              TRAILING:{} \
                              SLIDINGWINDOW:{} \
                              MINLEN:{} \
                              &>> {}/log/trimmomatic.log'.format(self.args.threads,
                                                                 self.config['trimmomatic']['phred'],
                                                                 fastq1,
                                                                 fastq2,
                                                                 trim1,
                                                                 unpaired1,
                                                                 trim2,
                                                                 unpaired2,
                                                                 self.config['trimmomatic']['ILLUMINACLIP'],
                                                                 self.config['trimmomatic']['LEADING'],
                                                                 self.config['trimmomatic']['TRAILING'],
                                                                 self.config['trimmomatic']['SLIDINGWINDOW'],
                                                                 self.config['trimmomatic']['MINLEN'],
                                                                 self.out)
        cmd = clean_cmd(cmd)
        sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

        print(time_stamp(),
              'trimming for {} and {} successfully finished.'.format(fastq1, fastq2),
              flush=True)

        aln = Alignment(self.args, self.config)
        aln.run(trim1, trim2, index)
