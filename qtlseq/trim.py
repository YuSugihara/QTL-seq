import sys
import os
import subprocess as sbp
from qtlseq.alignment import Alignment
from qtlseq.utils import time_stamp, clean_cmd, call_log


class Trim(object):

    def __init__(self, args):
        self.out = args.out
        self.args = args
        self.trim_params = self.params_parser(args.trim_params)

    def params_parser(self, trim_params):
        params_list = trim_params.split(',')
        trim_params = {}
        trim_params['phred'] = params_list[0]
        trim_params['ILLUMINACLIP'] = params_list[1]
        trim_params['LEADING'] = params_list[2]
        trim_params['TRAILING'] = params_list[3]
        trim_params['SLIDINGWINDOW'] = params_list[4]
        trim_params['MINLEN'] = params_list[5]

        if self.args.adapter is not None:
            ILLUMINACLIP = trim_params['ILLUMINACLIP'].split(':')
            ILLUMINACLIP[0] = self.args.adapter
            trim_params['ILLUMINACLIP'] = ':'.join(ILLUMINACLIP)
            
        return trim_params

    def run(self, fastq1, fastq2, index):
        print(time_stamp(),
        'Trimming reads for {} and {}.'.format(fastq1, fastq2),
        flush=True)

        trim1 = '{}/00_fastq/{}.1.trim.fastq.gz'.format(self.out, index)
        trim2 = '{}/00_fastq/{}.2.trim.fastq.gz'.format(self.out, index)
        unpaired1 = '{}/00_fastq/{}.1.unpaired.fastq.gz'.format(self.out, index)
        unpaired2 = '{}/00_fastq/{}.2.unpaired.fastq.gz'.format(self.out, index)

        if (len(self.trim_params['ILLUMINACLIP']) == 0) or \
           ('<ADAPTER_FASTA>' in self.trim_params['ILLUMINACLIP']):
            cmd = 'trimmomatic PE -threads {} \
                                  -phred{} {} {} {} {} {} {} \
                                  LEADING:{} \
                                  TRAILING:{} \
                                  SLIDINGWINDOW:{} \
                                  MINLEN:{} \
                                  >> {}/log/trimmomatic.log \
                                  2>&1'.format(self.args.threads,
                                               self.trim_params['phred'],
                                               fastq1,
                                               fastq2,
                                               trim1,
                                               unpaired1,
                                               trim2,
                                               unpaired2,
                                               self.trim_params['LEADING'],
                                               self.trim_params['TRAILING'],
                                               self.trim_params['SLIDINGWINDOW'],
                                               self.trim_params['MINLEN'],
                                               self.out)
        else:
            cmd = 'trimmomatic PE -threads {} \
                                  -phred{} {} {} {} {} {} {} \
                                  ILLUMINACLIP:{} \
                                  LEADING:{} \
                                  TRAILING:{} \
                                  SLIDINGWINDOW:{} \
                                  MINLEN:{} \
                                  >> {}/log/trimmomatic.log \
                                  2>&1'.format(self.args.threads,
                                               self.trim_params['phred'],
                                               fastq1,
                                               fastq2,
                                               trim1,
                                               unpaired1,
                                               trim2,
                                               unpaired2,
                                               self.trim_params['ILLUMINACLIP'],
                                               self.trim_params['LEADING'],
                                               self.trim_params['TRAILING'],
                                               self.trim_params['SLIDINGWINDOW'],
                                               self.trim_params['MINLEN'],
                                               self.out)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'trimmomatic', cmd)
            sys.exit(1)

        print(time_stamp(),
              'Trimming of {} and {} completed successfully.'.format(fastq1, fastq2),
              flush=True)

        aln = Alignment(self.args)
        aln.run(trim1, trim2, index)
