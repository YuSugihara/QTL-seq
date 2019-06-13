
import subprocess as sbp
from qtlseq.utils import time_stamp
from qtlseq.utils import clean_cmd


class Alignment(object):

    def __init__(self, args, config):
        self.args = args
        self.config = config
        self.out = args.out

    def run(self, fastq1, fastq2, index):
        print(time_stamp(),
              'start to align reads of {} by BWA.'.format(index),
              flush=True)

        cmd = 'bwa mem -t {0} {1} {2} {3} | \
                samtools view -b \
                              -o {4}/20_bam/{5}.bam \
                              &>> {4}/log/bwa.log'.format(self.args.threads,
                                                          self.args.ref,
                                                          fastq1,
                                                          fastq2,
                                                          self.out,
                                                          index)

        cmd = clean_cmd(cmd)

        sbp.run(cmd,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        print(time_stamp(),
              'alignment {} successfully finished.'.format(index),
              flush=True)
