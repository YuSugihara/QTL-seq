import sys
import subprocess as sbp
from qtlseq.utils import time_stamp, clean_cmd, call_log


class Alignment(object):

    def __init__(self, args):
        self.args = args
        self.out = args.out

    def run(self, fastq1, fastq2, index):
        print(time_stamp(),
              'start to align reads of {} by BWA.'.format(index),
              flush=True)

        cmd = 'bwa mem -t {0} \
                       {1} {2} {3} | \
               samtools fixmate -m \
                                - \
                                - | \
               samtools sort -m {4} \
                             -@ {0} | \
               samtools markdup -r \
                                - \
                                - | \
               samtools view -b \
                             -f 3 \
                             -o {5}/20_bam/{6}.bam \
                             >> {5}/log/alignment.log \
                             2>&1'.format(self.args.threads,
                                          self.args.ref,
                                          fastq1,
                                          fastq2,
                                          self.args.mem,
                                          self.out,
                                          index)
        

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'alignment', cmd)
            sys.exit(1)

        print(time_stamp(),
              'alignment {} successfully finished.'.format(index),
              flush=True)
