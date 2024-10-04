import sys
import subprocess as sbp
from qtlseq.utils import time_stamp, clean_cmd, call_log


class RefIndex(object):

    def __init__(self, args):
        self.args = args
        self.out = args.out

    def run(self):
        print(time_stamp(),
              'start to index reference fasta.',
              flush=True)

        cmd1 = 'bwa index {} \
                >> {}/log/bwa.log \
                2>&1'.format(self.args.ref, self.out)

        cmd2 = 'samtools faidx {} \
                >> {}/log/samtools.log \
                2>&1'.format(self.args.ref, self.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'bwa', cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'samtools', cmd1)
            sys.exit(1)

        print(time_stamp(),
              'indexing of reference successfully finished.',
              flush=True)