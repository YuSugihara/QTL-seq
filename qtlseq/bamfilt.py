
import os
import glob
import subprocess as sbp
from multiprocessing import Pool
from qtlseq.utils import time_stamp, clean_cmd, call_log


class BamFilt(object):

    def __init__(self, args):
        self.out = args.out
        self.args = args

    def get_labels(self):
        samples = glob.glob('{}/20_bam/*.bam'.format(self.out))
        labels = []
        label_with_flags = []
        for sample in samples:
            sample = os.path.basename(sample)
            label, ext = os.path.splitext(sample)
            labels.append(label)
            for flag in [83, 99, 147, 163]:
                label_with_flags.append((flag, label))
        return labels, label_with_flags

    def filt(self, label_with_flags):
        flag = label_with_flags[0]
        label = label_with_flags[1]
        cmd = 'samtools view -b \
                             -f {0} \
                             -o {1}/20_bam/{2}.f{0}.bam \
                             {1}/20_bam/{2}.bam \
                             >> {1}/log/samtools.{2}.f{0}.log \
                             2>&1'.format(flag,
                                          self.out,
                                          label)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'samtools', cmd)
            sys.exit()

    def merge(self, label):
        cmd1 = 'cat {0}/log/samtools.{1}.f83.log \
                    {0}/log/samtools.{1}.f99.log \
                    {0}/log/samtools.{1}.f147.log \
                    {0}/log/samtools.{1}.f163.log \
                    > {0}/log/samtools.{1}.log'.format(self.out, label)

        cmd2 = 'rm -f {0}/20_bam/{1}.bam'.format(self.out, label)
        cmd3 = 'rm -f {0}/log/samtools.{1}.f*.log'.format(self.out, label)

        cmd4 = 'samtools merge -f {0}/20_bam/{1}.filt.bam \
                                  {0}/20_bam/{1}.f83.bam \
                                  {0}/20_bam/{1}.f99.bam \
                                  {0}/20_bam/{1}.f147.bam \
                                  {0}/20_bam/{1}.f163.bam \
                                  >> {0}/log/samtools.{1}.log \
                                  2>&1'.format(self.out, label)

        cmd5 = 'rm -f {}/20_bam/{}.f83.bam'.format(self.out, label)
        cmd6 = 'rm -f {}/20_bam/{}.f99.bam'.format(self.out, label)
        cmd7 = 'rm -f {}/20_bam/{}.f147.bam'.format(self.out, label)
        cmd8 = 'rm -f {}/20_bam/{}.f163.bam'.format(self.out, label)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)
        cmd5 = clean_cmd(cmd5)
        cmd6 = clean_cmd(cmd6)
        cmd7 = clean_cmd(cmd7)
        cmd8 = clean_cmd(cmd8)

        sbp.run(cmd1, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd3, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

        try:
            sbp.run(cmd4, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'samtools', cmd4)
            sys.exit()

        sbp.run(cmd5, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd6, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd7, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd8, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def clean_log(self):
        cmd1 = 'cat {0}/log/samtools.*.log >> {0}/log/samtools.log'.format(self.out)
        cmd2 = 'rm -f {}/log/samtools.*.log'.format(self.out)
        sbp.run(cmd1, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def run(self):
        labels, label_with_flags = self.get_labels()

        print(time_stamp(), 'start to filter reads.', flush=True)
        p = Pool(self.args.threads)
        p.map(self.filt, label_with_flags)
        p.close()

        p = Pool(self.args.threads)
        p.map(self.merge, labels)
        p.close()

        self.clean_log()
        print(time_stamp(), 'filtering process successfully finished.', flush=True)
