
import re
import os
import glob
import subprocess as sbp
from multiprocessing import Pool
from qtlseq.utils import time_stamp
from qtlseq.utils import clean_cmd
from qtlseq.vcf2index import Vcf2Index


class Mpileup(object):

    def __init__(self, args):
        self.out = args.out
        self.args = args

    def get_bams(self, label):
        bams = glob.glob('{}/20_bam/{}*.filt.bam'.format(self.out, label))
        return bams

    def merge(self):
        for label in ['parent', 'bulk1', 'bulk2']:
            bams = self.get_bams(label)
            if len(bams) == 1:
                path_to_bam = os.path.abspath(bams[0])
                cmd1 = 'ln -s {} {}/20_bam/{}.unsorted.filt.bam'.format(path_to_bam,
                                                                        self.out,
                                                                        label)
            else:
                cmd1 = 'samtools merge -f {0}/20_bam/{1}.unsorted.filt.bam \
                                          {0}/20_bam/{1}*.filt.bam'.format(self.out,
                                                                           label)

            cmd2 = 'samtools sort -m {0} \
                                  -@ {1} \
                                  -o {2}/20_bam/{3}.filt.bam \
                                  {2}/20_bam/{3}.unsorted.filt.bam \
                                  &>> {2}/log/samtools.log'.format(self.args.mem,
                                                                   self.args.threads,
                                                                   self.out,
                                                                   label)

            cmd3 = 'samtools index {0}/20_bam/{1}.filt.bam \
                                   &>> {0}/log/samtools.log'.format(self.out,
                                                                    label)

            cmd4 = 'rm -f {}/20_bam/{}.*.filt.bam'.format(self.out, label)

            cmd1 = clean_cmd(cmd1)
            cmd2 = clean_cmd(cmd2)
            cmd3 = clean_cmd(cmd3)
            cmd4 = clean_cmd(cmd4)

            sbp.run(cmd1, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
            sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
            sbp.run(cmd3, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
            sbp.run(cmd4, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def get_header(self):
        ref = open(self.args.ref, 'r')
        pattern = re.compile('>')
        chr_names = []
        for line in ref:
            if pattern.match(line):
                line = line.rstrip('\n')
                chr_name = re.split('[> ]',line)[1]
                chr_names.append(chr_name)
        return chr_names

    def mpileup(self, chr_name):
        cmd1 = 'bcftools mpileup -a AD,ADF,ADR \
                                 -B \
                                 -q {0}\
                                 -Q {1} \
                                 -C {2} \
                                 -O u \
                                 -r {3} \
                                 -f {4} \
                                 {5}/20_bam/parent.filt.bam \
                                 {5}/20_bam/bulk1.filt.bam \
                                 {5}/20_bam/bulk2.filt.bam | \
                bcftools call -vm \
                              -f GQ,GP \
                              -O u | \
                bcftools filter -i "INFO/MQ>={0}" \
                                -O z \
                                -o {5}/30_vcf/qtlseq.{3}.vcf.gz \
                                &>> {5}/log/bcftools.{3}.log'.format(self.args.min_MQ,
                                                                     self.args.min_BQ,
                                                                     self.args.adjust_MQ,
                                                                     chr_name,
                                                                     self.args.ref,
                                                                     self.out)

        cmd2 = 'tabix -f \
                      -p vcf \
                      {0}/30_vcf/qtlseq.{1}.vcf.gz \
                      &>> {0}/log/tabix.{1}.log'.format(self.out, chr_name)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        sbp.run(cmd1, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def concat(self):
        cmd1 = 'cat {0}/log/bcftools.*.log > {0}/log/bcftools.log'.format(self.out)
        cmd2 = 'cat {0}/log/tabix.*.log > {0}/log/tabix.log'.format(self.out)

        cmd3 = 'bcftools concat -a \
                                -O z \
                                -o {0}/30_vcf/qtlseq.vcf.gz \
                                {0}/30_vcf/qtlseq.*.vcf.gz \
                                &>> {0}/log/bcftools.log'.format(self.out)

        cmd4 = 'rm -f {}/30_vcf/qtlseq.*.vcf.gz'.format(self.out)
        cmd5 = 'rm -f {}/30_vcf/qtlseq.*.vcf.gz.tbi'.format(self.out)
        cmd6 = 'rm -f {}/log/bcftools.*.log'.format(self.out)
        cmd7 = 'rm -f {}/log/tabix.*.log'.format(self.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)
        cmd5 = clean_cmd(cmd5)
        cmd6 = clean_cmd(cmd6)
        cmd7 = clean_cmd(cmd7)

        sbp.run(cmd1, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd3, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd4, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd5, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd6, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd7, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def mkindex(self):
        cmd = 'tabix -f \
                     -p vcf \
                     {0}/30_vcf/qtlseq.vcf.gz \
                     &>> {0}/log/tabix.log'.format(self.out)

        cmd = clean_cmd(cmd)
        sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def run(self):
        print(time_stamp(), 'start to merge BAMs.', flush=True)
        self.merge()
        print(time_stamp(), 'merging process successfully finished.', flush=True)

        print(time_stamp(), 'start to call variants.', flush=True)
        chr_names = self.get_header()

        p = Pool(self.args.threads)
        p.map(self.mpileup, chr_names)
        p.close()

        self.concat()
        print(time_stamp(), 'variant calling successfully finished.', flush=True)

        print(time_stamp(), 'start to index VCF.', flush=True)
        self.mkindex()
        print(time_stamp(), 'indexing VCF successfully finished.', flush=True)
