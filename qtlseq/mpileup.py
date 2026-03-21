import sys
import re
import os
import glob
import subprocess as sbp
from multiprocessing import Pool
from qtlseq.utils import time_stamp, clean_cmd, call_log


class Mpileup(object):

    def __init__(self, args):
        self.out = args.out
        self.args = args

    def get_bams(self, label):
        bams = sorted(glob.glob('{}/20_bam/{}*.bam'.format(self.out, label)))
        return bams

    def merge(self):
        for label in ['parent', 'bulk1', 'bulk2']:
            bams = self.get_bams(label)
            if len(bams) == 0:
                sys.stderr.write(('{} Error: No BAM files were found for "{}".\n'
                                  '{} Expected files like "{}/20_bam/{}*.bam".\n'
                                  '{} Please check your input options and previous alignment/merge steps.\n\n').format(time_stamp(),
                                                                                                                        label,
                                                                                                                        time_stamp(),
                                                                                                                        self.out,
                                                                                                                        label,
                                                                                                                        time_stamp()))
                sys.exit(1)

            if len(bams) == 1:
                path_to_bam = os.path.abspath(bams[0])
                cmd1 = 'ln -s {} {}/20_bam/{}.unsorted.bam'.format(path_to_bam,
                                                                   self.out,
                                                                   label)
            else:
                cmd1 = 'samtools merge -f {0}/20_bam/{1}.unsorted.bam \
                                          {0}/20_bam/{1}*.bam \
                                          >> {0}/log/samtools.log \
                                          2>&1'.format(self.out, label)

            cmd2 = 'samtools collate -O \
                                     -@ {3} \
                                     {0}/20_bam/{1}.unsorted.bam \
                                     2>> {0}/log/samtools.log | \
                    samtools fixmate -m \
                                     - \
                                     - \
                                     2>> {0}/log/samtools.log | \
                    samtools sort -m {2} \
                                  -@ {3} \
                                  2>> {0}/log/samtools.log | \
                    samtools markdup -r \
                                     - \
                                     - \
                                     2>> {0}/log/samtools.log | \
                    samtools view -b \
                                  -f 2 \
                                  -F 2048 \
                                  -o {0}/20_bam/{1}.bam \
                                  >> {0}/log/samtools.log \
                                  2>&1'.format(self.out,
                                               label,
                                               self.args.mem,
                                               self.args.threads)

            cmd3 = 'samtools index {0}/20_bam/{1}.bam \
                                   >> {0}/log/samtools.log \
                                   2>&1'.format(self.out, label)

            cmd4 = 'samtools index -c \
                                   {0}/20_bam/{1}.bam \
                                   >> {0}/log/samtools.log \
                                   2>&1'.format(self.out, label)

            cmd5 = 'rm -f {}/20_bam/{}.*.bam'.format(self.out, label)

            cmd1 = clean_cmd(cmd1)
            cmd2 = clean_cmd(cmd2)
            cmd3 = clean_cmd(cmd3)
            cmd4 = clean_cmd(cmd4)
            cmd5 = clean_cmd(cmd5)

            try:
                sbp.run(cmd1,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)
            except sbp.CalledProcessError as exc:
                call_log(self.out, 'samtools', cmd1, returncode=exc.returncode)
                sys.exit(1)

            try:
                sbp.run(cmd2,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)
            except sbp.CalledProcessError as exc:
                call_log(self.out, 'samtools', cmd2, returncode=exc.returncode)
                sys.exit(1)

            try:
                sbp.run(cmd3,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)
            except sbp.CalledProcessError as exc:
                try:
                    sbp.run(cmd4,
                            stdout=sbp.DEVNULL,
                            stderr=sbp.DEVNULL,
                            shell=True,
                            check=True)
                except sbp.CalledProcessError as exc2:
                    call_log(self.out, 'samtools', cmd4, returncode=exc2.returncode)
                    sys.exit(1)

            sbp.run(cmd5,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

    def get_header(self):
        ref = open(self.args.ref, 'r')
        pattern = re.compile('>')
        chr_names = []
        for line in ref:
            if pattern.match(line):
                line = line.rstrip('\n')
                chr_name = re.split('[> ]',line)[1]
                chr_names.append(chr_name)
        ref.close()
        if len(chr_names) == 0:
            sys.stderr.write(('{} Error: No contig headers were found in reference FASTA "{}".\n'
                              '{} Please verify FASTA format and that each contig starts with ">".\n\n').format(time_stamp(),
                                                                                                                  self.args.ref,
                                                                                                                  time_stamp()))
            sys.exit(1)
        return chr_names

    def mpileup(self, chr_name):
        cmd = 'bcftools mpileup -a AD,ADF,ADR \
                                -B \
                                -q {0}\
                                -Q {1} \
                                -C {2} \
                                -O u \
                                -r {3} \
                                -f {4} \
                                --ignore-RG \
                                {5}/20_bam/parent.bam \
                                {5}/20_bam/bulk1.bam \
                                {5}/20_bam/bulk2.bam \
                                2>> {5}/log/bcftools.{3}.log | \
               bcftools call -vm \
                             -f GQ,GP \
                             -O u \
                             2>> {5}/log/bcftools.{3}.log | \
               bcftools filter -i "INFO/MQ>={0}" \
                               -O z \
                               -o {5}/30_vcf/qtlseq.{3}.vcf.gz \
                               >> {5}/log/bcftools.{3}.log \
                               2>&1'.format(self.args.min_MQ,
                                            self.args.min_BQ,
                                            self.args.adjust_MQ,
                                            chr_name,
                                            self.args.ref,
                                            self.out)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError as exc:
            call_log(self.out,
                     'bcftools.{}'.format(chr_name),
                     cmd,
                     returncode=exc.returncode,
                     hint='If exit status is 127, please install bcftools and check PATH.')
            raise RuntimeError('bcftools mpileup failed for "{}"'.format(chr_name))

    def concat(self):

        cmd1 = 'bcftools concat -O z \
                                -o {0}/30_vcf/qtlseq.vcf.gz \
                                {0}/30_vcf/qtlseq.*.vcf.gz \
                                >> {0}/log/bcftools.log \
                                2>&1'.format(self.out)
        cmd2 = 'cat {0}/log/bcftools.*.log > {0}/log/bcftools.log'.format(self.out)
        cmd3 = 'rm -f {}/30_vcf/qtlseq.*.vcf.gz'.format(self.out)
        cmd4 = 'rm -f {}/log/bcftools.*.log'.format(self.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)

        vcf_chunks = sorted(glob.glob('{}/30_vcf/qtlseq.*.vcf.gz'.format(self.out)))
        if len(vcf_chunks) == 0:
            sys.stderr.write(('{} Error: No per-contig VCF files were generated.\n'
                              '{} Please check "{}/log/bcftools.*.log".\n\n').format(time_stamp(),
                                                                                      time_stamp(),
                                                                                      self.out))
            sys.exit(1)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError as exc:
            call_log(self.out, 'bcftools', cmd1, returncode=exc.returncode)
            sys.exit(1)

        try:
            sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
            sbp.run(cmd3, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
            sbp.run(cmd4, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        except sbp.CalledProcessError as exc:
            call_log(self.out, 'bcftools', cmd2, returncode=exc.returncode)
            sys.exit(1)

    def mkindex(self):
        cmd1 = 'tabix -f \
                      -p vcf \
                      {0}/30_vcf/qtlseq.vcf.gz \
                      >> {0}/log/tabix.log \
                      2>&1'.format(self.out)
        
        cmd2 = 'tabix -f \
                      -C \
                      -p vcf \
                      {0}/30_vcf/qtlseq.vcf.gz \
                      >> {0}/log/tabix.log \
                      2>&1'.format(self.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError as exc:
            try:
                sbp.run(cmd2, 
                        stdout=sbp.DEVNULL, 
                        stderr=sbp.DEVNULL, 
                        shell=True, 
                        check=True)
            except sbp.CalledProcessError as exc2:
                call_log(self.out, 'tabix', cmd2, returncode=exc2.returncode)
                sys.exit(1)

    def run(self):
        print(time_stamp(), 'Sorting BAM files.', flush=True)
        self.merge()
        print(time_stamp(), 'BAM file sorting completed successfully.', flush=True)

        print(time_stamp(), 'Calling variants.', flush=True)
        chr_names = self.get_header()

        try:
            with Pool(self.args.threads) as p:
                p.map(self.mpileup, chr_names)
        except Exception:
            sys.stderr.write(('{} Error: Variant calling failed during per-contig mpileup.\n'
                              '{} Please check "{}/log/bcftools.*.log".\n\n').format(time_stamp(),
                                                                                      time_stamp(),
                                                                                      self.out))
            sys.exit(1)

        self.concat()
        print(time_stamp(), 'Variant calling completed successfully.', flush=True)

        print(time_stamp(), 'Indexing the VCF file.', flush=True)
        self.mkindex()
        print(time_stamp(), 'VCF indexing completed successfully.', flush=True)
