#!/usr/bin/env python3

from qtlseq.params import Params


pm = Params('qtlseq')
args = pm.set_options()


import os
import sys
import subprocess as sbp
from multiprocessing import set_start_method
try:
    set_start_method('fork')
except (RuntimeError, ValueError):
    pass
# These environment variables can avoid the potential errors
# in python and java
if os.environ.get('OMP_NUM_THREADS') is None:
    os.environ['OMP_NUM_THREADS'] = str(args.threads)
if os.environ.get('USE_SIMPLE_THREADED_LEVEL3') is None:
    os.environ['USE_SIMPLE_THREADED_LEVEL3'] = str(args.threads)
if os.environ.get('JAVA_TOOL_OPTIONS') is None:
    os.environ['JAVA_TOOL_OPTIONS'] = '-Xmx4000m'
try:
    from qtlseq.refindex import RefIndex
    from qtlseq.trim import Trim
    from qtlseq.alignment import Alignment
    from qtlseq.mpileup import Mpileup
    from qtlseq.utils import time_stamp, check_required_commands
except ModuleNotFoundError as exc:
    sys.stderr.write(('Error: Missing Python dependency "{}".\n'
                      'Please reinstall QTL-seq with "pip install -e ." or use the conda package.\n\n').format(exc.name))
    sys.exit(1)


class QTLseq(object):

    def __init__(self, args):
        args = pm.check_max_threads(args)
        self.N_fastq = pm.check_args(args)
        self.args = args
        self.out = args.out
        self.check_dependencies()
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

    def check_dependencies(self):
        required = ['bwa', 'samtools', 'bcftools', 'tabix']
        if self.args.trim:
            required.append('trimmomatic')
        if self.args.snpEff is not None:
            required.append('snpEff')
        check_required_commands(required, 'QTL-seq')

    def mkdir(self):
        out_dirs = ['{}'.format(self.out),
                    '{}/log'.format(self.out),
                    '{}/10_ref'.format(self.out),
                    '{}/20_bam'.format(self.out)]
        if self.N_fastq > 0 and self.args.trim:
            out_dirs.insert(2, '{}/00_fastq'.format(self.out))

        for out_dir in out_dirs:
            try:
                os.mkdir(out_dir)
            except OSError as exc:
                sys.stderr.write(('{} Error: Failed to create directory "{}": {}\n'
                                  '{} Please check permissions and path.\n\n').format(time_stamp(),
                                                                                      out_dir,
                                                                                      exc,
                                                                                      time_stamp()))
                sys.exit(1)

    def link_ref(self):
        path_to_ref = os.path.abspath(self.args.ref)
        ref = os.path.basename(self.args.ref)
        sym_ref = '{}/10_ref/{}'.format(self.out, ref)
        try:
            os.symlink(path_to_ref, sym_ref)
        except OSError as exc:
            sys.stderr.write(('{} Error: Failed to link reference FASTA.\n'
                              '{} source: "{}"\n'
                              '{} target: "{}"\n'
                              '{} reason: {}\n\n').format(time_stamp(),
                                                           time_stamp(),
                                                           path_to_ref,
                                                           time_stamp(),
                                                           sym_ref,
                                                           time_stamp(),
                                                           exc))
            sys.exit(1)
        self.args.ref = sym_ref

    def link_bam(self, label, input_files):
        for i, input_name in enumerate(input_files):
            path_to_bam = os.path.abspath(input_name)
            index = '{}.1{:0>3}'.format(label, i)
            sym_bam = '{}/20_bam/{}.bam'.format(self.out, index)
            try:
                os.symlink(path_to_bam, sym_bam)
            except OSError as exc:
                sys.stderr.write(('{} Error: Failed to link BAM file.\n'
                                  '{} source: "{}"\n'
                                  '{} target: "{}"\n'
                                  '{} reason: {}\n\n').format(time_stamp(),
                                                               time_stamp(),
                                                               path_to_bam,
                                                               time_stamp(),
                                                               sym_bam,
                                                               time_stamp(),
                                                               exc))
                sys.exit(1)

    def refindex(self):
        ri = RefIndex(self.args)
        ri.run()

    def trimming(self):
        tm = Trim(self.args)
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
        aln = Alignment(self.args)
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
        vcf_dir = '{}/30_vcf'.format(self.out)
        try:
            os.mkdir(vcf_dir)
        except OSError as exc:
            sys.stderr.write(('{} Error: Failed to create VCF output directory "{}": {}\n'
                              '{} Please check permissions and path.\n\n').format(time_stamp(),
                                                                                  vcf_dir,
                                                                                  exc,
                                                                                  time_stamp()))
            sys.exit(1)
        mp = Mpileup(self.args)
        mp.run()

    def qtlplot(self):
        cmd = [sys.executable,
               '-m',
               'qtlseq.qtlplot',
               '-v',
               '{}/30_vcf/qtlseq.vcf.gz'.format(self.out),
               '-F',
               str(self.args.filial),
               '-n1',
               str(self.args.N_bulk1),
               '-n2',
               str(self.args.N_bulk2),
               '-t',
               str(self.args.threads),
               '-w',
               str(self.args.window),
               '-s',
               str(self.args.step),
               '-N',
               str(self.args.N_rep),
               '-D',
               str(self.args.max_depth),
               '-d',
               str(self.args.min_depth),
               '--dot-colors',
               self.args.dot_colors,
               '--line-colors',
               self.args.line_colors,
               '-o',
               '{}/40_qtlseq'.format(self.out)]
        if self.args.snpEff is not None:
            cmd.extend(['-e', self.args.snpEff])

        p = sbp.Popen(cmd,
                      stdout=sbp.PIPE,
                      stderr=sbp.STDOUT,
                      text=True)

        for line in iter(p.stdout.readline, ''):
            print(line.rstrip(), flush=True)

        status = p.wait()
        if status != 0:
            sys.stderr.write(('{} Error: qtlplot failed with exit code {}.\n\n').format(time_stamp(), status))
            sys.stderr.flush()
            sys.exit(status)

    def run(self):
        self.refindex()
        if self.args.trim:
            self.trimming()
        else:
            self.alignment()
        self.mpileup()
        self.qtlplot()

def main():
    print(time_stamp(), 'Starting QTL-seq analysis.', flush=True)
    try:
        QTLseq(args).run()
    except KeyboardInterrupt:
        sys.stderr.write('\n{} Error: Interrupted by user.\n\n'.format(time_stamp()))
        sys.exit(1)
    except SystemExit:
        raise
    except Exception as exc:
        sys.stderr.write(('{} Error: Unexpected failure in QTL-seq: {}\n'
                          '{} Please check log files under "{}/log" if they were generated.\n\n').format(time_stamp(),
                                                                                                           exc,
                                                                                                           time_stamp(),
                                                                                                           args.out))
        sys.exit(1)
    print(time_stamp(), 'QTL-seq analysis completed successfully.\n', flush=True)

if __name__ == '__main__':
    main()
