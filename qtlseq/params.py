import os
import sys
import argparse
from multiprocessing import Pool
import multiprocessing as multi
from qtlseq.utils import time_stamp
from qtlseq.__init__ import __version__


class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'qtlseq':
            parser = self.qtlseq_options()
        elif self.program_name == 'qtlplot':
            parser = self.qtlplot_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def qtlseq_options(self):
        parser = argparse.ArgumentParser(description='QTL-seq version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('qtlseq -r <FASTA> -p <BAM|FASTQ> -b1 <BAM|FASTQ>\n'
                        '              -b2 <BAM|FASTQ> -n1 <INT> -n2 <INT> -o <OUT_DIR>\n'
                        '              [-F <INT>] [-T] [-e <DATABASE>]')

        # set options
        parser.add_argument('-r',
                            '--ref',
                            action='store',
                            required=True,
                            type=str,
                            help='Reference fasta.',
                            metavar='')

        parser.add_argument('-p',
                            '--parent',
                            action='append',
                            required=True,
                            type=str,
                            help=('fastq or bam of parent. If you specify\n'
                                  'fastq, please separate pairs by comma,\n'
                                  'e.g. -p fastq1,fastq2. You can use this\n'
                                  'optiion multiple times'),
                            metavar='')

        parser.add_argument('-b1',
                            '--bulk1',
                            action='append',
                            required=True,
                            type=str,
                            help=('fastq or bam of bulk 1. If you specify\n'
                                  'fastq, please separate pairs by comma,\n'
                                  'e.g. -b1 fastq1,fastq2. You can use this\n'
                                  'optiion multiple times'),
                            metavar='')

        parser.add_argument('-b2',
                            '--bulk2',
                            action='append',
                            required=True,
                            type=str,
                            help=('fastq or bam of bulk 2. If you specify\n'
                                  'fastq, please separate pairs by comma,\n'
                                  'e.g. -b2 fastq1,fastq2. You can use this\n'
                                  'optiion multiple times'),
                            metavar='')

        parser.add_argument('-n1',
                            '--N-bulk1',
                            action='store',
                            required=True,
                            type=int,
                            help='Number of individuals in bulk 1.',
                            metavar='')

        parser.add_argument('-n2',
                            '--N-bulk2',
                            action='store',
                            required=True,
                            type=int,
                            help='Number of individuals in bulk 2.',
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help=('Output directory. Specified name must not\n'
                                  'exist.'),
                            metavar='')

        parser.add_argument('-F',
                            '--filial',
                            action='store',
                            default=2,
                            type=int,
                            help=('Filial generation. This parameter must be\n'
                                  'more than 1. [2]'),
                            metavar='')

        parser.add_argument('-t',
                            '--threads',
                            action='store',
                            default=2,
                            type=int,
                            help=('Number of threads. If you specify the number\n'
                                  'below one, then QTL-seq will use the threads\n'
                                  'as many as possible. [2]'),
                            metavar='')

        parser.add_argument('-w',
                            '--window',
                            action='store',
                            default=2000,
                            type=int,
                            help='Window size (kb). [2000]',
                            metavar='')

        parser.add_argument('-s',
                            '--step',
                            action='store',
                            default=100,
                            type=int,
                            help='Step size (kb). [100]',
                            metavar='')

        parser.add_argument('-D',
                            '--max-depth',
                            action='store',
                            default=250,
                            type=int,
                            help='Maximum depth of variants which will be used. [250]',
                            metavar='')

        parser.add_argument('-d',
                            '--min-depth',
                            action='store',
                            default=8,
                            type=int,
                            help='Minimum depth of variants which will be used. [8]',
                            metavar='')

        parser.add_argument('-N',
                            '--N-rep',
                            action='store',
                            default=5000,
                            type=int,
                            help=('Number of replicates for simulation to make \n'
                                  'null distribution. [5000]'),
                            metavar='')

        parser.add_argument('-T',
                            '--trim',
                            action='store_true',
                            default=False,
                            help='Trim fastq using trimmomatic.')

        parser.add_argument('-a',
                            '--adapter',
                            action='store',
                            type=str,
                            help=('FASTA of adapter sequences. This will be used\n'
                                  'when you specify "-T" for trimming.'),
                            metavar='')

        parser.add_argument('--trim-params',
                            action='store',
                            type=str,
                            help=('Parameters for trimmomatic. Input parameters\n'
                                  'must be separated by comma with following\n'
                                  'order: phred, ILLUMINACLIP, LEADING, TRAILING,\n'
                                  'SLIDINGWINDOW, MINLEN. If you want to remove\n'
                                  'adapters of illumina, please specify FASTA of\n'
                                  'the adapter sequences with "--adapter". Specified\n'
                                  'name will be inserted into <ADAPTER_FASTA>. If you\n'
                                  "don't specify it, adapter trimming will be skipped.\n"
                                  '[33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75]'),
                            metavar='')

        parser.add_argument('-e',
                            '--snpEff',
                            action='store',
                            type=str,
                            help=('Predict causal variant using SnpEff. Please\n'
                                  'check available databases in SnpEff.'),
                            metavar='')

        parser.add_argument('--mem',
                            action='store',
                            default='1G',
                            type=str,
                            help=('Maximum memory per thread when bam sorted;\n'
                                  'suffix K/M/G recognized. [1G]'),
                            metavar='')

        parser.add_argument('-q',
                            '--min-MQ',
                            action='store',
                            default=40,
                            type=int,
                            help='Minimum mapping quality in mpileup. [40]',
                            metavar='')

        parser.add_argument('-Q',
                            '--min-BQ',
                            action='store',
                            default=18,
                            type=int,
                            help='Minimum base quality in mpileup. [18]',
                            metavar='')

        parser.add_argument('-C',
                            '--adjust-MQ',
                            action='store',
                            default=50,
                            type=int,
                            help=('"adjust-MQ" in mpileup. Default parameter\n'
                                  'is suited for BWA. [50]'),
                            metavar='')

        # set version
        parser.add_argument('-v',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    def qtlplot_options(self):
        parser = argparse.ArgumentParser(description='QTL-plot version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('qtlplot -v <VCF> -n1 <INT> -n2 <INT> -o <OUT_DIR>\n'
                        '               [-F <INT>] [-t <INT>] [-w <INT>] [-s <INT>] [-D <INT>]\n'
                        '               [-d <INT>] [-N <INT>] [-m <FLOAT>] [-S <INT>] [-e <DATABASE>]\n'
                        '               [--igv] [--corr] [--indel]')

        # set options
        parser.add_argument('-v',
                            '--vcf',
                            action='store',
                            required=True,
                            type=str,
                            help=('VCF which contains parent, bulk1 and bulk2\n'
                                  'in this order.'),
                            metavar='')

        parser.add_argument('-n1',
                            '--N-bulk1',
                            action='store',
                            required=True,
                            type=int,
                            help='Number of individuals in bulk 1.',
                            metavar='')

        parser.add_argument('-n2',
                            '--N-bulk2',
                            action='store',
                            required=True,
                            type=int,
                            help='Number of individuals in bulk 2.',
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help='Output directory. Specified name can exist.',
                            metavar='')

        parser.add_argument('-F',
                            '--filial',
                            action='store',
                            default=2,
                            type=int,
                            help=('Filial generation. This parameter must be\n'
                                  'more than 1. [2]'),
                            metavar='')

        parser.add_argument('-t',
                            '--threads',
                            action='store',
                            default=2,
                            type=int,
                            help=('Number of threads. If you specify the number\n'
                                  'below one, then QTL-plot will use the threads\n'
                                  'as many as possible. [2]'),
                            metavar='')

        parser.add_argument('-w',
                            '--window',
                            action='store',
                            default=2000,
                            type=int,
                            help='Window size (kb). [2000]',
                            metavar='')

        parser.add_argument('-s',
                            '--step',
                            action='store',
                            default=100,
                            type=int,
                            help='Step size (kb). [100]',
                            metavar='')

        parser.add_argument('-D',
                            '--max-depth',
                            action='store',
                            default=250,
                            type=int,
                            help='Maximum depth of variants which will be used. [250]',
                            metavar='')

        parser.add_argument('-d',
                            '--min-depth',
                            action='store',
                            default=8,
                            type=int,
                            help='Minimum depth of variants which will be used. [8]',
                            metavar='')

        parser.add_argument('-N',
                            '--N-rep',
                            action='store',
                            default=5000,
                            type=int,
                            help=('Number of replicates for simulation to make \n'
                                  'null distribution. [5000]'),
                            metavar='')

        parser.add_argument('-m',
                            '--min-SNPindex',
                            action='store',
                            default=0.3,
                            type=float,
                            help='Cutoff of minimum SNP-index for clear results. [0.3]',
                            metavar='')

        parser.add_argument('-S',
                            '--strand-bias',
                            action='store',
                            default=7,
                            type=int,
                            help=('Filter spurious homo genotypes in cultivar using\n'
                                  'strand bias. If ADF (or ADR) is higher than this\n'
                                  'cutoff when ADR (or ADF) is 0, that SNP will be\n'
                                  'filtered out. If you want to supress this filtering,\n'
                                  'please set this cutoff to 0. [7]\n'),
                            metavar='')

        parser.add_argument('-e',
                            '--snpEff',
                            action='store',
                            type=str,
                            help=('Predict causal variant using SnpEff. Please\n'
                                  'check available databases in SnpEff.'),
                            metavar='')

        parser.add_argument('--igv',
                            action='store_true',
                            default=False,
                            help='Output IGV format file to check results on IGV.')

        parser.add_argument('--corr',
                            action='store_true',
                            default=False,
                            help=('Use corrected threshold in Huang et al. (2019).\n'
                                  'If you specify this option, the options related to\n'
                                  'the simulation ("-N" and "-F") will be ignored.'))

        parser.add_argument('--indel',
                            action='store_true',
                            default=False,
                            help='Plot SNP-index with INDEL.')

        parser.add_argument('--fig-width',
                            action='store',
                            default=7.5,
                            type=float,
                            help='Width allocated in chromosome figure. [7.5]',
                            metavar='')

        parser.add_argument('--fig-height',
                            action='store',
                            default=4.0,
                            type=float,
                            help='Height allocated in chromosome figure. [4.0]',
                            metavar='')

        parser.add_argument('--white-space',
                            action='store',
                            default=0.6,
                            type=float,
                            help=('White space between figures. (This option\n'
                                  'only affect vertical direction.) [0.6]'),
                            metavar='')

        # set version
        parser.add_argument('--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    def check_max_threads(self, args):
        max_cpu = multi.cpu_count()
        print(time_stamp(),
              'maximum number of threads which you can use is up to {}.'.format(max_cpu),
              flush=True)
        if max_cpu <= args.threads:
            sys.stderr.write(('!!WARNING!! You can use up to {0} threads. '
                              'This program will use {0} threads.\n').format(max_cpu))
            sys.stderr.flush()
            args.threads = max_cpu
        elif args.threads < 1:
            args.threads = max_cpu
        return args

    def check_args(self, args):
        if os.path.isdir(args.out):
            sys.stderr.write(('  Output directory already exist.\n'
                              '  Please rename output directory or '
                                'remove existing directory\n\n'))
            sys.exit(1)

        if args.filial < 2:
            sys.stderr.write('  Finial generation must be more than 1.\n')
            sys.exit(1)

        if args.adapter is not None:
            if not args.trim:
                sys.stderr.write(('!!WARNING!! You specified  "--adapter", but you　'
                                  'did not spesify "--trim".\n'
                                  '!!WARNING!! "--trim" was added.\n'))
                sys.stderr.flush()

        if args.trim:
            if args.trim_params is None:
                args.trim_params = '33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75'
        else:
            if args.trim_params is not None:
                sys.stderr.write(('!!WARNING!! You specified  "--trim-params", but you '
                                  'did not spesify "--trim".\n'
                                  '!!WARNING!! "--trim" was added.\n'))
                sys.stderr.flush()

        N_fastq = 0

        for input_name in args.parent:
            n_comma = input_name.count(',')
            if n_comma == 0:
                root, ext = os.path.splitext(input_name)
                if ext != '.bam':
                    sys.stderr.write(('  Please check "{}".\n'
                                      '  The extension of this file is not "bam".\n'
                                      '  If you wanted to specify fastq, please '
                                        'input them as paired-end reads which is separated '
                                        'by comma. e.g. -p fastq1,fastq2\n\n').format(input_name))
                    sys.exit(1)
            elif n_comma == 1:
                fastqs = input_name.split(',')
                for fastq in fastqs:
                    root, ext = os.path.splitext(fastq)
                    if ext == '.bam':
                        sys.stderr.write(('  Please check "{}".\n'
                                          '  The extension must not be "bam".\n'
                                          '  If you wanted to specify bam, please '
                                            'input them separately. e.g. -p bam1 '
                                            '-p bam2\n\n').format(input_name))
                        sys.exit(1)
                N_fastq += 1
            else:
                sys.stderr.write(('  Please check "{}".\n'
                                  '  You specified too much files, or '
                                    'your file name includes ",".\n\n').format(input_name))
                sys.exit(1)

        for input_name in  args.bulk1:
            n_comma = input_name.count(',')
            if n_comma == 0:
                root, ext = os.path.splitext(input_name)
                if ext != '.bam':
                    sys.stderr.write(('  Please check "{}".\n'
                                      '  The extension of this file is not "bam".\n'
                                      '  If you wanted to specify fastq, please '
                                        'input them as paired-end reads which is separated '
                                        'by comma. e.g. -b1 fastq1,fastq2\n\n').format(input_name))
                    sys.exit(1)
            elif n_comma == 1:
                fastqs = input_name.split(',')
                for fastq in fastqs:
                    root, ext = os.path.splitext(fastq)
                    if ext == '.bam':
                        sys.stderr.write(('  Please check "{}".\n'
                                          '  The extension must not be "bam".\n'
                                          '  If you wanted to specify bam, please '
                                            'input them separately. e.g. -b1 bam1 '
                                            '-b1 bam2\n\n').format(input_name))
                        sys.exit(1)
                N_fastq += 1
            else:
                sys.stderr.write(('  Please check "{}".\n'
                                  '  You specified too much files, or '
                                    'your file name includes ",".\n\n').format(input_name))
                sys.exit(1)

        for input_name in  args.bulk2:
            n_comma = input_name.count(',')
            if n_comma == 0:
                root, ext = os.path.splitext(input_name)
                if ext != '.bam':
                    sys.stderr.write(('  Please check "{}".\n'
                                      '  The extension of this file is not "bam".\n'
                                      '  If you wanted to specify fastq, please '
                                        'input them as paired-end reads which is separated '
                                        'by comma. e.g. -b2 fastq1,fastq2\n\n').format(input_name))
                    sys.exit(1)
            elif n_comma == 1:
                fastqs = input_name.split(',')
                for fastq in fastqs:
                    root, ext = os.path.splitext(fastq)
                    if ext == '.bam':
                        sys.stderr.write(('  Please check "{}".\n'
                                          '  The extension must not be "bam".\n'
                                          '  If you wanted to specify bam, please '
                                            'input them separately. e.g. -b2 bam1 '
                                            '-b2 bam2\n\n').format(input_name))
                        sys.exit(1)
                N_fastq += 1
            else:
                sys.stderr.write(('  Please check "{}".\n'
                                  '  You specified too much files, or '
                                    'your file name includes ",".\n\n').format(input_name))
                sys.exit(1)

        if N_fastq == 0 and args.trim:
            sys.stderr.write(('  You can specify "--trim" only when '
                                 'you input fastq.\n\n'))
            sys.exit(1)

        return N_fastq
