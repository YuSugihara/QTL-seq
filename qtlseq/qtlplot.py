#!/usr/bin/env python3

from qtlseq.params import Params


pm = Params('qtlplot')
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
    import pandas as pd
    from qtlseq.vcf2index import Vcf2Index
    from qtlseq.plot import Plot
    from qtlseq.utils import time_stamp, check_required_commands
except ModuleNotFoundError as exc:
    sys.stderr.write(('Error: Missing Python dependency "{}".\n'
                      'Please reinstall QTL-seq with "pip install -e ." or use the conda package.\n\n').format(exc.name))
    sys.exit(1)


class QTLplot(object):

    def __init__(self, args):
        args = pm.check_max_threads(args)
        pm.check_qtlplot_args(args)
        self.check_dependencies(args)
        self.out = args.out
        self.vcf = args.vcf
        self.snpEff = args.snpEff
        self.args = args

    def check_dependencies(self, args):
        required = []
        if args.snpEff is not None:
            required.append('snpEff')
        if len(required) > 0:
            check_required_commands(required, 'QTL-plot')

    def check_outdir(self):
        if os.path.isdir(self.out):
            print(time_stamp(),
                  'Output directory already exists.'.format(self.out),
                  flush=True)
        else:
            try:
                os.mkdir(self.out)
            except OSError as exc:
                sys.stderr.write(('{} Error: Failed to create output directory "{}": {}\n'
                                  '{} Please check permissions and path.\n\n').format(time_stamp(),
                                                                                      self.out,
                                                                                      exc,
                                                                                      time_stamp()))
                sys.exit(1)

    def run_snpEff(self):
        summary_file = '{}/snpEff_summary.html'.format(self.out)
        snpeff_vcf = '{}/qtlseq.snpEff.vcf'.format(self.out)
        snpeff_log = '{}/snpEff.log'.format(self.out)
        cmd = ['snpEff',
               'ann',
               '-Xmx4g',
               '-s',
               summary_file,
               self.snpEff,
               self.vcf]

        print(time_stamp(), 'Running SnpEff annotation.', flush=True)
        try:
            with open(snpeff_vcf, 'w') as out_vcf, open(snpeff_log, 'w') as out_log:
                sbp.run(cmd, stdout=out_vcf, stderr=out_log, check=True)
        except FileNotFoundError:
            sys.stderr.write(('{} Error: "snpEff" command was not found.\n'
                              '{} Please install SnpEff or remove "-e/--snpEff".\n\n').format(time_stamp(), time_stamp()))
            sys.exit(1)
        except sbp.CalledProcessError:
            sys.stderr.write(('{} Error: SnpEff annotation failed.\n'
                              '{} Please check "{}".\n\n').format(time_stamp(), time_stamp(), snpeff_log))
            sys.exit(1)

        print(time_stamp(), 'SnpEff annotation completed successfully.', flush=True)

        self.args.vcf = snpeff_vcf

    def has_non_empty_result(self):
        snp_index = '{}/snp_index.tsv'.format(self.out)
        sliding_window = '{}/sliding_window.tsv'.format(self.out)
        return os.path.isfile(snp_index) and os.path.getsize(snp_index) > 0 and \
               os.path.isfile(sliding_window) and os.path.getsize(sliding_window) > 0

    def get_outlier_SNPindex(self):
        if self.snpEff is None:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                sep='\t',
                                names=['CHROM',
                                       'POSI',
                                       'variant',
                                       'bulk1_depth',
                                       'bulk2_depth',
                                       'p99',
                                       'p95',
                                       'bulk1_SNPindex',
                                       'bulk2_SNPindex',
                                       'delta_SNPindex'])
        else:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                sep='\t',
                                names=['CHROM',
                                       'POSI',
                                       'variant',
                                       'impact',
                                       'bulk1_depth',
                                       'bulk2_depth',
                                       'p99',
                                       'p95',
                                       'bulk1_SNPindex',
                                       'bulk2_SNPindex',
                                       'delta_SNPindex'])

        snp_index[abs(snp_index['p99']) <= \
                  abs(snp_index['delta_SNPindex'])].to_csv('{}/snp_index.p99.tsv'.format(self.out),
                                                           sep='\t',
                                                           index=False)

        snp_index[abs(snp_index['p95']) <= \
                  abs(snp_index['delta_SNPindex'])].to_csv('{}/snp_index.p95.tsv'.format(self.out),
                                                           sep='\t',
                                                           index=False)

    def get_outlier_windows(self):
        sliding_window = pd.read_csv('{}/sliding_window.tsv'.format(self.out),
                                     sep='\t',
                                     names=['CHROM',
                                            'POSI',
                                            'mean_p99',
                                            'mean_p95',
                                            'mean_SNPindex1',
                                            'mean_SNPindex2',
                                            'mean_SNPindex'])

        sliding_window[abs(sliding_window['mean_p99']) <= \
                       abs(sliding_window['mean_SNPindex'])].to_csv('{}/sliding_window.p99.tsv'.format(self.out),
                                                                    sep='\t',
                                                                    index=False)

        sliding_window[abs(sliding_window['mean_p95']) <= \
                       abs(sliding_window['mean_SNPindex'])].to_csv('{}/sliding_window.p95.tsv'.format(self.out),
                                                                    sep='\t',
                                                                    index=False)

    def make_igv_file(self):
        if self.snpEff is None:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                sep='\t',
                                names=['CHROM',
                                       'POSI',
                                       'variant',
                                       'bulk1_depth',
                                       'bulk2_depth',
                                       'p99',
                                       'p95',
                                       'bulk1_SNPindex',
                                       'bulk2_SNPindex',
                                       'delta_SNPindex'])
        else:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                sep='\t',
                                names=['CHROM',
                                       'POSI',
                                       'variant',
                                       'impact',
                                       'bulk1_depth',
                                       'bulk2_depth',
                                       'p99',
                                       'p95',
                                       'bulk1_SNPindex',
                                       'bulk2_SNPindex',
                                       'delta_SNPindex'])

        snp_index['Start'] = snp_index['POSI'] - 1
        snp_index['End'] = snp_index['POSI']
        snp_index['Feature'] = snp_index['CHROM'].astype('str') \
                               + ':' + snp_index['POSI'].astype('str') \
                               + ':' + snp_index['variant'].astype('str') \
                               + ':' + snp_index['p99'].astype('str') \
                               + ':' + snp_index['p95'].astype('str')

        snp_index = snp_index[['CHROM',
                               'Start',
                               'End',
                               'Feature',
                               'delta_SNPindex']]

        snp_index.to_csv('{}/snp_index.igv'.format(self.out),
                         sep='\t',
                         index=False)

        sliding_window = pd.read_csv('{}/sliding_window.tsv'.format(self.out),
                                     sep='\t',
                                     names=['CHROM',
                                            'POSI',
                                            'mean_p99',
                                            'mean_p95',
                                            'mean_SNPindex1',
                                            'mean_SNPindex2',
                                            'mean_SNPindex'])

        sliding_window['Start'] = sliding_window['POSI'] - 1
        sliding_window['End'] = sliding_window['POSI']
        sliding_window['Feature'] = sliding_window['CHROM'].astype('str') \
                                    + ':' + sliding_window['POSI'].astype('str') \
                                    + ':' + sliding_window['mean_p99'].astype('str') \
                                    + ':' + sliding_window['mean_p95'].astype('str')

        sliding_window = sliding_window[['CHROM',
                                         'Start',
                                         'End',
                                         'Feature',
                                         'mean_SNPindex']]

        sliding_window.to_csv('{}/sliding_window.igv'.format(self.out),
                              sep='\t',
                              index=False)

    def run(self):
        self.check_outdir()

        if self.snpEff is not None:
            self.run_snpEff()

        v2i = Vcf2Index(self.args)
        v2i.run()

        if not self.has_non_empty_result():
            print(time_stamp(),
                  'No variants passed filtering; skipping plot generation.',
                  flush=True)
            return

        print(time_stamp(), 'Generating plot...', flush=True)
        try:
            pt = Plot(self.args)
            pt.run()
        except Exception as exc:
            sys.stderr.write(('{} Error: Plot generation failed: {}\n'
                              '{} Please check "{}" and "{}".\n'
                              '{} You can also try changing "--fig-width", "--fig-height", or "--format".\n\n').format(time_stamp(),
                                                                                                                         exc,
                                                                                                                         time_stamp(),
                                                                                                                         '{}/snp_index.tsv'.format(self.out),
                                                                                                                         '{}/sliding_window.tsv'.format(self.out),
                                                                                                                         time_stamp()))
            sys.exit(1)

        self.get_outlier_SNPindex()
        self.get_outlier_windows()

        if self.args.igv:
            self.make_igv_file()

def main():
    print(time_stamp(), 'Starting QTL-plot.', flush=True)
    try:
        qp = QTLplot(args)
        qp.run()
    except KeyboardInterrupt:
        sys.stderr.write('\n{} Error: Interrupted by user.\n\n'.format(time_stamp()))
        sys.exit(1)
    except SystemExit:
        raise
    except Exception as exc:
        sys.stderr.write(('{} Error: Unexpected failure in QTL-plot: {}\n'
                          '{} Please check your input VCF and output directory.\n\n').format(time_stamp(),
                                                                                              exc,
                                                                                              time_stamp()))
        sys.exit(1)
    print(time_stamp(), 'QTL-plot completed successfully.', flush=True)

if __name__ == '__main__':
    main()
