import os
import sys
import re
import gzip
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 13


class Plot(object):

    def __init__(self, args):
        self.args = args
        self.out = args.out
        self.vcf = args.vcf
        self.snpEff = args.snpEff
        self.line_colors = args.line_colors.split(',') #SNP-index, p95, and p99
        self.dot_colors = args.dot_colors.split(',') #bulk1, bulk2, and delta
        self.fig_width = args.fig_width
        self.fig_height = args.fig_height
        self.white_space = args.white_space
        self.plot_with_indel = args.indel
        self.snp_index, self.sliding_window = self.read_files()
        N_chr = len(self.sliding_window['CHROM'].unique())

        if N_chr > 50:
            print(('!!WARNING!! Your reference genome has too many contigs (>50). '
                   'Therefore, only the 50 longest contigs will be used for plotting.'), file=sys.stderr)
            N_chr, self.snp_index, self.sliding_window = self.get_50_contigs(N_chr, 
                                                                             self.snp_index, 
                                                                             self.sliding_window)

        self.N_col, self.N_raw = self.set_plot_style(N_chr)
        self.xmax = self.snp_index['POSI'].max()
        self.sliding_window = self.sliding_window.groupby('CHROM')

        try:
            import seaborn as sns
            sns.set_style('whitegrid')
        except ModuleNotFoundError:
            pass

    def read_files(self):
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

        sliding_window = pd.read_csv('{}/sliding_window.tsv'.format(self.out),
                                     sep='\t',
                                     names=['CHROM',
                                            'POSI',
                                            'mean_p99',
                                            'mean_p95',
                                            'mean_bulk1_SNPindex',
                                            'mean_bulk2_SNPindex',
                                            'mean_delta_SNPindex'])

        snp_index['CHROM'] = snp_index['CHROM'].astype('str')
        snp_index['POSI'] = snp_index['POSI']/1000000
        sliding_window['CHROM'] = sliding_window['CHROM'].astype('str')
        sliding_window['POSI'] = sliding_window['POSI']/1000000

        return snp_index, sliding_window

    def read_contig_length(self):
        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')
        contig_length = {}
        for line in vcf:
            if re.match(r'^##contig', line):
                contig = line.split('=')[2].split(',')[0]
                length = int(line.split('=')[3].replace('>', ''))
                contig_length[contig] = length
        return contig_length

    def get_50_contigs(self, N_chr, snp_index, sliding_window):
        contig_length = self.read_contig_length()
        contig_50_names = [k for k, v in sorted(contig_length.items(), key=lambda x:x[1], reverse=True)[:50]]
        N_chr = len(contig_50_names)
        snp_index = snp_index[snp_index['CHROM'].isin(contig_50_names)]
        sliding_window = sliding_window[sliding_window['CHROM'].isin(contig_50_names)]
        return N_chr, snp_index, sliding_window

    def set_plot_style(self, N_chr):
        if N_chr == 1:
            style = (1, 1)
        elif N_chr == 2:
            style = (1, 2)
        else:
            if N_chr % 2 == 0:
                style = (2, int(N_chr/2))
            else:
                style = (2, int(N_chr/2) + 1)
        return style

    def plot_bulk1_SNPindex(self):
        fig = plt.figure(figsize=(self.fig_width*self.N_col, 
                                  self.fig_height*self.N_raw))
        plt.subplots_adjust(hspace=self.white_space)

        for i, (chr_name, chr_sliding_window) in enumerate(self.sliding_window):
            chr_snp_index = self.snp_index[self.snp_index['CHROM']==chr_name]
            if not self.plot_with_indel:
                chr_snp_index = chr_snp_index[chr_snp_index['variant']=='snp']

            ax = fig.add_subplot(self.N_raw, self.N_col, i+1)

            ax.plot(chr_sliding_window['POSI'], 
                    chr_sliding_window['mean_bulk1_SNPindex'], 
                    color=self.line_colors[0], 
                    linewidth=3)

            if self.snpEff is None:
                ax.scatter(chr_snp_index['POSI'], 
                           chr_snp_index['bulk1_SNPindex'], 
                           marker='.', 
                           color=self.dot_colors[0])
            else:
                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODIFIER']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODIFIER']['bulk1_SNPindex'],
                           marker='.',
                           color=self.dot_colors[0])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='LOW']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='LOW']['bulk1_SNPindex'],
                           marker='.',
                           color=self.dot_colors[0])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODERATE']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODERATE']['bulk1_SNPindex'],
                           marker='+',
                           color=self.dot_colors[0])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='HIGH']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='HIGH']['bulk1_SNPindex'],
                           marker='x',
                           color='#772D8B')

            ax.hlines([0.5], 0, self.xmax, linestyles='dashed', color='black')
            ax.set_xlabel('Position (Mbp)')
            ax.set_ylabel('SNP-index 1')
            ax.set_xlim(0, self.xmax)
            ax.set_ylim(0, 1.05)
            ax.set_title(chr_name, fontsize=17)

        plt.savefig('{}/bulk1_SNPindex.{}'.format(self.out, self.args.format))

    def plot_bulk2_SNPindex(self):
        fig = plt.figure(figsize=(self.fig_width*self.N_col, 
                                  self.fig_height*self.N_raw))
        plt.subplots_adjust(hspace=self.white_space)

        for i, (chr_name, chr_sliding_window) in enumerate(self.sliding_window):
            chr_snp_index = self.snp_index[self.snp_index['CHROM']==chr_name]
            if not self.plot_with_indel:
                chr_snp_index = chr_snp_index[chr_snp_index['variant']=='snp']

            ax = fig.add_subplot(self.N_raw, self.N_col, i+1)

            ax.plot(chr_sliding_window['POSI'], 
                    chr_sliding_window['mean_bulk2_SNPindex'], 
                    color=self.line_colors[0], 
                    linewidth=3)

            if self.snpEff is None:
                ax.scatter(chr_snp_index['POSI'], 
                           chr_snp_index['bulk2_SNPindex'], 
                           marker='.', 
                           color=self.dot_colors[1])
            else:
                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODIFIER']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODIFIER']['bulk2_SNPindex'],
                           marker='.',
                           color=self.dot_colors[1])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='LOW']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='LOW']['bulk2_SNPindex'],
                           marker='.',
                           color=self.dot_colors[1])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODERATE']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODERATE']['bulk2_SNPindex'],
                           marker='+',
                           color=self.dot_colors[1])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='HIGH']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='HIGH']['bulk2_SNPindex'],
                           marker='x',
                           color='#772D8B')

            ax.hlines([0.5], 0, self.xmax, linestyles='dashed', color='black')
            ax.set_xlabel('Position (Mbp)')
            ax.set_ylabel('SNP-index 2')
            ax.set_xlim(0, self.xmax)
            ax.set_ylim(0, 1.05)
            ax.set_title(chr_name, fontsize=17)

        plt.savefig('{}/bulk2_SNPindex.{}'.format(self.out, self.args.format))

    def plot_delta_SNPindex(self):
        fig = plt.figure(figsize=(self.fig_width*self.N_col, 
                                  self.fig_height*self.N_raw))
        plt.subplots_adjust(hspace=self.white_space)

        for i, (chr_name, chr_sliding_window) in enumerate(self.sliding_window):
            chr_snp_index = self.snp_index[self.snp_index['CHROM']==chr_name]
            if not self.plot_with_indel:
                chr_snp_index = chr_snp_index[chr_snp_index['variant']=='snp']

            ax = fig.add_subplot(self.N_raw, self.N_col, i+1)

            ax.plot(chr_sliding_window['POSI'], 
                    chr_sliding_window['mean_p99'], 
                    color=self.line_colors[2], 
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'], 
                    chr_sliding_window['mean_p95'], 
                    color=self.line_colors[1], 
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'], 
                    - chr_sliding_window['mean_p99'], 
                    color=self.line_colors[2], 
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'], 
                    - chr_sliding_window['mean_p95'],
                    color=self.line_colors[1], 
                    linewidth=3)

            ax.plot(chr_sliding_window['POSI'], 
                    chr_sliding_window['mean_delta_SNPindex'], 
                    color=self.line_colors[0], 
                    linewidth=3)

            if self.snpEff is None:
                ax.scatter(chr_snp_index['POSI'], 
                           chr_snp_index['delta_SNPindex'], 
                           marker='.', 
                           color=self.dot_colors[2])
            else:
                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODIFIER']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODIFIER']['delta_SNPindex'],
                           marker='.',
                           color=self.dot_colors[2])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='LOW']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='LOW']['delta_SNPindex'],
                           marker='.',
                           color=self.dot_colors[2])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='MODERATE']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='MODERATE']['delta_SNPindex'],
                           marker='+',
                           color=self.dot_colors[2])

                ax.scatter(chr_snp_index[chr_snp_index['impact']=='HIGH']['POSI'],
                           chr_snp_index[chr_snp_index['impact']=='HIGH']['delta_SNPindex'],
                           marker='x',
                           color='#772D8B')

            ax.hlines([0], 0, self.xmax, linestyles='dashed', color='black')
            ax.set_xlabel('Position (Mbp)')
            ax.set_ylabel('$\Delta$SNP-index')
            ax.set_xlim(0, self.xmax)
            ax.set_ylim(-1.05, 1.05)
            ax.set_title(chr_name, fontsize=17)

        plt.savefig('{}/delta_SNPindex.{}'.format(self.out, self.args.format))


    def run(self):
        self.plot_bulk1_SNPindex()
        self.plot_bulk2_SNPindex()
        self.plot_delta_SNPindex()
