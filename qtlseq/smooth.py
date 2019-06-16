
import pandas as pd


class Smooth(object):

    def __init__(self, args):
        self.out = args.out
        self.snpEff = args.snpEff
        self.window_size = args.window*1000
        self.step_size = args.step*1000
        self.include_indel = not args.indel

    def calc_sliding_window(self, chr_name, chrom):
        cols = chrom.loc[:, ['POSI', 
                             'p99', 
                             'p95', 
                             'bulk1_SNPindex', 
                             'bulk2_SNPindex', 
                             'delta_SNPindex']].values

        max_posi = int(cols[:, 0].max())
        win_edge = self.window_size
        if max_posi <= win_edge:
            window = cols[(cols[:, 0] < win_edge) & (cols[:, 0] >= (win_edge - self.window_size))]
            N_SNPs_in_window = window.shape[0]
            if N_SNPs_in_window == 0:
                self.sliding_window.write('{}\t{}\t{}\t{}\t'
                                          '{}\t{}\t{}\n'.format(chr_name,
                                                                int(max_posi/2),
                                                                'nan',
                                                                'nan',
                                                                'nan',
                                                                'nan',
                                                                'nan'))
            else:
                mean_p99 = window[:, 1].mean()
                mean_p95 = window[:, 2].mean()
                mean_bulk1_SNPindex = window[:, 3].mean()
                mean_bulk2_SNPindex = window[:, 4].mean()
                mean_delta_SNPindex = window[:, 5].mean()
                self.sliding_window.write('{}\t{}\t{:.4f}\t{:.4f}\t'
                                          '{:.4f}\t{:.4f}\t{:.4f}\n'.format(chr_name,
                                                                            int(max_posi/2),
                                                                            mean_p99,
                                                                            mean_p95,
                                                                            mean_bulk1_SNPindex,
                                                                            mean_bulk2_SNPindex,
                                                                            mean_delta_SNPindex))
        else:
            win_edges = [n for n in range(self.window_size, max_posi+1, self.step_size)]
            for win_edge in win_edges:
                window = cols[(cols[:, 0] < win_edge) & (cols[:, 0] >= (win_edge - self.window_size))]
                N_SNPs_in_window = window.shape[0]
                if N_SNPs_in_window == 0:
                    self.sliding_window.write('{}\t{}\t{}\t{}\t'
                                              '{}\t{}\t{}\n'.format(chr_name,
                                                                    int(win_edge - self.window_size/2),
                                                                    'nan',
                                                                    'nan',
                                                                    'nan',
                                                                    'nan',
                                                                    'nan'))
                else:
                    mean_p99 = window[:, 1].mean()
                    mean_p95 = window[:, 2].mean()
                    mean_bulk1_SNPindex = window[:, 3].mean()
                    mean_bulk2_SNPindex = window[:, 4].mean()
                    mean_delta_SNPindex = window[:, 5].mean()
                    self.sliding_window.write('{}\t{}\t{:.4f}\t{:.4f}\t'
                                              '{:.4f}\t{:.4f}\t{:.4f}\n'.format(chr_name,
                                                                                int(win_edge - self.window_size/2),
                                                                                mean_p99,
                                                                                mean_p95,
                                                                                mean_bulk1_SNPindex,
                                                                                mean_bulk2_SNPindex,
                                                                                mean_delta_SNPindex))

    def run(self):
        self.sliding_window = open('{}/sliding_window.tsv'.format(self.out), 'w')
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

        if self.include_indel:
            snp_index = snp_index[snp_index['variant'] == 'snp']

        snp_index['CHROM'] = snp_index['CHROM'].astype('category')
        grouped = snp_index.groupby('CHROM')
        for chr_name, chrom in grouped:
            self.calc_sliding_window(chr_name, chrom)
        self.sliding_window.close()
