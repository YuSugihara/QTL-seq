
import os
import re
import sys
import gzip
import threading
from multiprocessing import Pool
from multiprocessing import Manager
lock = threading.Lock()
cache = Manager().dict()
import numpy as np
import pandas as pd
from functools import lru_cache
from qtlseq.snpfilt import SnpFilt
from qtlseq.smooth import Smooth
from qtlseq.utils import time_stamp


class Vcf2Index(object):

    def __init__(self, args):
        self.out = args.out
        self.vcf = args.vcf
        self.snpEff = args.snpEff
        self.filial = args.filial
        self.species = args.species
        self.N_bulk1 = args.N_bulk1
        self.N_bulk2 = args.N_bulk2
        self.N_replicates = args.N_rep
        self.min_SNPindex = args.min_SNPindex
        self.snp_index = '{}/snp_index.tsv'.format(self.out)
        self.sf = SnpFilt(args)
        self.args = args

        if self.snpEff is not None:
            self.ANN_re = re.compile(';ANN=(.*);*')

        if self.species is not None:
            self.u99, self.u95 = self.correct_threshold()

    def correct_threshold(self):
        if self.species == 'Arabidopsis':
            if self.filial == 2:
                u99 = 3.82
                u95 = 3.41
            elif self.filial == 3:
                u99 = 3.91
                u95 = 3.50
            elif self.filial == 4:
                u99 = 3.94
                u95 = 3.54
            else:
                u99 = 3.97
                u95 = 3.57

        elif self.species == 'Cucumber':
            if self.filial == 2:
                u99 = 4.02
                u95 = 3.62
            elif self.filial == 3:
                u99 = 4.10
                u95 = 3.72
            elif self.filial == 4:
                u99 = 4.14
                u95 = 3.75
            else:
                u99 = 4.16
                u95 = 3.78

        elif self.species == 'Maize':
            if self.filial == 2:
                u99 = 4.11
                u95 = 3.72
            elif self.filial == 3:
                u99 = 4.19
                u95 = 3.81
            elif self.filial == 4:
                u99 = 4.23
                u95 = 3.85
            else:
                u99 = 4.25
                u95 = 3.87

        elif self.species == 'Rapeseed':
            if self.filial == 2:
                u99 = 4.16
                u95 = 3.78
            elif self.filial == 3:
                u99 = 4.24
                u95 = 3.87
            elif self.filial == 4:
                u99 = 4.27
                u95 = 3.90
            else:
                u99 = 4.30
                u95 = 3.93

        elif self.species == 'Rice':
            if self.filial == 2:
                u99 = 4.05
                u95 = 3.65
            elif self.filial == 3:
                u99 = 4.13
                u95 = 3.74
            elif self.filial == 4:
                u99 = 4.16
                u95 = 3.78
            else:
                u99 = 4.19
                u95 = 3.80

        elif self.species == 'Tobacco':
            if self.filial == 2:
                u99 = 4.22
                u95 = 3.84
            elif self.filial == 3:
                u99 = 4.30
                u95 = 3.92
            elif self.filial == 4:
                u99 = 4.33
                u95 = 3.96
            else:
                u99 = 4.35
                u95 = 3.98

        elif self.species == 'Tomato':
            if self.filial == 2:
                u99 = 4.04
                u95 = 3.65
            elif self.filial == 3:
                u99 = 4.12
                u95 = 3.73
            elif self.filial == 4:
                u99 = 4.15
                u95 = 3.77
            else:
                u99 = 4.18
                u95 = 3.80

        elif self.species == 'Wheat':
            if self.filial == 2:
                u99 = 4.21
                u95 = 3.83
            elif self.filial == 3:
                u99 = 4.29
                u95 = 3.92
            elif self.filial == 4:
                u99 = 4.32
                u95 = 3.95
            else:
                u99 = 4.35
                u95 = 3.98

        elif self.species == 'Yeast':
            u99 = 4.31
            u95 = 3.93
            if self.filial >= 2:
                print('!!WARNING!! Filial generation must be 2 in yeast.', file=sys.stderr)

        else:
            print('You specified not supported species.', file=sys.stderr)
            sys.exit(1)

        return u99, u95


    def get_field(self):
        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')
        for line in vcf:
            if re.match(r'[^#]', line):
                fields = line.split()[8].split(':')
                try:
                    GT_pos = fields.index('GT')
                except ValueError:
                    sys.stderr.write(('{} No GT field'
                                      ' in your VCF!!\n').format(time_stamp()))
                    sys.exit(1)
                try:
                    AD_pos = fields.index('AD')
                except ValueError:
                    sys.stderr.write(('{} No AD field'
                                      ' in your VCF!!\n').format(time_stamp()))
                    sys.exit(1)

                if 'ADF' in fields and 'ADR' in fields:
                    ADF_pos = fields.index('ADF')
                    ADR_pos = fields.index('ADR')
                else:
                    ADF_pos = None
                    ADR_pos = None
                    sys.stderr.write(('{} no ADF or ADR field'
                                      ' in your VCF.\n').format(time_stamp()))
                    sys.stderr.write(('{} strand bias filter '
                                      'will be skipped.\n').format(time_stamp()))
                break
        vcf.close()
        return GT_pos, AD_pos, ADF_pos, ADR_pos

    def get_variant_impact(self, annotation):
        ANN = self.ANN_re.findall(annotation)[0]
        genes = ANN.split(',')
        impacts = [gene.split('|')[2] for gene in genes]
        if 'HIGH' in impacts:
            impact = 'HIGH'
        elif 'MODERATE' in impacts:
            impact = 'MODERATE'
        elif 'LOW' in impacts:
            impact = 'LOW'
        else:
            impact = 'MODIFIER'
        return impact

    def check_variant_type(self, REF, ALT):
        if (len(REF) == 1) and (len(ALT) == 1):
            variant = 'snp'
        else:
            variant = 'indel'
        return variant

    def check_depth(self, bulk1_depth, bulk2_depth):
        if bulk1_depth < bulk2_depth:
            depth1 = bulk1_depth
            depth2 = bulk2_depth
        else:
            depth1 = bulk2_depth
            depth2 = bulk1_depth
        return depth1, depth2

    def Fn_simulation(self, depth1, depth2):
        GT = [0, 1]
        replicates = []
        n = 0
        while n < self.N_replicates:
            bulk1_Fn_GT = np.ones(self.N_bulk1)
            bulk2_Fn_GT = np.ones(self.N_bulk2)
            for _ in range(self.args.filial - 1):
                bulk1_random_gt = np.random.choice([0, 1, 1, 2],
                                                   bulk1_Fn_GT.shape,
                                                   replace=True)

                bulk2_random_gt = np.random.choice([0, 1, 1, 2],
                                                   bulk2_Fn_GT.shape,
                                                   replace=True)

                bulk1_Fn_GT = np.where(bulk1_Fn_GT==1,
                                       bulk1_random_gt,
                                       bulk1_Fn_GT)

                bulk2_Fn_GT = np.where(bulk2_Fn_GT==1,
                                       bulk2_random_gt,
                                       bulk2_Fn_GT)

            bulk1_freq = sum(bulk1_Fn_GT)/(2*self.N_bulk1)
            bulk2_freq = sum(bulk2_Fn_GT)/(2*self.N_bulk2)

            bulk1_AD = np.random.choice([0, 1],
                                        depth1,
                                        p=[1 - bulk1_freq, bulk1_freq],
                                        replace=True)

            bulk2_AD = np.random.choice([0, 1],
                                        depth2,
                                        p=[1 - bulk2_freq, bulk2_freq],
                                        replace=True)

            bulk1_SNPindex = bulk1_AD.sum()/depth1
            bulk2_SNPindex = bulk2_AD.sum()/depth2

            if bulk1_SNPindex < self.min_SNPindex and \
               bulk2_SNPindex < self.min_SNPindex:
                continue
            else:
                delta_SNPindex = bulk2_SNPindex - bulk1_SNPindex
                replicates.append(abs(delta_SNPindex))
                n += 1

        replicates.sort()
        p99 = replicates[int(0.99*self.N_replicates) - 1]
        p95 = replicates[int(0.95*self.N_replicates) - 1]
        return p99, p95

    def calculate_SNPindex_sub(self, line):
        if re.match(r'[^#]', line):
            cols = line.split('\t')
            CHR = cols[0]
            POS = cols[1]
            REF = cols[3]
            ALT = cols[4]
            annotation = cols[7]
            GT_pos = self.field_pos[0]
            AD_pos = self.field_pos[1]
            ADF_pos = self.field_pos[2]
            ADR_pos = self.field_pos[3]

            parent_GT = cols[9].split(':')[GT_pos]
            parent_AD = cols[9].split(':')[AD_pos]
            bulk1_AD = cols[10].split(':')[AD_pos]
            bulk2_AD = cols[11].split(':')[AD_pos]

            if ADF_pos != None and ADR_pos != None:
                parent_ADF = cols[9].split(':')[ADF_pos]
                parent_ADR = cols[9].split(':')[ADR_pos]
                ADFR = (parent_ADF, parent_ADR)
            else:
                ADFR = None

            record = self.sf.filt(parent_GT, parent_AD, bulk1_AD, bulk2_AD, ADFR)
            if record['type'] == 'keep':
                variant = self.check_variant_type(REF, ALT)
                depth1, depth2 = self.check_depth(record['bulk1_depth'],
                                                  record['bulk2_depth'])

                if (depth1, depth2) in cache:
                    p99, p95 = cache[(depth1, depth2)]
                else:
                    if self.species is None:
                        p99, p95 = self.Fn_simulation(depth1, depth2)
                    else:
                        var_fix = 0.25*(depth1 + depth2) / \
                                       (depth1*depth2)
                        p99 = self.u99*var_fix
                        p95 = self.u95*var_fix

                    cache[(depth1, depth2)] = (p99, p95)

                lock.acquire()
                snp_index = open(self.snp_index + ".temp", 'a')
                if self.snpEff is None:
                    snp_index.write(('{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t'
                                     '{:.4f}\t{:.4f}\t{:.4f}\n').format(CHR,
                                                                        POS,
                                                                        variant,
                                                                        record['bulk1_depth'],
                                                                        record['bulk2_depth'],
                                                                        p99,
                                                                        p95,
                                                                        record['bulk1_SNPindex'],
                                                                        record['bulk2_SNPindex'],
                                                                        record['delta_SNPindex']))
                else:
                    impact = self.get_variant_impact(annotation)
                    snp_index.write(('{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t'
                                     '{:.4f}\t{:.4f}\t{:.4f}\n').format(CHR,
                                                                        POS,
                                                                        variant,
                                                                        impact,
                                                                        record['bulk1_depth'],
                                                                        record['bulk2_depth'],
                                                                        p99,
                                                                        p95,
                                                                        record['bulk1_SNPindex'],
                                                                        record['bulk2_SNPindex'],
                                                                        record['delta_SNPindex']))
                snp_index.close()
                lock.release()

    def table_sort(self):
        if self.snpEff is None:
            snp_index = pd.read_csv('{}/snp_index.tsv.temp'.format(self.out),
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
            snp_index = pd.read_csv('{}/snp_index.tsv.temp'.format(self.out),
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

        snp_index = snp_index.sort_values(by=['CHROM', 'POSI'])

        snp_index.to_csv('{}/snp_index.tsv'.format(self.out),
                         sep='\t',
                         index=False,
                         header=False)

        os.remove('{}/snp_index.tsv.temp'.format(self.out))

    def calculate_SNPindex(self):
        if os.path.exists('{}/snp_index.tsv.temp'.format(self.out)):
            os.remove('{}/snp_index.tsv.temp'.format(self.out))

        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')

        p = Pool(self.args.threads)
        p.map(self.calculate_SNPindex_sub, vcf)
        p.close()

        vcf.close()

        self.table_sort()

    def run(self):
        print(time_stamp(), 'start to calculate SNP-index.', flush=True)
        self.field_pos = self.get_field()
        self.calculate_SNPindex()
        print(time_stamp(), 'SNP-index successfully finished.', flush=True)

        print(time_stamp(), 'start to smooth SNP-index.', flush=True)
        sm = Smooth(self.args)
        sm.run()
        print(time_stamp(), 'smoothing process successfully finished.', flush=True)
