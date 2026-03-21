
import os
import re
import sys
import gzip
from multiprocessing import Pool
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
        self.N_bulk1 = args.N_bulk1
        self.N_bulk2 = args.N_bulk2
        self.N_replicates = args.N_rep
        self.min_SNPindex = args.min_SNPindex
        self.snp_index = '{}/snp_index.tsv'.format(self.out)
        self.sf = SnpFilt(args)
        self.args = args

        if self.snpEff is not None:
            self.ANN_re = re.compile(';ANN=(.*);*')

        self.p99_index = int(0.99*self.N_replicates) - 1
        self.p95_index = int(0.95*self.N_replicates) - 1

    def get_field(self):
        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')
        GT_pos = None
        AD_pos = None
        ADF_pos = None
        ADR_pos = None
        sample_count = None
        has_column_header = False
        for line in vcf:
            if line.startswith('#CHROM'):
                has_column_header = True
                columns = line.rstrip('\n').split('\t')
                if len(columns) < 12:
                    sys.stderr.write(('{} Error: Your VCF must include at least 3 samples (parent, bulk1, bulk2).\n'
                                      '{} Found {} sample column(s).\n'
                                      '{} Please check VCF column order: parent, bulk1, bulk2.\n\n').format(time_stamp(),
                                                                                                              time_stamp(),
                                                                                                              max(len(columns) - 9, 0),
                                                                                                              time_stamp()))
                    sys.exit(1)
                sample_count = len(columns) - 9
                if sample_count > 3:
                    sys.stderr.write(('{} Warning: Your VCF has {} sample columns.\n'
                                      '{} QTL-seq will use only the first three samples as parent, bulk1, bulk2.\n').format(time_stamp(),
                                                                                                                              sample_count,
                                                                                                                              time_stamp()))
            if re.match(r'[^#]', line):
                fields = line.split()[8].split(':')
                try:
                    GT_pos = fields.index('GT')
                except ValueError:
                    sys.stderr.write(('{} Error: No GT (Genotype) field found in your VCF file!\n').format(time_stamp()))
                    sys.exit(1)
                try:
                    AD_pos = fields.index('AD')
                except ValueError:
                    sys.stderr.write(('{} Error: No AD (Allelic Depth) field found in your VCF file!\n').format(time_stamp()))
                    sys.exit(1)

                if 'ADF' in fields and 'ADR' in fields:
                    ADF_pos = fields.index('ADF')
                    ADR_pos = fields.index('ADR')
                else:
                    ADF_pos = None
                    ADR_pos = None
                    sys.stderr.write(('{} Warning: No ADF (Allelic Depth Forward) or ADR (Allelic Depth Reverse) field found in your VCF file.\n').format(time_stamp()))
                    sys.stderr.write(('{} Warning: Strand bias filtering will be skipped.\n').format(time_stamp()))
                break
        vcf.close()
        if not has_column_header:
            sys.stderr.write(('{} Error: No #CHROM header line found in "{}".\n'
                              '{} Please check whether this is a valid VCF file.\n\n').format(time_stamp(),
                                                                                               self.vcf,
                                                                                               time_stamp()))
            sys.exit(1)
        if GT_pos is None or AD_pos is None:
            sys.stderr.write(('{} Error: No variant records were found in "{}".\n').format(time_stamp(), self.vcf))
            sys.exit(1)
        return GT_pos, AD_pos, ADF_pos, ADR_pos

    def get_variant_impact(self, annotation):
        ann_matches = self.ANN_re.findall(annotation)
        if len(ann_matches) == 0:
            return 'MODIFIER'

        ANN = ann_matches[0]
        genes = ANN.split(',')
        # The impact of the variant is determined by the most severe impact
        # of the genes that the variant is located in.
        # The impact order is HIGH > MODERATE > LOW > MODIFIER
        # The impacts of variants are oredered by severity in the ANN field.
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

    @lru_cache(maxsize=None)
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
        p99 = replicates[self.p99_index]
        p95 = replicates[self.p95_index]
        return p99, p95

    def calculate_SNPindex_sub(self, line):
        if re.match(r'[^#]', line):
            cols = line.split('\t')
            if len(cols) < 12:
                return None

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

            if ADF_pos is not None and ADR_pos is not None:
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
                if self.snpEff is None:
                    p99, p95 = self.Fn_simulation(depth1, depth2)
                    return ('{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t'
                            '{:.4f}\t{:.4f}\t{:.4f}\n').format(CHR,
                                                               POS,
                                                               variant,
                                                               record['bulk1_depth'],
                                                               record['bulk2_depth'],
                                                               p99,
                                                               p95,
                                                               record['bulk1_SNPindex'],
                                                               record['bulk2_SNPindex'],
                                                               record['delta_SNPindex'])
                else:
                    p99, p95 = self.Fn_simulation(depth1, depth2)
                    impact = self.get_variant_impact(annotation)
                    return ('{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t'
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
                                                               record['delta_SNPindex'])
        return None

    def table_sort(self):
        temp_path = '{}/snp_index.tsv.temp'.format(self.out)
        out_path = '{}/snp_index.tsv'.format(self.out)
        if not os.path.exists(temp_path) or os.path.getsize(temp_path) == 0:
            open(out_path, 'w').close()
            if os.path.exists(temp_path):
                os.remove(temp_path)
            return False

        if self.snpEff is None:
            names = ['CHROM',
                     'POSI',
                     'variant',
                     'bulk1_depth',
                     'bulk2_depth',
                     'p99',
                     'p95',
                     'bulk1_SNPindex',
                     'bulk2_SNPindex',
                     'delta_SNPindex']
        else:
            names = ['CHROM',
                     'POSI',
                     'variant',
                     'impact',
                     'bulk1_depth',
                     'bulk2_depth',
                     'p99',
                     'p95',
                     'bulk1_SNPindex',
                     'bulk2_SNPindex',
                     'delta_SNPindex']

        try:
            snp_index = pd.read_csv(temp_path, sep='\t', names=names)
        except pd.errors.EmptyDataError:
            open(out_path, 'w').close()
            os.remove(temp_path)
            return False

        snp_index = snp_index.sort_values(by=['CHROM', 'POSI'])

        snp_index.to_csv(out_path,
                         sep='\t',
                         index=False,
                         header=False)

        os.remove(temp_path)
        return True

    def calculate_SNPindex(self):
        temp_path = '{}/snp_index.tsv.temp'.format(self.out)
        if os.path.exists(temp_path):
            os.remove(temp_path)

        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')

        try:
            with open(temp_path, 'w') as snp_index_temp:
                with Pool(self.args.threads) as p:
                    for record in p.imap(self.calculate_SNPindex_sub, vcf, chunksize=200):
                        if record is not None:
                            snp_index_temp.write(record)
        except Exception as exc:
            vcf.close()
            sys.stderr.write(('{} Error: Failed to parse VCF records while calculating SNP-index: {}\n'
                              '{} Please verify FORMAT/sample fields (GT, AD) and sample order (parent, bulk1, bulk2).\n\n').format(time_stamp(),
                                                                                                                                      exc,
                                                                                                                                      time_stamp()))
            sys.exit(1)

        vcf.close()
        return self.table_sort()

    def run(self):
        print(time_stamp(), 'Calculating SNP-index.', flush=True)
        self.field_pos = self.get_field()
        has_result = self.calculate_SNPindex()
        print(time_stamp(), 'SNP-index calculation completed successfully.', flush=True)

        print(time_stamp(), 'Smoothing SNP-index.', flush=True)
        sm = Smooth(self.args)
        sm.run()
        print(time_stamp(), 'SNP-index smoothing completed successfully.', flush=True)
        if not has_result:
            print(time_stamp(),
                  'Warning: No variants passed filtering. Empty result files were generated.',
                  file=sys.stderr,
                  flush=True)
