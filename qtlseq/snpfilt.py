

class SnpFilt(object):

    def __init__(self, args):
        self.min_SNPindex = args.min_SNPindex
        self.maxDP = args.max_depth
        self.minDP = args.min_depth
        self.strand_bias = args.strand_bias

    def filt_parent_GT(self, parent_GT, bulk1_AD, bulk2_AD):
        record = {}
        # only use homo in cultivar
        if parent_GT in ['0/0', '0|0', '1/1', '1|1']:
            bulk1_ADs = bulk1_AD.split(',')
            bulk2_ADs = bulk2_AD.split(',')
            # check biallele or multi-allele
            if len(bulk1_ADs) == 2:
                # filter missing
                if (not '.' in bulk1_ADs) and (not '.' in bulk2_ADs):
                    # check whether REF homo or ALT homo.
                    if parent_GT in ['0/0', '0|0']:
                        record['bulk1_ref_AD'] = int(bulk1_ADs[0])
                        record['bulk1_alt_AD'] = int(bulk1_ADs[1])
                        record['bulk2_ref_AD'] = int(bulk2_ADs[0])
                        record['bulk2_alt_AD'] = int(bulk2_ADs[1])
                        # if depths of REF or ALT are zero in both bulks,
                        # it will be discarded because SNP-index will be zero.
                        if (record['bulk1_ref_AD'] != 0 or record['bulk2_ref_AD'] != 0) and \
                           (record['bulk1_alt_AD'] != 0 or record['bulk2_alt_AD'] != 0):
                            record['parent_GT'] = '0/0'
                            record['type'] = 'keep'
                        else:
                            record['type'] = 'discard'
                    else:
                        record['bulk1_ref_AD'] = int(bulk1_ADs[1])
                        record['bulk1_alt_AD'] = int(bulk1_ADs[0])
                        record['bulk2_ref_AD'] = int(bulk2_ADs[1])
                        record['bulk2_alt_AD'] = int(bulk2_ADs[0])
                        # if depth of REF is zero in bulk,
                        # it will be discarded because SNP-index will be zero.
                        if (record['bulk1_ref_AD'] != 0 or record['bulk2_ref_AD'] != 0) and \
                           (record['bulk1_alt_AD'] != 0 or record['bulk2_alt_AD'] != 0):
                            record['parent_GT'] = '1/1'
                            record['type'] = 'keep'
                        else:
                            record['type'] = 'discard'
                else:
                    record['type'] = 'discard'
            # check ALT homo in cultivar
            elif len(bulk1_ADs) == 3:
                # filter missing
                if (not '.' in bulk1_ADs) and (not '.' in bulk2_ADs):
                    if parent_GT in ['1/1', '1|1']:
                        if (bulk1_ADs[0] == '0') and (bulk2_ADs[0] == '0'):
                            record['bulk1_ref_AD'] = int(bulk1_ADs[1])
                            record['bulk1_alt_AD'] = int(bulk1_ADs[2])
                            record['bulk2_ref_AD'] = int(bulk2_ADs[1])
                            record['bulk2_alt_AD'] = int(bulk2_ADs[2])
                            record['parent_GT'] = '1/1'
                            record['type'] = 'keep'
                        else:
                            record['type'] = 'discard'
                    else:
                        record['type'] = 'discard'
                else:
                    record['type'] = 'discard'
            else:
                record['type'] = 'discard'
        else:
            record['type'] = 'discard'
        return record

    def filt_depth(self, record, parent_AD):
        record['parent_depth'] = sum([int(AD) for AD in parent_AD.split(',')])
        record['bulk1_depth'] = record['bulk1_ref_AD'] + record['bulk1_alt_AD']
        record['bulk2_depth'] = record['bulk2_ref_AD'] + record['bulk2_alt_AD']
        if record['parent_depth'] < self.minDP or record['parent_depth'] > self.maxDP:
            record['type'] = 'discard'
        elif record['bulk1_depth'] < self.minDP or record['bulk1_depth'] > self.maxDP:
            record['type'] = 'discard'
        elif record['bulk2_depth'] < self.minDP or record['bulk2_depth'] > self.maxDP:
            record['type'] = 'discard'
        return record

    def filt_index(self, record):
        record['bulk1_SNPindex'] = record['bulk1_alt_AD']/record['bulk1_depth']
        record['bulk2_SNPindex'] = record['bulk2_alt_AD']/record['bulk2_depth']
        if (record['bulk1_SNPindex'] < self.min_SNPindex) and \
           (record['bulk2_SNPindex'] < self.min_SNPindex):
            record['type'] = 'discard'
        else:
            record['delta_SNPindex'] = record['bulk2_SNPindex'] - record['bulk1_SNPindex']
        return record

    def filt_strand_bias(self, record, ADFR):
        ADF = ADFR[0]
        ADR = ADFR[1]

        if record['parent_GT'] == '0/0':
            cultivar_ADF = int(ADF.split(',')[0])
            cultivar_ADR = int(ADR.split(',')[0])
        else:
            cultivar_ADF = int(ADF.split(',')[1])
            cultivar_ADR = int(ADR.split(',')[1])

        if cultivar_ADF == 0 and cultivar_ADR > self.strand_bias:
            record['type'] = 'discard'
        elif cultivar_ADR == 0 and cultivar_ADF > self.strand_bias:
            record['type'] = 'discard'
        return record


    def filt(self, parent_GT, parent_AD, bulk1_AD, bulk2_AD, ADFR):
        record = self.filt_parent_GT(parent_GT, bulk1_AD, bulk2_AD)
        if record['type'] == 'keep':
            record = self.filt_depth(record, parent_AD)
            if record['type'] == 'keep':
                record = self.filt_index(record)
                if record['type'] == 'keep':
                    if ADFR != None:
                        record = self.filt_strand_bias(record, ADFR)

        return record
