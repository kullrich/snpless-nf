#!/usr/bin/env python3


import sys
import operator
import numpy as np
from cyvcf2 import VCF

#mean_coverage input file should have <sample_id><tab><meanDepth>

vcf = VCF(sys.argv[1])
mean_coverage = sys.argv[2]
mean_coverage_dict = {}
with open(mean_coverage, 'rt') as inhandle:
    for lines in inhandle:
        sample, depth = lines.strip().split('\t')
        mean_coverage_dict[sample.replace('.breseq','').replace('.minimap2','').replace('.bwa','')] = float(depth)
samples_split = [x.split("_") for x in vcf.samples]
samples_split = [x for x in samples_split if x[3] != "c"]
splitted = [[x[0],int(x[1]),int(x[2]),x[3]] for x in samples_split]
splitted.sort(key = operator.itemgetter(1, 2))
samples_ordered = ['_'.join([str(x) for x in x]) for x in splitted]
samples_index = [samples_ordered.index(x) for x in vcf.samples if x in samples_ordered]
samples_timepoints = ' '.join([str(x) for x in range(len(samples_ordered))])
samples_zero = ' '.join([str(0) for x in range(len(samples_ordered))])
samples_coverage = ' '.join([str(mean_coverage_dict[x]) for x in samples_ordered])
#chromosome, pos, depth, timepoints <0 1 2 3 ...>, alternative depth <int int int ...>, mean depth <int int int ...> 
var_out_header = ', '.join(['chromosome', '0', 'Depth', samples_timepoints, samples_zero, samples_coverage])
var_out = []
for variant in vcf:
    #chromosome, pos, depth, timepoints <0 1 2 3 ...>, alternative depth <int int int ...>, total depth <int int int ...> 
    if variant.var_type == 'snp':
        var_out.append(', '.join([variant.CHROM, str(variant.POS), variant.REF+'->'+variant.ALT[0],samples_timepoints,' '.join([str(x) for x in variant.gt_alt_depths]),' '.join([str(x) for x in np.array(variant.gt_depths) - np.array(variant.gt_alt_depths)])]))
    if variant.var_type == 'deletion':
        var_out.append(', '.join([variant.CHROM, str(variant.POS), 'deletion',samples_timepoints,' '.join([str(x) for x in variant.gt_alt_depths]),' '.join([str(x) for x in np.array(variant.gt_depths) - np.array(variant.gt_alt_depths)])]))
    if variant.var_type == 'indel':
        var_out.append(', '.join([variant.CHROM, str(variant.POS), 'indel',samples_timepoints,' '.join([str(x) for x in variant.gt_alt_depths]),' '.join([str(x) for x in np.array(variant.gt_depths) - np.array(variant.gt_alt_depths)])]))
sys.stdout.write(var_out_header+'\n')
[sys.stdout.write(x+'\n') for x in var_out]