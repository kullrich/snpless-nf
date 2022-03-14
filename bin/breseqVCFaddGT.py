#!/usr/bin/env python3


import sys
import numpy as np


vcf = sys.argv[1]
sample_id = sys.argv[2]
with open(vcf, 'rt') as inhandle:
    for lines in inhandle:
        if lines != '##fileDate\n':
            if lines[:6] != '#CHROM' and lines[0] == '#':
                sys.stdout.write(lines)
            if lines[:6] == '#CHROM' and lines[0] == '#':
                sys.stdout.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
                sys.stdout.write('##FORMAT=<ID=AD,Number=1,Type=Float,Description="Allele Depth (avg read count)">\n')
                sys.stdout.write('##FORMAT=<ID=DP,Number=1,Type=Float,Description="Total Depth (avg read count)">\n')
                sys.stdout.write(lines.strip() + '\tFORMAT\t' + sample_id + '\n')
            if lines[:6] != '#CHROM' and lines[0] != '#':
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = lines.strip().split('\t')
                AF, AD, DP = INFO.split(';')
                AF = np.float64(AF.replace('AF=',''))
                AD = np.float64(AD.replace('AD=',''))
                DP = np.float64(DP.replace('DP=',''))
                if AF == 1:
                    sys.stdout.write(lines.strip() + '\tGT:AF:AD:DP\t1/1:' + str(AF) + ':' + str(AD) + ':' + str(DP) + '\n')
                if AF != 1:
                    sys.stdout.write(lines.strip() + '\tGT:AF:AD:DP\t0/1:' + str(AF) + ':' + str(AD) + ':' + str(DP) + '\n')
