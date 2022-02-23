#!/usr/bin/env python3


import sys
import operator


covlist = sys.argv[1]


for f in open(covlist, 'rt'):
    foo = open(f.strip(),'rt')
    next(foo)
    fcov = next(foo).strip().split('\t')[6]
    foo.close()
    sys.stdout.write(f.strip().split('/')[-1].replace('.mean.coverage','').replace('.breseq', '').replace('.minimap2', '').replace('.bwa', '')+'\t'+fcov+'\n')
