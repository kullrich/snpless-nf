#!/usr/bin/env python3


import sys
import operator


filepath = sys.argv[1]


files = filepath.split()
for f in files:
    foo = open(f,'rt')
    next(foo)
    fcov = next(foo).strip().split('\t')[6]
    foo.close()
    sys.stdout.write(f.replace('.mean.coverage','').replace('.breseq', '').replace('.minimap2', '').replace('.bwa', '')+'\t'+fcov+'\n')
