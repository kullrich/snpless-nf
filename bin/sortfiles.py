#!/usr/bin/env python3


import sys
import glob
import operator


filepath = sys.argv[1]
fileending = str(sys.argv[2])
filesplit = str(sys.argv[3])


files = glob.glob(filepath+'/*.'+fileending)
splitted = [x.replace(filepath+'/','').replace('.'+fileending,'').strip().split(filesplit) for x in files]
splitted = [[x[0],int(x[1]),int(x[2]),x[3]] for x in splitted]
splitted.sort(key = operator.itemgetter(1, 2))
filesout = ['_'.join([str(x) for x in x]) for x in splitted]
filesout = [filepath+'/'+x+'.'+fileending for x in filesout]
[sys.stdout.write(x+'\n') for x in filesout]
