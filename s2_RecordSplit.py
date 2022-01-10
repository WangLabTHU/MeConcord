# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 20:06:00 2021

@author: xianglin
"""

#%% read pars from environment

import sys,getopt
opts,args = getopt.getopt(sys.argv[1:],"hi:o:g:")
chroms = ['chr'+str(x) for x in range(1,23)]
for op,val in opts:
    if op == '-i':
        infile = val
    elif op == '-o':
        outfile = val
    elif op == '-g':
        chroms = val.split(',')
    elif op == '-h':
        print('Usage: python s2_RecordSplit.py -i ./test_ReadsMethyAndMuts.txt -o ./test -g chr1,chr2,chr3,chr4,chr5\n\
              i,The path to s1 output. ( end with _ReadsMethyAndMuts.txt);\n\
              o,Output prefix;\n\
              g,Chromosomes used; (default chromsome 1-22); chromosomes shoud be seperated by comma;\n\
              h,Help information')
        sys.exit()

if __name__ == '__main__':
    for chrom in chroms:
        cc = 0
        out_tmp = open(outfile +'_ReadsMethyAndMuts_'+chrom+'.txt','w')
        with open(infile,'r') as f:
            for line in f.readlines():
                if cc == 0:
                    cc = 1
                else:
                    line1 = line.strip().split('\t')
                    if line1[0] == chrom and len(line1[6] + line1[7]) >= 2:
                        out_tmp.write(line)
            f.close()
        out_tmp.close()
    