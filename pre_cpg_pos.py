# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:55:52 2021

@author: xianglin.zhang
"""

import pandas as pd
import math

import sys,getopt
opts,args = getopt.getopt(sys.argv[1:],"hi:o:")
for op,val in opts:
    if op == '-i':
        infile = val
    elif op == '-o':
        outfile = val
    elif op == '-h':
        print('Usage: python pre_cpg_pos.py -i hg38.fa -o ./cpg_pos/\n\
              i,The path to reference sequences (.fa);\n\
              o,The path that you want to deposit the positions of CpG sites, each chromosome has a seperate file;\n\
              h,Help information')
        sys.exit()

##
def deriveseq(data):
    chrom_seq = ''
    chrom_name = []
    while True:
        f = data.readline()
        if not f:
            break
        if f[0] == '>':
            #call the pos of cpgs
            if chrom_seq == '':
                cpg_pos = []
            else:
                chrom_seq = chrom_seq.upper()
                cpgpos = detercpgpos(chrom_seq, final_name)
                cpg_pos = cpg_pos + cpgpos
                chrom_seq = ''
            #call next chromosome info
            #chromosome name
            raw_name = f[1:].strip()
            pos1 = raw_name.find(' ')
            pos2 = raw_name.find('\t')
            if pos1 == -1 and pos2 != -1:
                final_name = raw_name[0:pos2]
            elif pos1 !=-1 and pos2 == -1:
                final_name = raw_name[0:pos1]
            elif pos1 == -1 and pos2 == -1:
                final_name = raw_name
            elif pos1 != -1 and pos2 != -1:
                final_name = raw_name[0:min(pos1,pos2)]
            chrom_name.append(final_name)
            print('We are going to process '+final_name +' now!')
        else:
            chrom_seq += f.strip()
    if chrom_seq != '':
        chrom_seq = chrom_seq.upper()
        cpgpos = detercpgpos(chrom_seq, final_name)
        cpg_pos = cpg_pos + cpgpos
        chrom_seq = ''
    return cpg_pos,chrom_name

def detercpgpos(chrom_seq,final_name):
    chrom_len = len(chrom_seq)
    cpgpos = []
    for i in range(0,chrom_len-1):
        if math.ceil(i/20000000)*20000000 == i:
            print('Processing '+str(i/1000000)+ ' Mb...')
        if chrom_seq[i:(i+2)] == 'CG':
            tmp_pos = [[final_name,i+1]]
            cpgpos += tmp_pos
        else:
            pass
    return cpgpos

if __name__ == "__main__":
    data = open(infile,'r')
    cpg_pos,chrom_name = deriveseq(data)
    out_pos = pd.DataFrame(cpg_pos,columns=['chrom','pos_C'])
    for item in chrom_name:
        chr_pos = out_pos.loc[out_pos.iloc[:,0] == item,:].copy()
        chr_pos.to_csv(outfile+'/'+'cpgpos_'+item+'.pos',sep = '\t',header = False,index = False)