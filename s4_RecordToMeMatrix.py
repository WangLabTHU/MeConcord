# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 19:00:19 2021
@author: xianglin
"""

import numpy as np
import pandas as pd
import math
import sys,getopt
basedon0 = 0
max_dis = 600
chrom_used = ['chr'+str(x) for x in range(1,23)]
opts,args = getopt.getopt(sys.argv[1:],"hi:o:r:m:c:z:g:")
for op,val in opts:
    if op == '-i':
        infile = val
    elif op == '-o':
        outpre = val
    elif op == '-r':
        binfile = val
    elif op == '-c':
        cpgpos_path = val
    elif op == '-m':
        max_dis = int(val)
    elif op == '-z':
        basedon0 = int(val)
    elif op == '-g':
        chrom_used = val.split(',')
    elif op == '-h':
        print('Usage: python s4_RecordToMeMatrix.py -i ./test -o ./test -r ./region.bed -c ./cpgpos/ -m 600 -z 0 -g chr1,chr2\n\
              i,The path to s2_RecordSplit.py output, with prefix file name;\n\
              o,Output prefix;\n\
              r,The files with genomic regions for computation, chrom, start, end seperated by tab;\n\
              c,Cpg position folder, output of pre_cpg_pos.py;\n\
              z,Whether is the genomic file based on 0; 0 (default) or 1; output is same to input bins; if -r is a bed file, -z should be 1;\n\
              g,Chromosomes used; (default chromsome 1-22); chromosomes shoud be seperated by comma;\n\
              m,Maximum of reads length (default 600bp for paired-end reads). if there are single-end reads,\
m should be set length of reads, if not sure, default will work for most cases;')
        sys.exit()
def extracting_methy(reads1,reads2,pos1,pos2,cpg_pos_wanted,methy_pattern):
    bin_reads_me_end = np.zeros([1,cpg_pos_wanted.shape[0]],int)
    bin_reads_unme_end = np.zeros([1,cpg_pos_wanted.shape[0]],int)
    methy_me = np.array([int(x) for x in methy_pattern.replace('Z','1').replace('z','0').replace('N','0')]).reshape(1,-1)
    methy_unme = np.array([int(x) for x in methy_pattern.replace('Z','0').replace('z','1').replace('N','0')]).reshape(1,-1)
    if reads1>=pos1 and reads1<=pos2:
        relative_ind1 = sum(reads1>cpg_pos_wanted.iloc[:,1])
        relative_ind2 = min(len(methy_pattern)+relative_ind1,cpg_pos_wanted.shape[0])
        bin_reads_me_end[0,relative_ind1:relative_ind2] = methy_me[0,0:(relative_ind2-relative_ind1)]
        bin_reads_unme_end[0,relative_ind1:relative_ind2] = methy_unme[0,0:(relative_ind2-relative_ind1)]
    elif reads2>=pos1 and reads2<=pos2:
        relative_ind1 = 0
        relative_ind2 = sum(reads2-1>=cpg_pos_wanted.iloc[:,1])
        bin_reads_me_end[0,relative_ind1:relative_ind2] = methy_me[0,-relative_ind2:]
        bin_reads_unme_end[0,relative_ind1:relative_ind2] = methy_unme[0,-relative_ind2:]
    return bin_reads_me_end,bin_reads_unme_end

if __name__ == '__main__':
    min_cpg = 2
    min_read = 2
    filter_reads_cpgnum = 2
    bindata = pd.read_csv(binfile, sep = '\t',header = None,)
    for chrom in chrom_used:
        bin_chrom = bindata.loc[bindata.iloc[:,0]==chrom,:].copy()
        if bin_chrom.shape[0]>=1:
            bin_sort = bin_chrom.sort_values(by=[2],ascending = [True])
            print('We are processing '+chrom +' now!')
            reads_chrom0 = pd.read_csv(infile+'_ReadsMethyAndMuts_'+chrom+'.txt',sep = '\t',header = None,)
            reads_chrom0.iloc[:,6].fillna('',inplace = True)
            reads_chrom0.iloc[:,7].fillna('',inplace = True)
            reads_chrom = reads_chrom0.loc[reads_chrom0.iloc[:,6].map(len) + reads_chrom0.iloc[:,7].map(len)>=filter_reads_cpgnum,:].copy()
            reads_sort = reads_chrom.sort_values(by=[2,3],ascending = [True,True])
            if basedon0 == 1:
                bin_sort.iloc[:,1] = bin_chrom.iloc[:,1]+1 ##has been converted to 1-based
            cpg_pos = pd.read_csv(cpgpos_path+'/cpgpos_'+chrom+'.pos',sep = '\t',header = None,)
            idx = 0
            for i in range(0, bin_sort.shape[0]):
                print('We are processing genomic intervals of '+chrom+':'+str(bin_sort.iloc[i,1])+\
                      '-'+str(bin_sort.iloc[i,2])+'...')
                pos1 = bin_sort.iloc[i,1]
                pos2 = bin_sort.iloc[i,2]
                cpg_pos_wanted = cpg_pos.loc[(cpg_pos.iloc[:,1]>=pos1) & (cpg_pos.iloc[:,1]<=pos2),:]
                cpg_pos_help = cpg_pos.loc[(cpg_pos.iloc[:,1]>=(pos1-1000000)) & (cpg_pos.iloc[:,1]<=pos2),:]
                if cpg_pos_wanted.shape[0]>=min_cpg:##bin has at least cpgs
                    bin_reads_me = np.zeros([0,cpg_pos_wanted.shape[0]],int)
                    bin_reads_unme = np.zeros([0,cpg_pos_wanted.shape[0]],int)
                    while idx >0 and pos1 - reads_sort.iloc[idx,3]< max_dis:
                        idx -= 1
                    while idx < (reads_sort.shape[0]-1) and pos1 - reads_sort.iloc[idx,3] > max_dis:
                        idx += 1
                    while idx < (reads_sort.shape[0]-1) and reads_sort.iloc[idx,2] < pos2:
                        bin_reads_me_read = np.zeros([1,cpg_pos_wanted.shape[0]],int)
                        bin_reads_unme_read = np.zeros([1,cpg_pos_wanted.shape[0]],int)
                        if len(reads_sort.iloc[idx,6]) >=1:
                            if sum((cpg_pos_wanted.iloc[:,1]>=reads_sort.iloc[idx,2]) & (cpg_pos_wanted.iloc[:,1]<reads_sort.iloc[idx,3])) >= 1:
                                if reads_sort.iloc[idx,2] < pos1 and reads_sort.iloc[idx,2] >= cpg_pos_help.iloc[0,1]:
                                    lost_num = sum((reads_sort.iloc[idx,2]<=cpg_pos_help.iloc[:,1]) & (cpg_pos_wanted.iloc[0,1]>cpg_pos_help.iloc[:,1]))
                                    methy_reads = reads_sort.iloc[idx,6][lost_num:]
                                    read_start = cpg_pos_wanted.iloc[0,1]
                                    read_end = reads_sort.iloc[idx,3]
                                else:
                                    methy_reads = reads_sort.iloc[idx,6]
                                    read_start = reads_sort.iloc[idx,2]
                                    read_end = reads_sort.iloc[idx,3]
                                bin_reads_me_end,bin_reads_unme_end = extracting_methy(read_start,read_end,\
                                                                                       pos1,pos2,cpg_pos_wanted,methy_reads)
                                bin_reads_me_read += bin_reads_me_end
                                bin_reads_unme_read += bin_reads_unme_end
                        if len(reads_sort.iloc[idx,7]) >= 1:
                            if sum((cpg_pos_wanted.iloc[:,1]>=reads_sort.iloc[idx,4]) & (cpg_pos_wanted.iloc[:,1]<reads_sort.iloc[idx,5])) >= 1:
                                if reads_sort.iloc[idx,4] < pos1 and reads_sort.iloc[idx,4] >= cpg_pos_help.iloc[0,1]:
                                    lost_num = sum((reads_sort.iloc[idx,4]<=cpg_pos_help.iloc[:,1]) & (cpg_pos_wanted.iloc[0,1]>cpg_pos_help.iloc[:,1]))
                                    methy_reads = reads_sort.iloc[idx,7][lost_num:]
                                    read_start = cpg_pos_wanted.iloc[0,1]
                                    read_end = reads_sort.iloc[idx,5]
                                else:
                                    methy_reads = reads_sort.iloc[idx,7]
                                    read_start = reads_sort.iloc[idx,4]
                                    read_end = reads_sort.iloc[idx,5]
                                bin_reads_me_end,bin_reads_unme_end = extracting_methy(read_start,read_end,\
                                                                                       pos1,pos2,cpg_pos_wanted,methy_reads)
                                bin_reads_me_read += bin_reads_me_end
                                bin_reads_unme_read += bin_reads_unme_end
                        if bin_reads_me_read.sum() + bin_reads_unme_read.sum()>=filter_reads_cpgnum:##at least 3 cpgs for each fragment
                            bin_reads_me = np.vstack([bin_reads_me, bin_reads_me_read])
                            bin_reads_unme = np.vstack([bin_reads_unme, bin_reads_unme_read])
                        idx += 1
                    #calculating bin 3 metrics
                    if basedon0 == 1:
                        pos1 = pos1-1 #to 0-based
                    name1 = outpre + '_'+chrom+'_'+str(pos1)+'_'+str(pos2)+'_me.txt'
                    name2 = outpre + '_'+chrom+'_'+str(pos1)+'_'+str(pos2)+'_unme.txt'
                    np.savetxt(name1,bin_reads_me,fmt='%d',delimiter='\t')
                    np.savetxt(name2,bin_reads_unme,fmt='%d',delimiter='\t')

