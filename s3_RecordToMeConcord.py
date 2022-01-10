# -*- coding: utf-8 -*-
"""
@author: xianglin
"""

import numpy as np
import pandas as pd
import math
import sys,getopt
import random
from scipy import stats
import multiprocessing
import os
opts,args = getopt.getopt(sys.argv[1:],"hi:o:r:m:c:b:p:z:g:")
interval_size = 150
cores = 4
basedon0 = 0
chrom_used = ['chr'+str(x) for x in range(1,23)]
for op,val in opts:
    if op == '-i':
        infile = val
    elif op == '-o':
        outpre = val
    elif op == '-r':
        binfile = val
    elif op == '-c':
        cpgpos_path = val
    elif op == '-b':
        interval_size = int(val)
    elif op == '-m':
        max_dis = int(val)
    elif op == '-p':
        cores = int(val)
    elif op == '-z':
        basedon0 = int(val)
    elif op == '-g':
        chrom_used = val.split(',')
    elif op == '-h':
        print('Usage: python s3_RecordToMeConcord.py -p 4 -i ./test -o ./test -r ./region.bed -c ./cpgpos/ -b 150 -m 600 -z 0 -g chr1,chr2,chr3\n\
              i,The path to s2_RecordSplit.py output, with prefix file name;\n\
              p,Threads used for parallel computation; default is 4;\n\
              o,Output prefix;\n\
              r,The files with genomic regions for computation, chrom, start, end seperated by tab;\n\
              c,Cpg position folder, output of pre_cpg_pos.py;\n\
              b,Bin size (default 150bp);\n\
              z,Whether is the genomic file based on 0; 0 (default) or 1; output is same to input bins; if -r is a bed file, -z should be 1;\n\
              g,Chromosomes used; (default chromsome 1-22); chromosomes shoud be seperated by comma;\n\
              m,Maximum of fragement length in sequencing library(default 600bp for paired-end reads). if there are single-end reads,\
m should be set as the length of reads, if not sure, default will work for most cases;')
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

def each_chrom(outpre,infile,cpgpos_path,binfile,chrom):
    # inner pars
    min_cpg = 2
    min_read = 2
    min_overlap = 2
    max_read = 500
    filter_reads_cpgnum = 2
    #inner loading and saving
    outdata = open(outpre + '_' + chrom + '_MeConcord.txt','w')
    reads_chrom0 = pd.read_csv(infile+'_ReadsMethyAndMuts_'+chrom+'.txt',sep = '\t',header = None,)
    reads_chrom0.iloc[:,6].fillna('',inplace = True)
    reads_chrom0.iloc[:,7].fillna('',inplace = True)
    reads_chrom = reads_chrom0.loc[reads_chrom0.iloc[:,6].map(len) + reads_chrom0.iloc[:,7].map(len)>=filter_reads_cpgnum,:].copy()
    reads_sort = reads_chrom.sort_values(by=[2,3],ascending = [True,True])
    bindata = pd.read_csv(binfile, sep = '\t',header = None,)
    bin_chrom = bindata.loc[bindata.iloc[:,0]==chrom,:].copy()
    bin_sort = bin_chrom.sort_values(by=[2],ascending = [True])
    if basedon0 == 1:
        bin_sort.iloc[:,1] = bin_chrom.iloc[:,1]+1 ##has been converted to 1-based
    cpg_pos = pd.read_csv(cpgpos_path+'/cpgpos_'+chrom+'.pos',sep = '\t',header = None,)
    #processing
    idx = 0
    for i in range(0, bin_sort.shape[0]):
        split_num = int(math.ceil((bin_sort.iloc[i,2]-bin_sort.iloc[i,1]+1)/interval_size))
        if basedon0 == 1:
            interval_name = chrom+'_'+str(bin_sort.iloc[i,1]-1)+'_'+str(bin_sort.iloc[i,2])#to 0-based
        else:
            interval_name = chrom+'_'+str(bin_sort.iloc[i,1])+'_'+str(bin_sort.iloc[i,2])
        for j in range(0,split_num):
            pos1 = bin_sort.iloc[i,1]+j*interval_size
            pos2 = min(bin_sort.iloc[i,2], bin_sort.iloc[i,1]+(j+1)*interval_size-1)
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
                    if bin_reads_me_read.sum() + bin_reads_unme_read.sum()>=2:
                        bin_reads_me = np.vstack([bin_reads_me, bin_reads_me_read])
                        bin_reads_unme = np.vstack([bin_reads_unme, bin_reads_unme_read])
                    idx += 1
                #calculating bin 3 metrics
                if bin_reads_me.shape[0] >= 1:
                    bin_reads_total = bin_reads_me+bin_reads_unme
                    methy_level = round(bin_reads_me.sum()*1.0/bin_reads_total.sum(),3)
                    total_me = bin_reads_me.sum()
                    total_cpg = bin_reads_total.sum()
                else:
                    methy_level = np.nan
                    total_me = 0
                    total_cpg = 0
                if bin_reads_me.shape[0]>=min_read and bin_reads_me.shape[1]>=min_cpg and bin_reads_total.sum(axis = 0).max() >= min_overlap:
                    if bin_reads_me.shape[0]>max_read:
                        f1 = random.sample(range(bin_reads_me.shape[0]),max_read)
                        bin_reads_me = bin_reads_me[f1,:]
                        bin_reads_unme = bin_reads_unme[f1,:]
                        bin_reads_total = bin_reads_total[f1,:]
                        downsample = 1
                    else:
                        downsample = 0
                    #concordant of reads
                    reads_me = bin_reads_me.dot(bin_reads_me.T)
                    reads_unme = bin_reads_unme.dot(bin_reads_unme.T)
                    reads_total = bin_reads_total.dot(bin_reads_total.T)
                    help_mat = np.ones([bin_reads_me.shape[0],bin_reads_me.shape[0]],int) - np.eye(bin_reads_me.shape[0],dtype = int)
                    concordant_reads = round((reads_me*help_mat+reads_unme*help_mat).sum()*1.0/((reads_total*help_mat).sum()),3)
                    #pvals for concordant of reads
                    total_pair = (reads_total*help_mat).sum()
                    wanted_pair = (reads_me*help_mat+reads_unme*help_mat).sum()
                    pair_me_count = ((bin_reads_me.dot(bin_reads_total.T))*help_mat).sum()
                    pair_unme_count = ((bin_reads_unme.dot(bin_reads_total.T))*help_mat).sum()
                    pair_me_frac = pair_me_count*1.0/(pair_me_count+pair_unme_count)
                    bino_rat = pair_me_frac*pair_me_frac+(1-pair_me_frac)*(1-pair_me_frac)
                    exp_reads = round(bino_rat,3)
                    if pair_me_count != 0 and pair_unme_count != 0:
                        if wanted_pair >= total_pair*bino_rat:
                            p_reads = 1-stats.binom.cdf(wanted_pair-1,total_pair,bino_rat)
                        else:
                            p_reads = stats.binom.cdf(wanted_pair,total_pair,bino_rat)
                    else:
                        p_reads = 1
                    #concordant of cpg site
                    site_me = bin_reads_me.T.dot(bin_reads_me)
                    site_unme = bin_reads_unme.T.dot(bin_reads_unme)
                    site_total = bin_reads_total.T.dot(bin_reads_total)
                    help_mat = np.ones([bin_reads_me.shape[1],bin_reads_me.shape[1]],int) - np.eye(bin_reads_me.shape[1],dtype = int)
                    concordant_sites = round((site_me*help_mat+site_unme*help_mat).sum()*1.0/((site_total*help_mat).sum()),3)
                    #pvals for concordant of CpGs
                    total_pair = (site_total*help_mat).sum()
                    wanted_pair = (site_me*help_mat+site_unme*help_mat).sum()
                    pair_me_count = ((bin_reads_me.T.dot(bin_reads_total))*help_mat).sum()
                    pair_unme_count = ((bin_reads_unme.T.dot(bin_reads_total))*help_mat).sum()
                    pair_me_frac = pair_me_count*1.0/(pair_me_count+pair_unme_count)
                    bino_rat = pair_me_frac*pair_me_frac+(1-pair_me_frac)*(1-pair_me_frac)
                    exp_cpgs = round(bino_rat,3)
                    if pair_me_count != 0 and pair_unme_count != 0:
                        if wanted_pair >= total_pair*bino_rat:
                            p_cpgs = 1-stats.binom.cdf(wanted_pair-1,total_pair,bino_rat)
                        else:
                            p_cpgs = stats.binom.cdf(wanted_pair,total_pair,bino_rat)
                    else:
                        p_cpgs = 1
                    
                else:
                    downsample = 0
                    concordant_reads = np.nan
                    concordant_sites = np.nan
                    p_reads = 1
                    exp_reads = np.nan
                    p_cpgs = 1
                    exp_cpgs = np.nan
                if basedon0 == 1:
                    pos1 = pos1-1 #to 0-based
                outdata.write(interval_name+ '\t'+\
                              chrom+'\t'+\
                              str(pos1)+'\t'+\
                              str(pos2)+'\t'+\
                              str(bin_reads_me.shape[0])+'\t'+\
                              str(bin_reads_me.shape[1])+'\t'+\
                              str(total_me)+'\t'+\
                              str(total_cpg)+'\t'+\
                              str(methy_level)+'\t'+\
                              str(concordant_reads)+'\t'+\
                              str(concordant_sites)+'\t'+\
                              str(concordant_reads-exp_reads)+'\t'+\
                              '%.3e' % p_reads+'\t'+\
                              str(concordant_sites-exp_cpgs)+'\t'+\
                              '%.3e' % p_cpgs+'\t'+\
                              str(downsample)+'\n')
            else:
                pass
    outdata.close()
    
def merge_each_chrom(args):
	each_chrom(*args)

if __name__ == '__main__':
    pars = [(outpre,infile,cpgpos_path,binfile,chrom) for chrom in chrom_used]
    pool = multiprocessing.Pool(processes=cores)
    pool.map(merge_each_chrom,pars)
    pool.close()
    #integrating
    outfile = outpre + '_MeConcord.txt'
    outdata = open(outfile,'w')
    outdata.write('interval'+'\t'+\
                  'chrom'+'\t'+\
                  'start'+'\t'+\
                  'end'+'\t'+\
                  'ReadNum'+'\t'+\
                  'CpGNum'+'\t'+\
                  'MeCpG'+'\t'+\
                  'TotalCpG'+'\t'+\
                  'DNAme'+'\t'+\
                  'RC'+'\t'+\
                  'CC'+'\t'+\
                  'NRC'+'\t'+\
                  'pr'+'\t'+\
                  'NCC'+'\t'+\
                  'pc'+'\t'+\
                  'IfSubsampled'+'\n')
    for chrom in chrom_used:
        tmpfile = open(outpre+ '_' + chrom + '_MeConcord.txt','r')
        for line in tmpfile.readlines():
            outdata.write(line)
        tmpfile.close()
        if (os.path.isfile(outpre+ '_' + chrom + '_MeConcord.txt')):
            os.remove(outpre+ '_' + chrom + '_MeConcord.txt')
    outdata.close()

