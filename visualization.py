# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 10:20:28 2022

@author: xianglin.Zhang
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

infile = './tmp_test/test_chr1_1287967_1288117'
cpgfile = './tmp_test/'
outfile = './'
zero_based = 0
import sys,getopt
opts,args = getopt.getopt(sys.argv[1:],"hi:o:c:z:")
for op,val in opts:
    if op == '-i':
        infile = val
    elif op == '-o':
        outfile = val
    elif op == '-c':
        cpgfile = val
    elif op == '-z':
        zero_based = int(val)
    elif op == '-h':
        print('Usage: python visualization.py -i ./test/test_chr1_1287967_1288117 -o ./test/test -c ./cpgpos/ \n\
              i,The prefix path to recordMat (prefix before _me.txt or _unme.txt);\n\
              o,The prefix path that you want to deposit pdf files;\n\
              c,Cpg position folder, output of pre_cpg_pos.py;\n\
              h,Help information')
        sys.exit()

unmemat = pd.read_csv(infile+'_unme.txt',sep = '\t',header = None,)
memat = pd.read_csv(infile+'_me.txt',sep = '\t',header = None,)
merge = memat - unmemat
fig = plt.figure(figsize = [merge.shape[1]/2,merge.shape[0]/3])
ax = fig.add_subplot(111)
for i in range(0, merge.shape[0]):
    x = np.array([0,merge.shape[1]+1])
    y = np.array([i+1,i+1])
    
    x2 = np.arange(1,merge.shape[1]+1)
    y2 = np.repeat(i+1,merge.shape[1])
    c2 = np.array(merge.iloc[merge.shape[0]-1-i,:])
    cz2 = np.zeros([merge.shape[1],3])
    cz2[c2==-1,:] = np.repeat([0.68,0.92,1],sum(c2==-1),axis = 0).reshape((sum(c2==-1),3),order = 'F')
    cz2[c2==0,:] = np.repeat([0.8,0.8,0.8],sum(c2==0),axis = 0).reshape((sum(c2==0),3),order = 'F')
    cz2[c2==1,:] = np.repeat([0.6,0.2,0],sum(c2==1),axis = 0).reshape((sum(c2==1),3),order = 'F')
    ax.scatter(x2,y2,s = 150,c = cz2,edgecolors = 'k',alpha = 1)
    ax.plot(x,y,'k',linewidth=1.0,zorder= 0)
ax.set(xlim=(0, merge.shape[1]+1), xticks=np.arange(1, merge.shape[1]+1),
       ylim=(0, merge.shape[0]+1), yticks=np.arange(1, merge.shape[0]+1))

names = infile.strip().split('_')
if zero_based == 0:
    anno = '(genomic position is 1-based)'
    pos1 = int(names[-2])
else:
    anno = '(genomic position is 0-based)'
    pos1 = int(names[-2])+1
ax.set_title(names[-3]+':'+names[-2]+'-'+names[-1]+'bp'+ '\n'+ anno)
ax.set_xlabel("# CpG sites")
ax.set_ylabel("# Reads")
#plt.show()
fig.savefig(outfile + '_'+names[-3]+'_'+names[-2]+'_'+names[-1]+'_cpgnum.pdf',dpi=300)


## with distance of genomic distances
cpg_pos = pd.read_csv(cpgfile+'/cpgpos_'+names[-3]+'.pos',sep = '\t', header = None)
posw = cpg_pos.loc[np.logical_and(cpg_pos.iloc[:,1]>=pos1,cpg_pos.iloc[:,1]<int(names[-1])),1]
relative_posw = posw-pos1
lens = int(names[-1]) - pos1+1

fig2 = plt.figure(figsize = [8,merge.shape[0]/3])
ax = fig2.add_subplot(111)
for i in range(0, merge.shape[0]):
    x = np.array([-5,lens+5])
    y = np.array([i+1,i+1])
    
    x2 = np.array(relative_posw)
    y2 = np.repeat(i+1,merge.shape[1])
    c2 = np.array(merge.iloc[merge.shape[0]-1-i,:])
    cz2 = np.zeros([merge.shape[1],3])
    cz2[c2==-1,:] = np.repeat([0.68,0.92,1],sum(c2==-1),axis = 0).reshape((sum(c2==-1),3),order = 'F')
    cz2[c2==0,:] = np.repeat([0.8,0.8,0.8],sum(c2==0),axis = 0).reshape((sum(c2==0),3),order = 'F')
    cz2[c2==1,:] = np.repeat([0.6,0.2,0],sum(c2==1),axis = 0).reshape((sum(c2==1),3),order = 'F')
    ax.scatter(x2,y2,s = 150,c = cz2,edgecolors = 'k',alpha = 1)
    ax.plot(x,y,'k',linewidth=1.0,zorder= 0)
ax.set(xlim=(-5, lens+5), xticks=np.arange(0, lens+5,50),
       ylim=(0, merge.shape[0]+1), yticks=np.arange(1, merge.shape[0]+1))
ax.set_title(names[-3]+':'+names[-2]+'-'+names[-1]+'bp'+ '\n'+ anno)
ax.set_xlabel("Relative genomic positions (bp)")
ax.set_ylabel("# Reads")
#plt.show()
fig2.savefig(outfile + '_'+names[-3]+'_'+names[-2]+'_'+names[-1]+'_cpgpos.pdf',dpi=300)
