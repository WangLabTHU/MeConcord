# -*- coding: utf-8 -*-
"""
python2 working!!!!!
to dissect reads from bismark. with this script, reads from bismark could be converted to methylation of
only CpG on reference watson strand and mutations

reminder1:only for Bismark, cannot be used for BS-seeker2

reminder2:this script could read both sam file and bam file(if read bam file, it requires pysam package) 
and even sam.gz (require gzip package)

reminder3:bam or sam file could contain both paired-end and single-end reads, but it requires two mates 
of paired-end reads are adjacent. Directed output of Bismark or reads-name-sorted files could meet its 
requirments.

reminder4:for mutation calling, insert could be more than 1 base and pos is first base of inserted sequence.
So for one postion on reference, there could be both mutation and insertion meaning first insertion and mutation afterward.

reminder5:mutations types include muations(M), insertion(I), and Deletion(D).

reminder6:ouput based are 1-based. 
    
@author: Xianglin Zhang
"""

#%% read pars from environment

import sys,getopt
opts,args = getopt.getopt(sys.argv[1:],"hi:o:c:")
clipping = 0
for op,val in opts:
    if op == '-i':
        infile = val
    elif op == '-o':
        outfile = val+'_ReadsMethyAndMuts.txt'
    elif op == '-c':
        clipping = int(val)
    elif op == '-h':
        print('Usage: python s1_bamToMeRecord.py -i test.bam -o test -c 0 \n\
              i,The path to input files (.bam or .sam or .sam.gz);\n\
              o,Output prefix;\n\
              c,Clipping read ends with such base number (defalut 0); can be used when sequencing quality of read ends is not good. such as -c 5 to remove 5 bases from the both ends of the reads.\n\
              h,Help information')
        sys.exit()

#%% tmp read test file
'''
infile = './tmp.sam'
outfile = './tmp'+'_ReadsMethyAndMuts.txt'
'''
#%% reading cigar scores
def flags_convert(score):
    par = [128,64,32,16,8,4,2]
    flags_info = []
    for item in par:
        flags_info.append(score/item)
        score = score%item
    flags_info.append(score)
    return flags_info
def cigar_convert(str1):
    left_str = ''
    lens = []
    types = []
    for item in str1:
        if ord(item)>=48 and ord(item)<=57:
            left_str += item
        else:
            lens.append(int(left_str))
            left_str = ''
            types.append(item)
    return lens,types
def seq_reads_clipping(seq_reads,cigar_len,cigar_symbol):
    seq_clipping = seq_reads
    seq_wanted = ''
    insert_pos_anno = ''
    passed_len = 0
    out_insert_seq = []
    out_insert_pass_bases = []
    for i in range(len(cigar_len)):
        element1,element2 = cigar_len[i],cigar_symbol[i]
        if element2 == 'S':
            seq_clipping = seq_clipping[element1:]
        elif element2 == 'M':
            seq_wanted += seq_clipping[:element1]
            insert_pos_anno += '*'*element1
            seq_clipping = seq_clipping[element1:]
            passed_len += element1
        elif element2 == 'D': # changed from bsseeker2 scripts
            seq_wanted += 'D'*element1
            insert_pos_anno += '*'*element1
            passed_len += element1
        elif element2 == 'I': 
            insert_seq = seq_clipping[:element1]
            out_insert_seq.append(insert_seq)
            out_insert_pass_bases.append(passed_len)
            seq_clipping = seq_clipping[element1:]
            if len(insert_pos_anno) > 0:
                insert_pos_anno = insert_pos_anno[:-1] + 'I' #assume ends has no insert/delete
    if len(insert_pos_anno) != len(seq_wanted):
        print('insert_pos_anno length is not equal with seq length!\n'+insert_pos_anno+'\n'+seq_wanted)
        raise KeyboardInterrupt
    return seq_wanted,out_insert_seq,out_insert_pass_bases,insert_pos_anno


def determine_ref(seq_wanted,mut_anno):
    out_del_seq = []
    out_del_pass_bases = []
    if len(mut_anno)>=1:
        ref_wanted = ''
        left_number = ''
        base_id = 0
        dele_sign = 0
        for item in mut_anno:
            if ord(item)>=48 and ord(item)<=57:
                if dele_sign == 1:
                    out_del_seq.append(del_seq)
                    del_seq = ''
                    dele_sign = 0
                left_number += item
            else:
                if item == '^':
                    dele_sign = 1
                    base_id += int(left_number)
                    out_del_pass_bases.append(base_id)
                    ref_wanted += seq_wanted[(base_id-int(left_number)):base_id]
                    left_number = ''
                    del_seq = ''
                elif dele_sign == 1:
                    del_seq += item
                    ref_wanted += item
                    base_id += 1
                else:
                    base_id += int(left_number)
                    ref_wanted += seq_wanted[(base_id-int(left_number)):base_id]
                    ref_wanted += item
                    base_id += 1
                    left_number = ''
        if len(left_number) > 0:
            base_id += int(left_number)
            ref_wanted += seq_wanted[(base_id-int(left_number)):base_id]
            left_number = ''
    else:
        ref_wanted = seq_wanted
    if len(ref_wanted) != len(seq_wanted):
        print('reference length is not equal with seq length!\n'+ref_wanted+'\n'+seq_wanted)
        raise KeyboardInterrupt
    return ref_wanted, out_del_seq, out_del_pass_bases

def find_methy_mut(ref_wanted,seq_wanted,insert_pos_anno,strand):
    if (len(ref_wanted)!=len(seq_wanted)) or (len(ref_wanted)!=len(insert_pos_anno)):
        print('Error in find methy_mut, length not equal!!'+ref_wanted+'\n'+seq_wanted+'\n'+insert_pos_anno)
        raise KeyboardInterrupt
    if strand == '+':
        seq1 = 'TG'
        seq2 = 'C'
        seq3 = 'T'
    else:
        seq1 = 'CA'
        seq2 = 'G'
        seq3 = 'A'
    cpg_methy = ''
    for i in range(0,len(ref_wanted)-1):
        if ref_wanted[i:i+2] == 'CG' and seq_wanted[i:i+2] == 'CG' and insert_pos_anno[i] != 'I':
            cpg_methy += 'Z' # methy
        elif ref_wanted[i:i+2] == 'CG' and seq_wanted[i:i+2] == seq1 and insert_pos_anno[i] != 'I':
            cpg_methy += 'z' # unmethy
        elif ref_wanted[i:i+2] == 'CG':
            cpg_methy += 'N' # no signal or structure of molecule was already not CpG.
    mut_seq = []
    ref_seq = []
    pos_seq = []
    for i in range(0,len(ref_wanted)-1):
        if (ref_wanted[i] != seq2 or seq_wanted[i] != seq3) and seq_wanted[i] != 'D' and ref_wanted[i] != seq_wanted[i]:
            mut_seq.append(seq_wanted[i])
            ref_seq.append(ref_wanted[i])
            pos_seq.append(i)
    return cpg_methy,mut_seq,ref_seq,pos_seq

def get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos):
    Mut_infos = {}
    seq_wanted,out_insert,out_insert_pass_bases,insert_pos_anno =\
    seq_reads_clipping(seq_reads, cigar_len,cigar_symbol)
    if len(out_insert) >=1:
        for num_insert in range(len(out_insert)):
            Mut_infos[str(pos+out_insert_pass_bases[num_insert])+';I;;'+out_insert[num_insert]] = 1
    ref_wanted, out_del_seq, out_del_pass_bases = determine_ref(seq_wanted,mut_anno)
    if len(out_del_seq) >=1:
        for num_del in range(len(out_del_seq)):
            Mut_infos[str(pos+out_del_pass_bases[num_del])+';D;'+out_del_seq[num_del]+';'] = 1
    pos_start = pos
    pos_end = pos + len(ref_wanted)-1 #here is true end. Output has a output end for minus 1 pos in the end.
    if clipping >0 and len(seq_wanted)>2*clipping:
        pos_start = pos_start + clipping
        pos_end = pos_end - clipping
        seq_wanted = seq_wanted[clipping:-clipping]
        ref_wanted = ref_wanted[clipping:-clipping]
        insert_pos_anno = insert_pos_anno[clipping:-clipping]
    return seq_wanted,insert_pos_anno,ref_wanted,pos_start,pos_end,Mut_infos

def update_muts(Mut_infos,mut_seq,ref_seq,pos_seq,pos_start):
    if len(mut_seq)>=1:
        for i in range(len(mut_seq)):
            Mut_infos[str(pos_start+pos_seq[i])+';M;'+ref_seq[i]+';'+mut_seq[i]] = 1
    return Mut_infos

def fine_deoverlap_mates(pos_start1,pos_end1,pos_start2,pos_end2,strand,ref_wanted1,seq_wanted1,insert_pos_anno1,ref_wanted2,seq_wanted2,insert_pos_anno2,Mut_infos1,Mut_infos2):
    Mut_infos = Mut_infos1.copy()
    Mut_infos.update(Mut_infos2)
    if strand == '+':
        x1 = pos_start1
        x2 = pos_end1
        if pos_end1 > pos_start2 and pos_end1+1 <= pos_end2:
            y1 = pos_end1
            y2 = pos_end2
            ref_wanted2 = ref_wanted2[(y1-pos_start2):]
            seq_wanted2 = seq_wanted2[(y1-pos_start2):]
            insert_pos_anno2 = insert_pos_anno2[(y1-pos_start2):]
            m1,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand)
            m2,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand)
            Mut_infos = update_muts(Mut_infos,mut_seq1,ref_seq1,pos_seq1,x1)
            Mut_infos = update_muts(Mut_infos,mut_seq2,ref_seq2,pos_seq2,y1)
        elif pos_end1 <= pos_start2:
            y1 = pos_start2
            y2 = pos_end2
            m1,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand)
            m2,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand)
            Mut_infos = update_muts(Mut_infos,mut_seq1,ref_seq1,pos_seq1,x1)
            Mut_infos = update_muts(Mut_infos,mut_seq2,ref_seq2,pos_seq2,y1)
        elif pos_end1+1 > pos_end2:
            m1,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand)
            Mut_infos = update_muts(Mut_infos,mut_seq1,ref_seq1,pos_seq1,x1)
            m2 = ''
            y1 = ''
            y2 = ''
    else:
        x1 = pos_start2 #R2 to be first!!
        x2 = pos_end2
        if pos_end2 > pos_start1 and pos_end2+1<=pos_end1:
            y1 = pos_end2
            y2 = pos_end1
            ref_wanted1 = ref_wanted1[(y1-pos_start1):]
            seq_wanted1 = seq_wanted1[(y1-pos_start1):]
            insert_pos_anno1 = insert_pos_anno1[(y1-pos_start1):]
            m2,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand)
            m1,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand)
            Mut_infos = update_muts(Mut_infos,mut_seq1,ref_seq1,pos_seq1,y1)
            Mut_infos = update_muts(Mut_infos,mut_seq2,ref_seq2,pos_seq2,x1)
        elif pos_end2<=pos_start1:
            y1 = pos_start1
            y2 = pos_end1
            m2,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand)
            m1,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand)
            Mut_infos = update_muts(Mut_infos,mut_seq1,ref_seq1,pos_seq1,y1)
            Mut_infos = update_muts(Mut_infos,mut_seq2,ref_seq2,pos_seq2,x1)
        elif pos_end2+1 > pos_end1:
            m1,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand)
            Mut_infos = update_muts(Mut_infos,mut_seq2,ref_seq2,pos_seq2,x1)
            m2 = ''
            y1 = ''
            y2 = ''
    x1 = str(x1)
    x2 = str(x2)
    y1 = str(y1)
    y2 = str(y2)
    return x1,x2,y1,y2,m1,m2,Mut_infos

#%% read files and processing each line
if __name__ == '__main__':
    import math
    outdata = open(outfile,'w')
    outdata.write('Chrom'+'\t'+\
              'Strand'+'\t'+\
              'R1Start'+'\t'+\
              'R1End'+'\t'+\
              'R2Start'+'\t'+\
              'R2End'+'\t'+\
              'R1MeCpG'+'\t'+\
              'R2MeCpG'+'\t'+\
              'QualR1' +'\t'+\
              'QualR2' + '\t'+\
              'Mut'+'\t'+'ReadID'+'\tReadType\n')#header
    if infile[-4:] == '.sam':
        indata = open(infile,'r')
    elif infile[-4:] == '.bam':
        import pysam
        indata = pysam.AlignmentFile(infile,'rb')
    elif infile[-7:] == '.sam.gz':
        import gzip
        indata = gzip.open(infile,'rb')
    
    Reads_id1 = 'xl.Z'
    Reads_id2 = 'xl.Z'
    chrom1 = 'Unknown'
    chrom2 = 'Unknown'
    previous_out_signal1 = 1
    previous_out_signal2 = 1
    reads_num_count = 0
    while True:
        reads_num_count += 1
        if int(math.ceil(reads_num_count/10000))*10000 == reads_num_count:
            print('We have processed '+str(reads_num_count)+ ' reads...')
        if reads_num_count< -6:##just test
            break
        try:
            line = indata.next()
            if infile[-4:] == '.bam':
                line = line.to_string()
        except:
            break
        line_split1 = line.strip().split('\t')
        read_name = line_split1[0]
        chrom_name = line_split1[2]
        flags_score = flags_convert(int(line_split1[1]))
        pos = int(line_split1[3])
        map_quality = line_split1[4]
        cigar_info = line_split1[5]
        cigar_len,cigar_symbol = cigar_convert(cigar_info)
        seq_reads = line_split1[9]
        mut_anno = line_split1[12].split(':')[2]
        if flags_score[6] == 1: #paired reads
            #%% mate1
            if flags_score[1]==1: #mate 1
                if previous_out_signal1 == 1 and previous_out_signal2 == 1:
                    #new_mate1_process
                    previous_out_signal1 = 0
                    Reads_id1 = read_name
                    chrom1 = chrom_name
                    mq1 = map_quality
                    if flags_score[2]==1 and flags_score[3]==0: #molecule +
                        strand1 = '+'
                    else:
                        strand1 = '-'
                    seq_wanted1,insert_pos_anno1,ref_wanted1,pos_start1,pos_end1,Mut_infos1 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
                elif previous_out_signal1 == 0:
                    cpg_methy1,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand1)
                    Mut_infos1 = update_muts(Mut_infos1,mut_seq1,ref_seq1,pos_seq1,pos_start1)
                    str_mut_info = ''
                    if len(Mut_infos1) >=1:
                        for mut_info in Mut_infos1.keys():
                            str_mut_info += mut_info
                            str_mut_info += '|'
                    if len(cpg_methy1)>=2:
                        outdata.write(chrom1+'\t'+\
                                      strand1+'\t'+\
                                      str(pos_start1)+'\t'+\
                                      str(pos_end1)+'\t'+\
                                      '\t'+\
                                      '\t'+\
                                      cpg_methy1+'\t'+\
                                      '\t'+\
                                      mq1 +'\t'+\
                                      '\t'+\
                                      str_mut_info[:-1]+'\t'+Reads_id1+'\tL\n') #L:lost the other mate
                    previous_out_signal1 = 1
                    #new_mate1_process
                    previous_out_signal1 = 0
                    Reads_id1 = read_name
                    chrom1 = chrom_name
                    mq1 = map_quality
                    if flags_score[2]==1 and flags_score[3]==0: #molecule +
                        strand1 = '+'
                    else:
                        strand1 = '-'
                    seq_wanted1,insert_pos_anno1,ref_wanted1,pos_start1,pos_end1,Mut_infos1 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
                elif previous_out_signal2 == 0:
                    if Reads_id2 == read_name:#paired reads
                        previous_out_signal1 = 1
                        previous_out_signal2 = 1
                        #new mate1 process
                        chrom1 = chrom_name
                        mq1 = map_quality
                        if flags_score[2]==1 and flags_score[3]==0: #molecule +
                            strand1 = '+'
                        else:
                            strand1 = '-'
                        if chrom1 == chrom2 and strand1 == strand2:
                            seq_wanted1,insert_pos_anno1,ref_wanted1,pos_start1,pos_end1,Mut_infos1 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
                            x1,x2,y1,y2,m1,m2,Mut_infos = fine_deoverlap_mates(pos_start1,pos_end1,pos_start2,pos_end2,strand1,ref_wanted1,seq_wanted1,insert_pos_anno1,ref_wanted2,seq_wanted2,insert_pos_anno2,Mut_infos1,Mut_infos2)
                            str_mut_info = ''
                            if len(Mut_infos) >=1:
                                for mut_info in Mut_infos.keys():
                                    str_mut_info += mut_info
                                    str_mut_info += '|'
                            if len(m1)+len(m2)>=2:
                                outdata.write(chrom1+'\t'+\
                                              strand1+'\t'+\
                                              x1+'\t'+\
                                              x2+'\t'+\
                                              y1+'\t'+\
                                              y2+'\t'+\
                                              m1+'\t'+\
                                              m2+'\t'+\
                                              mq1 +'\t'+\
                                              mq2 + '\t'+\
                                              str_mut_info[:-1]+'\t'+Reads_id2+'\tP\n')#paired
                        else:
                            print('Inconcordant strand or chrom! '+Reads_id2)
                            pass
                    else:
                        cpg_methy2,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand2)
                        Mut_infos2 = update_muts(Mut_infos2,mut_seq2,ref_seq2,pos_seq2,pos_start2)
                        str_mut_info = ''
                        if len(Mut_infos2) >=1:
                            for mut_info in Mut_infos2.keys():
                                str_mut_info += mut_info
                                str_mut_info += '|'
                        if len(cpg_methy2)>=2:
                            outdata.write(chrom2+'\t'+\
                                          strand2+'\t'+\
                                          str(pos_start2)+'\t'+\
                                          str(pos_end2)+'\t'+\
                                          '\t'+\
                                          '\t'+\
                                          cpg_methy2+'\t'+\
                                          '\t'+\
                                          mq2 +'\t'+\
                                          '\t'+\
                                          str_mut_info[:-1]+'\t'+Reads_id2+'\tL\n')#lost
                        previous_out_signal2 = 1
                        #new_mate1_process
                        previous_out_signal1 = 0
                        Reads_id1 = read_name
                        chrom1 = chrom_name
                        mq1 = map_quality
                        if flags_score[2]==1 and flags_score[3]==0: #molecule +
                            strand1 = '+'
                        else:
                            strand1 = '-'
                        seq_wanted1,insert_pos_anno1,ref_wanted1,pos_start1,pos_end1,Mut_infos1 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
            else:
                if previous_out_signal1 == 1 and previous_out_signal2 == 1:
                    #new_mate2_process
                    previous_out_signal2 = 0
                    Reads_id2 = read_name
                    chrom2 = chrom_name
                    mq2 = map_quality
                    if flags_score[2]==0 and flags_score[3]==1: #molecule +
                        strand2 = '+'
                    else:
                        strand2 = '-'
                    seq_wanted2,insert_pos_anno2,ref_wanted2,pos_start2,pos_end2,Mut_infos2 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
                elif previous_out_signal2 == 0:
                    cpg_methy2,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand2)
                    Mut_infos2 = update_muts(Mut_infos1,mut_seq2,ref_seq2,pos_seq2,pos_start2)
                    str_mut_info = ''
                    if len(Mut_infos2) >=1:
                        for mut_info in Mut_infos2.keys():
                            str_mut_info += mut_info
                            str_mut_info += '|'
                    if len(cpg_methy2)>=2:
                        outdata.write(chrom2+'\t'+\
                                      strand2+'\t'+\
                                      str(pos_start2)+'\t'+\
                                      str(pos_end2)+'\t'+\
                                      '\t'+\
                                      '\t'+\
                                      cpg_methy2+'\t'+\
                                      '\t'+\
                                      mq2 +'\t'+\
                                      '\t'+\
                                      str_mut_info[:-1]+'\t'+Reads_id2+'\tL\n')
                    previous_out_signal2 = 1
                    #new_mate1_process
                    previous_out_signal2 = 0
                    Reads_id2 = read_name
                    chrom2 = chrom_name
                    mq2 = map_quality
                    if flags_score[2]==0 and flags_score[3]==1: #molecule +
                        strand2 = '+'
                    else:
                        strand2 = '-'
                    seq_wanted2,insert_pos_anno2,ref_wanted2,pos_start2,pos_end2,Mut_infos2 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
                elif previous_out_signal1 == 0:
                    if Reads_id1 == read_name:#paired reads
                        previous_out_signal1 = 1
                        previous_out_signal2 = 1
                        #new mate1 process
                        chrom2 = chrom_name
                        mq2 = map_quality
                        if flags_score[2]==0 and flags_score[3]==1: #molecule +
                            strand2 = '+'
                        else:
                            strand2 = '-'
                        if chrom1 == chrom2 and strand1 == strand2:
                            seq_wanted2,insert_pos_anno2,ref_wanted2,pos_start2,pos_end2,Mut_infos2 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
                            x1,x2,y1,y2,m1,m2,Mut_infos = fine_deoverlap_mates(pos_start1,pos_end1,pos_start2,pos_end2,strand1,ref_wanted1,seq_wanted1,insert_pos_anno1,ref_wanted2,seq_wanted2,insert_pos_anno2,Mut_infos1,Mut_infos2)
                            str_mut_info = ''
                            if len(Mut_infos) >=1:
                                for mut_info in Mut_infos.keys():
                                    str_mut_info += mut_info
                                    str_mut_info += '|'
                            if len(m1)+len(m2) >=2:
                                outdata.write(chrom1+'\t'+\
                                              strand1+'\t'+\
                                              x1+'\t'+\
                                              x2+'\t'+\
                                              y1+'\t'+\
                                              y2+'\t'+\
                                              m1+'\t'+\
                                              m2+'\t'+\
                                              mq1 +'\t'+\
                                              mq2 +'\t'+\
                                              str_mut_info[:-1]+'\t'+Reads_id1+'\tP\n')
                        else:
                            print('Inconcordant strand or chrom! '+Reads_id1)
                            pass
                    else:
                        cpg_methy1,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand1)
                        Mut_infos1 = update_muts(Mut_infos1,mut_seq1,ref_seq1,pos_seq1,pos_start1)
                        str_mut_info = ''
                        if len(Mut_infos1) >=1:
                            for mut_info in Mut_infos1.keys():
                                str_mut_info += mut_info
                                str_mut_info += '|'
                        if len(cpg_methy1)>=2:
                            outdata.write(chrom1+'\t'+\
                                          strand1+'\t'+\
                                          str(pos_start1)+'\t'+\
                                          str(pos_end1)+'\t'+\
                                          '\t'+\
                                          '\t'+\
                                          cpg_methy1+'\t'+\
                                          '\t'+\
                                          mq1 +'\t'+\
                                          '\t'+\
                                          str_mut_info[:-1]+'\t'+Reads_id1+'\tL\n')
                        previous_out_signal1 = 1
                        #new_mate1_process
                        previous_out_signal2 = 0
                        Reads_id2 = read_name
                        chrom2 = chrom_name
                        mq2 = map_quality
                        if flags_score[2]==0 and flags_score[3]==1: #molecule +
                            strand2 = '+'
                        else:
                            strand2 = '-'
                        seq_wanted2,insert_pos_anno2,ref_wanted2,pos_start2,pos_end2,Mut_infos2 = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
        else:
            if previous_out_signal1 == 0:
                cpg_methy1,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand1)
                Mut_infos1 = update_muts(Mut_infos1,mut_seq1,ref_seq1,pos_seq1,pos_start1)
                str_mut_info = ''
                if len(Mut_infos1) >=1:
                    for mut_info in Mut_infos1.keys():
                        str_mut_info += mut_info
                        str_mut_info += '|'
                if len(cpg_methy1)>=2:
                    outdata.write(chrom1+'\t'+\
                                  strand1+'\t'+\
                                  str(pos_start1)+'\t'+\
                                  str(pos_end1)+'\t'+\
                                  '\t'+\
                                  '\t'+\
                                  cpg_methy1+'\t'+\
                                  '\t'+\
                                  mq1 +'\t'+\
                                  '\t'+\
                                  str_mut_info[:-1]+'\t'+Reads_id1+'\tL\n')
                previous_out_signal1 = 1
            elif previous_out_signal2 == 0:
                cpg_methy2,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand2)
                Mut_infos2 = update_muts(Mut_infos1,mut_seq2,ref_seq2,pos_seq2,pos_start2)
                str_mut_info = ''
                if len(Mut_infos2) >=1:
                    for mut_info in Mut_infos2.keys():
                        str_mut_info += mut_info
                        str_mut_info += '|'
                if len(cpg_methy2)>=2:
                    outdata.write(chrom2+'\t'+\
                                  strand2+'\t'+\
                                  str(pos_start2)+'\t'+\
                                  str(pos_end2)+'\t'+\
                                  '\t'+\
                                  '\t'+\
                                  cpg_methy2+'\t'+\
                                  '\t'+\
                                  mq2 +'\t'+\
                                  '\t'+\
                                  str_mut_info[:-1]+'\t'+Reads_id2+'\tL\n')
                previous_out_signal2 = 1                
            previous_out_signal1 == 1
            previous_out_signal2 == 1
            ##read single sequence reads
            Reads_id = read_name
            chrom = chrom_name
            mq = map_quality
            if flags_score[3]==0: #molecule +
                strand = '+'
            else:
                strand = '-'
            seq_wanted,insert_pos_anno,ref_wanted,pos_start,pos_end,Mut_infos = get_mate_info(seq_reads,cigar_len,cigar_symbol,mut_anno,pos)
            cpg_methy,mut_seq,ref_seq,pos_seq = find_methy_mut(ref_wanted,seq_wanted,insert_pos_anno,strand)
            Mut_infos = update_muts(Mut_infos,mut_seq,ref_seq,pos_seq,pos_start)
            str_mut_info = ''
            if len(Mut_infos) >=1:
                for mut_info in Mut_infos.keys():
                    str_mut_info += mut_info
                    str_mut_info += '|'
            if len(cpg_methy)>=2:
                outdata.write(chrom+'\t'+\
                              strand+'\t'+\
                              str(pos_start)+'\t'+\
                              str(pos_end)+'\t'+\
                              '\t'+\
                              '\t'+\
                              cpg_methy+'\t'+\
                              '\t'+\
                              mq +'\t'+\
                              '\t'+\
                              str_mut_info[:-1]+'\t'+Reads_id+'\tS\n') #single
    if previous_out_signal1 == 0:
        cpg_methy1,mut_seq1,ref_seq1,pos_seq1 = find_methy_mut(ref_wanted1,seq_wanted1,insert_pos_anno1,strand1)
        Mut_infos1 = update_muts(Mut_infos1,mut_seq1,ref_seq1,pos_seq1,pos_start1)
        str_mut_info = ''
        if len(Mut_infos1) >=1:
            for mut_info in Mut_infos1.keys():
                str_mut_info += mut_info
                str_mut_info += '|'
        if len(cpg_methy1)>=2:
            outdata.write(chrom1+'\t'+\
                          strand1+'\t'+\
                          str(pos_start1)+'\t'+\
                          str(pos_end1)+'\t'+\
                          '\t'+\
                          '\t'+\
                          cpg_methy1+'\t'+\
                          '\t'+\
                          mq1 +'\t'+\
                          '\t'+\
                          str_mut_info[:-1]+'\t'+Reads_id1+'\tL\n')
        previous_out_signal1 = 1
    elif previous_out_signal2 == 0:
        cpg_methy2,mut_seq2,ref_seq2,pos_seq2 = find_methy_mut(ref_wanted2,seq_wanted2,insert_pos_anno2,strand2)
        Mut_infos2 = update_muts(Mut_infos1,mut_seq2,ref_seq2,pos_seq2,pos_start2)
        str_mut_info = ''
        if len(Mut_infos2) >=1:
            for mut_info in Mut_infos2.keys():
                str_mut_info += mut_info
                str_mut_info += '|'
        if len(cpg_methy2)>=2:
            outdata.write(chrom2+'\t'+\
                          strand2+'\t'+\
                          str(pos_start2)+'\t'+\
                          str(pos_end2)+'\t'+\
                          '\t'+\
                          '\t'+\
                          cpg_methy2+'\t'+\
                          '\t'+\
                          mq2 +'\t'+\
                          '\t'+\
                          str_mut_info[:-1]+'\t'+Reads_id2+'\tL\n')
    outdata.close()
    indata.close()


