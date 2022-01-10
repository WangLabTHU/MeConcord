
# MeConcord
MeConcord is a method used to investigate local read-level DNA methylation patterns for intermediately methylated regions with bisulfite sequencing data.

# 1.Obtaining CpG positions across genome
Usage: python pre_cpg_pos.py -i hg38.fa -o ./cpg_pos/
* i,  The path to reference sequences (.fa);
* o,  The path that you want to deposit the positions of CpG sites, each chromosome has a seperate file;
* h,  Help information

# 2.Converting mapped Bam, Sam, Sam.gz files from Bismark To methylation recordings read-by-read
Usage: python s1_bamToMeRecord.py -i test.bam -o test -c 0
* i,  The path to input files (.bam or .sam or .sam.gz);
* o,  Output prefix;
* c,  Clipping read ends with such base number (defalut 0); can be used when sequencing quality of read ends is not good. such as -c 5 to remove 5 bases from the both ends of the reads.
* h,  Help information

# 3.Spliting the big MeRecord files into small files of each chromosome to redude memory requirement in next step
Usage: python s2_RecordSplit.py -i ./test_ReadsMethyAndMuts.txt -o ./test -g chr1,chr2,chr3,chr4,chr5
* i,  The path to s1 output. ( end with _ReadsMethyAndMuts.txt);
* o,  Output prefix;
* g,  Chromosomes used; (default chromsome 1-22); chromosomes shoud be seperated by comma;
* h,  Help information

# 4. Calculating concordance metrics
Usage: python s3_RecordToMeConcord.py -p 4 -i ./test -o ./test -r ./region.bed -c ./cpgpos/ -b 150 -m 600 -z 0 -g chr1,chr2,chr3
* i,  The path to s2_RecordSplit.py output, with prefix file name;
* p,  Threads used for parallel computation; default is 4;
* o,  Output prefix;
* r,  The files with genomic regions for computation, chrom, start, end seperated by tab;
* c,  Cpg position folder, output of pre_cpg_pos.py;
* b,  Bin size (default 150bp);
* z,  Whether is the genomic file based on 0; 0 (default) or 1; output is same to input bins; if -r is a bed file, -z should be 1;
* g,  Chromosomes used; (default chromsome 1-22); chromosomes shoud be seperated by comma;
* m,  Maximum of fragement length in sequencing library(default 600bp for paired-end reads). if there are single-end reads,m should be set as the length of reads, if not sure, default will work for most cases;

# 5. Methylation recordings to methylation matrix (optional)
Usage: python s4_RecordToMeMatrix.py -i ./test -o ./test -r ./p1.bed -c ./cpgpos/ -m 600 -z 0 -g chr1,chr2
* i,  The path to s2_RecordSplit.py output, with prefix file name;
* o,  Output prefix;
* r,  The files with genomic regions for computation, chrom, start, end seperated by tab;
* c,  Cpg position folder, output of pre_cpg_pos.py;
* z,  Whether is the genomic file based on 0; 0 (default) or 1; output is same to input bins; if -r is a bed file, -z should be 1;
* g,  Chromosomes used; (default chromsome 1-22); chromosomes shoud be seperated by comma;
* m,  Maximum of reads length (default 600bp for paired-end reads). if there are single-end reads,m should be set length of reads, if not sure, default will work for most cases;

# 6. Visualization of methylation matrix (optional)
Usage: visualization_Matlab.m 
Open this script and edit
* `path_to_matrix` as the path you deposit the MeMatrix;
* `path_to_cpgPos` as the path you deposit CpG positions of genome, which is the result of pre_cpg_pos.py;
* `name` as the name of MeMatrix, for example 'test_chr1_1287967_1288117';
Output:two lollipop plots, one without considering distance between CpGs, one considering distance between CpGs.
*unmethylated CpGs are labeled as light blue
*CpGs without signal are labeled as grey
*methylated CpGs are labeled as dark red

# Test for an example
*python s1_bamToMeRecord.py -i ./GM12878_chr1_1286017_1294783.bam -o ./test/test -c 2

*python s2_RecordSplit.py -i ./test/test_ReadsMethyAndMuts.txt -o ./test/test -g chr1

*python s3_RecordToMeConcord.py -p 1 -i ./test/test -o ./test/test -r ./test/tmp1.bed -c ./test/ -b 150 -m 600 -z 1 -g chr1

*python s4_RecordToMeMatrix.py -i ./test/test -o ./test/test -r ./test/tmp2.bed -c ./test/ -m 600 -z 1 -g chr1