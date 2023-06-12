"""
minjiezhang123@gmail.com    2023-03-16      python 3
This script can be used to calculate the RiboMethScore(ScoreC2) for small non-coding RNAs, especially for tRNAs.

tRNA molecules contain abundant modifications that will stop reverse transcription and generate heavy RT-stops. To avoid \
the false positive results resulted by RT-Stops caused by other modifications, the methylation score (MethScore) was \
computed with four neighboring nucleotides (+/-2 nucleotide window) from 5'-end and 3'-end reads at nucleotide resolution. \
The empirical rule with 1.5 standard deviation (u +/- 1.5*SD) was also used to filter the discrete values. If there are \
observations that fall outside the 1.5 standard deviation, only two neighboring nucleotides (+/-1 nucleotide window) from \
5'-end and 3'-end reads were used to calculate MethScore.

"""


#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
import sys, argparse, numpy, os, re, itertools, random
from datetime import datetime
from multiprocessing import Process, Lock, Manager
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np


if len(sys.argv) < 5:
    print('\n'+"Usage:  ")
    print("python  RiboMethScore_calling  5end.bedgraph  3end.bedgraph  gene.bed  5end/3end/both  outputprefix")
    print("")
    print("5end.bedgraph:   Input bedgraph file of 5end read coverage")
    print("3end.bedgraph:   Input bedgraph file of 3end read coverage")
    print("gene.bed:        gene bed files, 6 columns")
    print("5end/3end/both:  5end or 3end or both")
    print("                 5end: only use 5end read coverage to calculate the RiboMethScore")
    print("                 3end: only use 3end read coverage to calculate the RiboMethScore")
    print("                 both: use both 5end and 3end read coverage to calculate the RiboMethScore")
    print("outputprefix:    output prefix "+'\n')
    sys.exit()
    
bedgraph5end = sys.argv[1]
bedgraph3end = sys.argv[2]
genebed = sys.argv[3]
endcover = sys.argv[4]
outputprefix = sys.argv[5]
################################################################################


#2. function definitions
################################################################################
################################################################################
def timenow(): return str(datetime.now())[:-7]

##call RiboMethScore2(ScoreC2)    
def RiboMethScore2(l2,l1,x,r1,r2):
    #l2: reads count at position l2 (2nd left flanking regions of i)
    #l1: reads count at position l1 (1st left flanking regions of i)
    #x:  reads count at target position
    #r1: reads count at position r1 (1st right flanking regions of i)
    #r2: reads count at position r2 (2nd right flanking regions of i)
    if l2==0 and l1==0 and r1==0 and r2==0:
        ScoreC2 = 0
    else:
        ScoreC2 = 1-x/(0.5*( (l2*0.9+l1*1.0)/(1.0+0.9) + (r1*1.0+r2*0.9)/(1.0+0.9) ))
    if ScoreC2 > 0:
        return ScoreC2
    else:
        return 0
    

def RiboMethScore1(l1,x,r1):
    #l1: reads count at position l1 (1st left flanking regions of i)
    #x:  reads count at target position
    #r1: reads count at position r1 (1st right flanking regions of i)
    if l1==0 and r1==0:
        ScoreC1 = 0
    else:
        ScoreC1 = 1-x/(0.5*( (l1*1.0)/(1.0) + (r1*1.0)/(1.0) ))
    if ScoreC1 > 0:
        return ScoreC1
    else:
        return 0


def CallRiboMethforTRNA(l2,l1,x,r1,r2):
    list2 = [l2,l1,r1,r2]
    list_new = []
    # check the discrete value for list2
    outlier_ll = np.mean(list2) - 1.5*np.std(list2) # +/- 1SD: sigma 68.2%;  +/- 1.5SD: sigma 86.6%;  
    outlier_ul = np.mean(list2) + 1.5*np.std(list2) # +/- 1SD: sigma 68.2%;  +/- 1.5SD: sigma 86.6%;  
    for i in list2:
        if i <= outlier_ul: list_new.append(i)
        else: list_new.append('NA')
    if ('NA' in [list_new[0],list_new[3]]):   return RiboMethScore1(l1,x,r1)
    else:  return RiboMethScore2(l2,l1,x,r1,r2)
 

    
def getGenePos(genomepos, genestart, geneend, genestandard):
    genepos = ''
    if genestandard == "+":
        if int(genomepos) > int(geneend):       genepos = '+'+str(int(genomepos) - int(geneend))
        elif int(genomepos) >= int(genestart):  genepos = int(genomepos) - int(genestart) + 1
        elif int(genomepos) < int(genestart):   genepos = '-'+str(int(genestart) - int(genomepos))
    if genestandard == "-":
        if int(genomepos) > int(geneend):       genepos = '-'+str(int(genomepos) - int(geneend))
        elif int(genomepos) >= int(genestart):  genepos = int(geneend) - int(genomepos) + 1
        elif int(genomepos) < int(genestart):   genepos = '+'+str(int(genestart) - int(genomepos))
    return genepos
################################################################################



#3. Processing gene bed file
################################################################################
print(str(datetime.now())[:-7], "Reading gene bed file ...")
genebeddict = {} # {[gene],[gene_pos]}
end5_coverage = {};  end3_coverage = {}; 
genebedfile = open(genebed,'r')
# chr1    1   101 ACTB    1000    +
for line in genebedfile:    #(1 based)
    chr, start, end, genename, score, strand, biotype = line.rstrip('\n').split('\t')
    for i in range(int(start),int(end)+1,1):
        genebeddict[chr+'_'+str(i)] = [genename, getGenePos(i, start, end, strand)]
    for i in range(int(start)-10,int(end)+10,1):
        end5_coverage[chr+'_'+str(i)] = [0, strand]
        end3_coverage[chr+'_'+str(i)] = [0, strand]
genebedfile.close()


################ Process bam file to get 5'end #####################
if endcover == '5end' or endcover == 'both':
    print(str(datetime.now())[:-7], "Reading 5'end coverage bedgraph ...")
    # os.system("bedtools genomecov -ibam %s -5 -bg > %s" % (inputbam, outputprefix+'_5end.bedgraph'))
    # get 5end coverage
    MethySore_5end = open(bedgraph5end, 'r')
    # chr1  16520582    16520583    4  (0 based)
    for line in MethySore_5end:
        align = line.rstrip('\n').split('\t')
        chr, start, end, coverage = align[0], int(align[1]), int(align[2]), int(eval(align[3]))
        for pos in range(start+1, end+1, 1):
            if chr+'_'+str(pos) in end5_coverage:
                if end5_coverage[chr+'_'+str(pos)][1] == "+":
                    if chr+'_'+str(pos-1) in end5_coverage: end5_coverage[chr+'_'+str(pos-1)][0] += coverage
                if end5_coverage[chr+'_'+str(pos)][1] == "-":
                    if chr+'_'+str(pos+1) in end5_coverage: end5_coverage[chr+'_'+str(pos+1)][0] += coverage
    MethySore_5end.close()


################ Process bam file to get 3'end #####################
if endcover == '3end' or endcover == 'both':
    print(str(datetime.now())[:-7], "Reading 3'end coverage bedgraph ...")
    # os.system("bedtools genomecov -ibam %s -3 -bg > %s" % (inputbam, outputprefix+'_3end.bedgraph'))
    # get 3end coverage
    MethySore_3end = open(bedgraph3end, 'r')
    # chr1  16520582    16520583    4  (0 based)
    for line in MethySore_3end:
        align = line.rstrip('\n').split('\t')
        chr, start, end, coverage = align[0], int(align[1]), int(align[2]), int(eval(align[3]))
        for pos in range(start+1, end+1, 1):
            if chr+'_'+str(pos) in end3_coverage:
                end3_coverage[chr+'_'+str(pos)][0] += coverage
    MethySore_3end.close()


################ bedtools overlap aligned bedfile and genebed file #####################
output = open(outputprefix+'_'+str(endcover)+'_Nm.txt', 'w')
for i in genebeddict: #{chr_pos:[genename, genepos, strand]}
    chr = i.split('_')[0];  j = int(i.split('_')[1])
    genename = str(genebeddict[i][0]); genepos = str(genebeddict[i][1])
    pos,list_new = '',[]
    ## calcultae 5'end + 3'end
    l3,l2,l1,x,r1,r2,r3 = 0.01,0.01,0.01,0.01,0.01,0.01,0.01
    l2 = end5_coverage[chr+'_'+str(j-2)][0] + end3_coverage[chr+'_'+str(j-2)][0]
    l1 = end5_coverage[chr+'_'+str(j-1)][0] + end3_coverage[chr+'_'+str(j-1)][0]
    x  = end5_coverage[chr+'_'+str(j)][0]   + end3_coverage[chr+'_'+str(j)][0]
    r1 = end5_coverage[chr+'_'+str(j+1)][0] + end3_coverage[chr+'_'+str(j+1)][0]
    r2 = end5_coverage[chr+'_'+str(j+2)][0] + end3_coverage[chr+'_'+str(j+2)][0]
    #print(chr+'_'+str(j), l2,l1,r1,r2)
    if  numpy.mean([l1,x,r1]) >=20:
        if endcover == '5end':
            string = str(i)+'\t'+genename+'\t'+str(genepos)+'\t'+str(end5_coverage[chr+'_'+str(j)][0])+'\t'+str(CallRiboMethforTRNA(l2,l1,x,r1,r2))+'\n'
            output.write(string)
        if endcover == '3end':
            string = str(i)+'\t'+genename+'\t'+str(genepos)+'\t'+str(end3_coverage[chr+'_'+str(j)][0])+'\t'+str(CallRiboMethforTRNA(l2,l1,x,r1,r2))+'\n'
            output.write(string)
        if endcover == 'both':
            string = str(i)+'\t'+genename+'\t'+str(genepos)+'\t'+str(end5_coverage[chr+'_'+str(j)][0])+'\t'+str(end3_coverage[chr+'_'+str(j)][0])+'\t'+str(x)+'\t'+str(CallRiboMethforTRNA(l2,l1,x,r1,r2))+'\n'
            output.write(string)
    else:
        if endcover == '5end':
            string = str(i)+'\t'+genename+'\t'+str(genepos)+'\t'+str(end5_coverage[chr+'_'+str(j)][0])+'\t'+'NA\n'
            output.write(string)
        if endcover == '3end':
            string = str(i)+'\t'+genename+'\t'+str(genepos)+'\t'+str(end3_coverage[chr+'_'+str(j)][0])+'\t'+'NA\n'
            output.write(string)
        if endcover == 'both':
            string = str(i)+'\t'+genename+'\t'+str(genepos)+'\t'+str(end5_coverage[chr+'_'+str(j)][0])+'\t'+str(end3_coverage[chr+'_'+str(j)][0])+'\t'+str(x)+'\t'+'NA\n'
            output.write(string)

output.close()
#################################################################################



