#!/usr/bin python2.7
# -*- coding: utf-8 -*-
# @Author: fengyq
# @Date:   2020-10-26 11:17:08
# @Last Modified by:   fengyq
# @Last Modified time: 2021-06-11 15:21:11
## Python version 2.7

# split a merged read of a read-pair by the linkers.
from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
from collections import Counter
__author__ = 'FYQ'

# split two reads of a read-pair by the linker. and rejoin by neighbor fragments

def fq(file):
    ## check if the file exist in the folder
    if not os.path.isfile(file):
        raise FileNotFoundError("File %s not found!"%file)
    ## read fastq file
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]  ## read in 4 lines in 1 cycle.

def split_qual_by_linker(Readin,linker1):
    # split a read by linker, r1 is list with 4 lines (a read record)
    # get the fragment locations seperated by linker
    r1_pos=[(0,0)] # add start_pos
    for match in re.finditer(linker1, Readin[1]):
        r1_pos.append(match.span())
    r1_pos.append((len(Readin[3]),len(Readin[3]))) # add end_pos

    # assign quality by the fragment positions
    r1_qual=[]
    for x in range(0,len(r1_pos)-1,1):
        f_start=r1_pos[x][1]
        f_end=r1_pos[x+1][0]
        r1_qual.append(Readin[3][f_start:f_end])
    return r1_qual

# split paired reads by linker
def combine_frag(read1, read2, read1_out, read2_out, statfile, frag_size, linker_seq):
    # Create new r1 FASTQs
    r1_splited = open(read1_out, 'w')
    # Create new r2 FASTQs
    r2_splited = open(read2_out, 'w')
    linker_seq = linker_seq # linker_seq = '.CGCGATATCTTATCTGACAG.|.CTGTCAGATAAGATATCGCG.'
    m= int(frag_size)
    # count the linkers in merged read
    r1_linker_count=[] # stats of linkers in read1
    r2_linker_count=[] # stats of linkers in read2
    r1_left_linker=0  # counts of dangling read1
    r2_left_linker=0  # counts of dangling read2
    L1=0 # sum of reads in fq1
    L2=0 # sum of reads in fq2
    newcount=0

    for r1,r2 in zip(fq(read1), fq(read2)):  # iteration Python 3's zip function is Python 2's izip
        #print(r1,r2)
        # count linker in the read1.
        r1_fragments=re.split(linker_seq, r1[1])
        L1=L1+1 # sum of all read
        r1_linker_count.append(len(r1_fragments)-1)

        # count linker in the read2.
        r2_fragments=re.split(linker_seq, r2[1])
        L2=L2+1 # sum of all read
        r2_linker_count.append(len(r2_fragments)-1)

        # split two reads by linker
        r1_qual=split_qual_by_linker(r1,linker_seq)
        r2_qual=split_qual_by_linker(r2,linker_seq)

        # Count the dangling linker
        if r1_qual[0]=="":
            r1_left_linker=r1_left_linker+1
        if r2_qual[0]=="":
            r2_left_linker=r2_left_linker+1

        # remove R1 framents<20nt, m is the shortest fragment to keep
        if len(r1_fragments) == len(r1_qual):
            r1_fragments = [elem for elem in r1_fragments if len(elem) >= m]
            r1_qual = [elem for elem in r1_qual if len(elem) >= m]
        else:
             print ("check the read: {}, frag_count is not equal to qual_count \n".format(r1[0]))

        # remove R2 framents<20nt, m is the shortest fragment to keep
        if len(r2_fragments) == len(r2_qual):
            r2_fragments = [elem for elem in r2_fragments if len(elem) >= m]
            r2_qual = [elem for elem in r2_qual if len(elem) >= m]
        else:
             print ("check the read2: {}, frag_count is not equal to qual_count \n".format(r2[0]))
             

        # Merge two list as one
        r1_fragments.extend(r2_fragments)
        r1_qual.extend(r2_qual)

        # Re_pair the fragments and write new paired-read one by one
        n=0 # counts of new read-pair
        for f1,f2,q1,q2 in zip(r1_fragments[:-1],r1_fragments[1:],r1_qual[:-1],r1_qual[1:]):
            n=n+1  # counts of new Frag-pair in a PE
            newname=re.sub('GWNJ.*th:','',r1[0]).split()[0] + "_split" + str(n) # rename the framents
            r1_L=[newname,f1.rstrip(), "+", q1.rstrip()]
            r2_R=[newname,f2.rstrip(), "+", q2.rstrip()]
            newcount=newcount+1 # sum of new Frag-pair
            for line in r1_L:
                r1_splited.write(line+"\n")
            for line2 in r2_R:
                r2_splited.write(line2+"\n")

    # statistic of linker_count and frequency
    linker_stat = open(statfile, 'w')
    linker_stat.write("linker_in_read1"+"\t"+"freq_in_fq1"+"\t"+"linker_in_read2"+"\t"+"freq_in_fq2"+"\n")
    #for key, value in Counter(linker_count).iteritems(): python2
    #for key, value in Counter(linker_count).items(): python3 
    stat1=Counter(r1_linker_count).keys() # linker number in_a_read in fq1
    stat2=Counter(r1_linker_count).values()  # linker freq in fq1
    stat3=Counter(r2_linker_count).keys()
    stat4=Counter(r2_linker_count).values()
    for a,b,c,d in zip(stat1,stat2,stat3,stat4):
        linker_stat.write(str(a) + "\t" + str(b) + "\t" + str(c) + "\t" + str(d) +"\n")

    # double check the read numbers and print out
    print (' read1 has a linker at 5 end, counts of read are {}.\n'.format(r1_left_linker))
    print (' read2 has a linker at 5 end, counts of read are {}.\n'.format(r2_left_linker))
    print (' Counts of input reads are fastq1:{} and fastq2:{}.\n'.format(L1,L2))
    print (' Counts of new splited output read-paires are {}.\n'.format(newcount))

    r1_splited.close()
    r2_splited.close()
    linker_stat.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-L', '--linker_seq', required=True,help="eg:.CGCGATATCTTATCTGACAG.|.CTGTCAGATAAGATATCGCG.")
    parser.add_argument('-S', '--statfile', required=True,help="write statfile of linkers")
    parser.add_argument('-M', '--frag_size', required=True,help="remove read with length < m")
    parser.add_argument('-r1', '--read1', required=True,help="read1 input, fq or fq.gz")
    parser.add_argument('-r2', '--read2', required=True,help="read2 input, fq or fq.gz")
    parser.add_argument('-R1', '--read1_out', required=True,help="Read1 output, PE1")
    parser.add_argument('-R2', '--read2_out', required=True,help="Read2 output, PE2")
    args = vars(parser.parse_args())
    combine_frag( args['read1'], args['read2'], args['read1_out'], args['read2_out'],args['statfile'], args['frag_size'], args['linker_seq'])

if __name__ == '__main__':
    main()
