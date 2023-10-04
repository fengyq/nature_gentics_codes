#!/usr/bin python2.7
# -*- coding: utf-8 -*-
# @Author: fengyq
# @Date:   2020-10-26 11:17:08
# @Last Modified by:   fengpku
# @Last Modified time: 2021-10-03 12:17:42
## Python version 2.7

from __future__ import print_function
import os
import re
import gzip
import argparse
from collections import Counter
__author__ = 'FYQ'

# split a merged read by the linkers.

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

# split merged read by linker
def linker_remove(linker_seq, read1, read0_out, read1_out, read2_out, statfile, frag_size):

    # Create new r1 FASTQ
    r1_splited = open(read1_out, 'w')
    # Create new r2 FASTQ
    r2_splited = open(read2_out, 'w')
    # Create non-useful single read FASTQ
    r0_nouse = open(read0_out, 'w')
    linker_seq = linker_seq # linker_seq = '.CGCGATATCTTATCTGACAG.|.CTGTCAGATAAGATATCGCG.'
    m= int(frag_size)

    # count the linkers in merged read
    linker_count=[]
    left_linker=0
    right_linker=0
    L1=0 # input read counts
    L2=0 # output read pairs
    L3=0 # nouse read counts

    for r1 in fq(read1):  ## iteration using the predefined function fq()

        #split read by the linker
        r1_fragments=re.split(linker_seq,r1[1])
        L1=L1+1 # sum of all read

        # get the fragment locations seperated by linker
        r1_pos=[(0,0)] # add start_pos
        for match in re.finditer(linker_seq, r1[1]):
            r1_pos.append(match.span())
        r1_pos.append((len(r1[3]),len(r1[3]))) # add end_pos

        # assign quality by the fragment positions
        r1_qual=[]
        for x in range(0,len(r1_pos)-1,1):
            f_start=r1_pos[x][1]
            f_end=r1_pos[x+1][0]
            r1_qual.append(r1[3][f_start:f_end])

        # count linkers in the read
        linker_count.append(len(r1_fragments)-1)
        # Count the frament types
        if r1_qual[0]=="":
            left_linker=left_linker+1
        if r1_qual[-1]=="":
            right_linker=right_linker+1

        # remove framents<20nt, m is the shortest fragment to keep
        if len(r1_fragments) == len(r1_qual):
            r1_fragments = [elem for elem in r1_fragments if len(elem) >= m]
            r1_qual = [elem for elem in r1_qual if len(elem) >= m]
        else:
             print ("check the read: {}, frag_count is not equal to qual_count \n".format(r1[0]))
             
        ## Write read into files
        #1 the merged read contain no fragements>20nt or no linker
        if len(r1_fragments)<2:
            L3=L3+1
            for line in r1:
                r0_nouse.write(line)

        #2 if multiple fragmen >20nt in the splited read
        elif len(r1_fragments)>1 and len(r1_fragments)==len(r1_qual):
            L2=L2+len(r1_fragments)-1  # sum of all new splited read
            n=0
            # Report the fragments one by one
            for f1,f2,q1,q2 in zip(r1_fragments[:-1],r1_fragments[1:],r1_qual[:-1],r1_qual[1:]):
                n=n+1
                newname=re.sub('GWNJ.*th:','',r1[0]).split()[0] + "_split" + str(n)  # rename the framents
                r1_L=[newname,f1,"+",q1]
                r2_R=[newname,f2.rstrip(), "+", q2.rstrip()]
                for line in r1_L:
                    r1_splited.write(line+"\n")
                for line2 in r2_R:
                    r2_splited.write(line2+"\n")

        # other conditions
        else:
             print ("check the read {}, something is wrong \n".format(r1[0]))

    # statistic of linker_count
    # Create a file to write statistic
    linker_stat = open(statfile, 'w')
    linker_stat.write("linker_in_a_read"+"\t"+"read_count"+"\n")
    #for key, value in Counter(linker_count).iteritems(): python2
    #for key, value in Counter(linker_count).items(): python3
    for key, value in zip(Counter(linker_count).keys(),Counter(linker_count).values()):
        linker_stat.write(str(key) + "\t" + str(value) + "\n")
    #Counter(linker_count).keys() # equals to list(set(words))
    #Counter(linker_count).values() # counts the elements' frequency

    # double check the read numbers
    print (' Merged read has a linker at 5 end, counts of read are {}.\n'.format(left_linker))
    print (' Merged read has a linker at 3 end, counts of read are {}.\n'.format(right_linker))
    print (' Counts of input reads are {}.\n'.format(L1))
    print (' Counts of output read-paires are {}.\n'.format(L2))
    print (' Counts of no-linker merged reads are {}.\n'.format(L3))

    r1_splited.close()
    r2_splited.close()
    r0_nouse.close()
    linker_stat.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-L', '--linker_seq', required=True,help="eg:.CGCGATATCTTATCTGACAG.|.CTGTCAGATAAGATATCGCG.")
    parser.add_argument('-S', '--statfile', required=True,help="write statfile of linkers")
    parser.add_argument('-M', '--frag_size', required=True,help="remove read with length < m")
    parser.add_argument('-r1', '--read1', required=True,help="read1 input, fq or fq.gz")
    parser.add_argument('-R0', '--read0_out', required=True,help="Read0 output, no linker")
    parser.add_argument('-R1', '--read1_out', required=True,help="Read1 output, PE1")
    parser.add_argument('-R2', '--read2_out', required=True,help="Read2 output, PE2")
    args = vars(parser.parse_args())
    linker_remove(args['linker_seq'], args['read1'], args['read0_out'],args['read1_out'], args['read2_out'],args['statfile'], args['frag_size'])


if __name__ == '__main__':
    main()
