#!/usr/bin/env python


import glob
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqFeature import FeatureLocation
import sys
import os, errno
import argparse



parser = argparse.ArgumentParser(description='Trim missing/gappy position in alignment according to a given thresold')
parser.add_argument("-i","--infile", help="Fasta input file of alignment",
                    type=str)
parser.add_argument("-t","--thresold", help="minimal thresold value allowed for missing/gappy position (1.0 means to position with 100% of missing/gapy Nt were trimmed)",
                    type=float)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

file=args.infile
thr=args.thresold

align=AlignIO.read(file,"fasta")
len_seq=align.get_alignment_length()

'identify missing/gappy position'
trimpos=[]
for Nt in range(len_seq):
    count_N=align[:,Nt].count('N')
    count_gap=align[:,Nt].count('-')
    Nscore=(float(count_N)/float(len(align)))*100
    Gscore=(float(count_gap)/float(len(align)))*100
    if Nscore >= thr or Gscore >= thr:
        trimpos.append(int(Nt))

poslist=set(range(len_seq))-set(trimpos)

'extraction of non-missing/gappy positions'
for seq in range(len(align)):
    newseq="".join([align[seq].seq[x] for x in poslist])
    header=align[seq].description
    print (">"+header)
    print (newseq)
