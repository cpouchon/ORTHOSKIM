#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
import numpy
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


parser = argparse.ArgumentParser(description='')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-i","--infile", help="input file list of sequence files",
                    type=str)
parser.add_argument("-pfind","--pathfind", help="[mode] search all sequences files in given path (-p) otherwise parse a given list (-i)",
                    action="store_true")
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

if args.pathfind:
    'Search all fasta files files'
    path = args.path
    in_files = []
    for r, d, f in os.walk(path):
        for file in f:
            if '.fa' in file:
                in_files.append(os.path.join(r, file))
else:
    'We parse a list of files'
    input_list = args.infile
    with open(input_list) as f:
        in_files = f.readlines()

print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("gene","taxa","mean","min","max","std","pct25","pct50","pct75"))

for file in in_files:
    cf = file.rstrip()
    cfile_path = cf
    cfile_name = os.path.basename(cfile_path)
    cindex_of_dot = cfile_name.index('.')
    gene_name = cfile_name[:cindex_of_dot]

    cur_genome = SeqIO.parse(file, "fasta")

    taxa = 0
    'keep the seq object into new list'
    seqs = []
    for record in cur_genome:
        #info = record.description.split(" ")[1:len(record.description.split(" "))]
        #dicinfo = dict(item.split("=") for item in info)
        #dicinfo['seq']=str(record.seq)
        taxa = taxa+1
        seqs.append(len(record.seq))

    maxl=max(seqs)
    minl=min(seqs)
    mean=round(numpy.mean(seqs),2)
    std=round(numpy.std(seqs),2)
    pct25=round(numpy.percentile(seqs,25),2)
    pct50=round(numpy.percentile(seqs,50),2)
    pct75=round(numpy.percentile(seqs,75),2)
    print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (str(gene_name),int(taxa),int(mean),int(minl),int(maxl),int(std),int(pct25),int(pct50),int(pct75)))
