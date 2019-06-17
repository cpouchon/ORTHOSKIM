#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
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

print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("gene_name","#taxa","maxlen","#1.0","#0.75","#0.5","#0.25")

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
        seqs.append(str(record.seq))

    reflen=len(max(seqs))

    mseqs=[]
    for elem in seqs:
        match=elem.replace("N","")
        mseqs.append(float(len(match))/float(reflen))

    full=sum(i >= 1.0 for i in mseqs)
    perct25=sum(i >= 0.25 for i in mseqs)
    perct50=sum(i >= 0.50 for i in mseqs)
    perct75= sum(i >= 0.75 for i in mseqs)

    print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (str(gene_name),int(taxa),int(reflen),int(full),int(perct75),int(perct50),int(perct25))
