#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
#import numpy
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


parser = argparse.ArgumentParser(description='')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-i","--infile", help="input file list of sequence files",
                    type=str)
parser.add_argument("-t","--taxa", help="list of taxa",
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

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    #taxa_list.append(l)
    taxa_dict.setdefault(l, [])


#print ("%s\t%s\t%s\t%s" % ("gene","taxa","size","refsize"))

for file in in_files:
    cf = file.rstrip()
    cfile_path = cf
    cfile_name = os.path.basename(cfile_path)
    cindex_of_dot = cfile_name.index('.')
    gene_name = cfile_name[:cindex_of_dot]

    #cur_genome = SeqIO.parse(file, "fasta")
    cur_genome = SeqIO.to_dict(SeqIO.parse(file, "fasta"))

    for taxa in taxa_dict:
        taxaname=str(taxa+";")
        if taxaname in cur_genome.keys():
            if cur_genome[taxaname].description.split("; ")[2].split("=")[1]=="exon":
                rsiz=1.0
                siz=int(cur_genome[taxaname].description.split("; ")[3].split("=")[1])
                type=str(cur_genome[taxaname].description.split("; ")[2].split("=")[1])
            else:
                rsiz=1.0
                siz=int(cur_genome[taxaname].description.split("; ")[3].split("=")[1])
                type=str(cur_genome[taxaname].description.split("; ")[2].split("=")[1])
        else:
            rsiz=0.0
            siz=0
            type="NA"
        print ("%s\t%s\t%s\t%s\t%s" % (str(gene_name),str(taxa),int(siz),float(rsiz),str(type)))
