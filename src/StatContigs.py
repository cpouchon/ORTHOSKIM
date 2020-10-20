#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
#import numpy
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqFeature import FeatureLocation


parser = argparse.ArgumentParser(description='')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-trim", help="[mode] statistics made on trimmed alignments",
                    action="store_true")
parser.add_argument("-s","--size", help="taxa number in alignment",
                    type=int)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

taxa_size=args.size

path = args.path

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    #taxa_list.append(l)
    taxa_dict.setdefault(l, [])


for taxa in taxa_dict.keys():
        #chloro=list()
        try:
            with open(os.path.join(str(path+"chloroplast"), str(taxa+".contigs_chloro.info")), 'r') as out:
                chloro=out.read().splitlines()
            chloro_size=0
            chloro_count=len(chloro)
            cov=0
            for cont in chloro:
                chloro_size=chloro_size+int(cont.split("length_")[1].split("_")[0])
                cov=cov+float(cont.split("cov_")[1])*int(cont.split("length_")[1].split("_")[0])
            chloro_cov=round(cov/chloro_size,2)
        except:
            chloro_size=0
            chloro_count=0
            chloro_cov=0.0

        


    print ("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (str(gene_name),int(aln_size),int(mp1),int(mp2),int(mp3),int(mp4),int(info_sites)))
