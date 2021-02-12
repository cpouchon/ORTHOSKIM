#! /usr/bin/env python

import sys
import os, errno
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-t","--taxa", help="taxa list",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


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
            with open(os.path.join(str(path+"chloroplast"), str(taxa+".cont_cpdna.log")), 'r') as out:
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

        print ("%s\t%s\t%s\t%s" % (str(taxa),int(chloro_count),int(chloro_size),float(chloro_cov)))
