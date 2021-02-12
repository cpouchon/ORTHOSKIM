#! /usr/bin/env python

import sys
import os, errno
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-t","--taxa", help="taxa list",
                    type=str)
parser.add_argument("-m","--mode", help="taxa list",
                    type=str,choices=["all","chloroplast","mitochondrion","nucrdna"])
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


path = args.path
model=args.mode

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    #taxa_list.append(l)
    taxa_dict.setdefault(l, [])

if model=="all":
    for taxa in taxa_dict.keys():
        try:
            with open(os.path.join(str(path+"chloroplast"), str(taxa+".cont_cpdna.log")), 'r') as out:
                chloro=out.read().splitlines()
            chloro_size=0
            chloro_count=len(chloro)
            cov=0
            for cont in chloro:
                if cont=="None":
                    continue
                else:
                    chloro_size=chloro_size+int(cont.split("length_")[1].split("_")[0])
                    cov=cov+float(cont.split("cov_")[1])*int(cont.split("length_")[1].split("_")[0])
            chloro_cov=round(cov/chloro_size,2)
        except:
            chloro_size=0
            chloro_count=0
            chloro_cov=0.0

        try:
            with open(os.path.join(str(path+"mitochondrion"), str(taxa+".cont_mtdna.log")), 'r') as out:
                mito=out.read().splitlines()
            mito_size=0
            mito_count=len(mito)
            cov=0
            for cont in mito:
                if cont=="None":
                    continue
                else:
                    mito_size=mito_size+int(cont.split("length_")[1].split("_")[0])
                    cov=cov+float(cont.split("cov_")[1])*int(cont.split("length_")[1].split("_")[0])
            mito_cov=round(cov/mito_size,2)
        except:
            mito_size=0
            mito_count=0
            mito_cov=0.0

        try:
            with open(os.path.join(str(path+"nucrdna"), str(taxa+".cont_nucrdna.log")), 'r') as out:
                rdna=out.read().splitlines()
            rdna_size=0
            rdna_count=len(rdna)
            cov=0
            for cont in rdna:
                if cont=="None":
                    continue
                else:
                    rdna_size=rdna_size+int(cont.split("length_")[1].split("_")[0])
                    cov=cov+float(cont.split("cov_")[1])*int(cont.split("length_")[1].split("_")[0])
            rdna_cov=round(cov/rdna_size,2)
        except:
            rdna_size=0
            rdna_count=0
            rdna_cov=0.0
        print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (str(taxa),int(chloro_count),int(chloro_size),float(chloro_cov),int(mito_count),int(mito_size),float(mito_cov),int(rdna_count),int(rdna_size),float(rdna_cov)))
elif model=="chloroplast":
    for taxa in taxa_dict.keys():
        try:
            with open(os.path.join(str(path+"chloroplast"), str(taxa+".cont_cpdna.log")), 'r') as out:
                chloro=out.read().splitlines()
            chloro_size=0
            chloro_count=len(chloro)
            cov=0
            for cont in chloro:
                if cont=="None":
                    continue
                else:
                    chloro_size=chloro_size+int(cont.split("length_")[1].split("_")[0])
                    cov=cov+float(cont.split("cov_")[1])*int(cont.split("length_")[1].split("_")[0])
            chloro_cov=round(cov/chloro_size,2)
        except:
            chloro_size=0
            chloro_count=0
            chloro_cov=0.0
        print ("%s\t%s\t%s\t%s" % (str(taxa),int(chloro_count),int(chloro_size),float(chloro_cov)))
elif model=="mitochondrion":
    for taxa in taxa_dict.keys():
            try:
                with open(os.path.join(str(path+"mitochondrion"), str(taxa+".cont_mtdna.log")), 'r') as out:
                    mito=out.read().splitlines()
                mito_size=0
                mito_count=len(mito)
                cov=0
                for cont in mito:
                    if cont=="None":
                        continue
                    else:
                        mito_size=mito_size+int(cont.split("length_")[1].split("_")[0])
                        cov=cov+float(cont.split("cov_")[1])*int(cont.split("length_")[1].split("_")[0])
                mito_cov=round(cov/mito_size,2)
            except:
                mito_size=0
                mito_count=0
                mito_cov=0.0
            print ("%s\t%s\t%s\t%s" % (str(taxa),int(mito_count),int(mito_size),float(mito_cov)))
elif model=="nucrdna":
    for taxa in taxa_dict.keys():
        try:
            with open(os.path.join(str(path+"nucrdna"), str(taxa+".cont_nucrdna.log")), 'r') as out:
                rdna=out.read().splitlines()
            rdna_size=0
            rdna_count=len(rdna)
            cov=0
            for cont in rdna:
                if cont=="None":
                    continue
                else:
                    rdna_size=rdna_size+int(cont.split("length_")[1].split("_")[0])
                    cov=cov+float(cont.split("cov_")[1])*int(cont.split("length_")[1].split("_")[0])
            rdna_cov=round(cov/rdna_size,2)
        except:
            rdna_size=0
            rdna_count=0
            rdna_cov=0.0
        print ("%s\t%s\t%s\t%s" % (str(taxa),int(rdna_count),int(rdna_size),float(rdna_cov)))
else:
    pass
