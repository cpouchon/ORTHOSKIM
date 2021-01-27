#! /usr/bin/env python

import sys
import os, errno
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-t","--taxa", help="taxa list",
                    type=str)
parser.add_argument("-d","--database", help="path to database files (RFAM, SILVA and DBFAM)",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


path = args.path
outdirect = args.outdir
dbpath= args.database
giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    #taxa_list.append(l)
    taxa_dict.setdefault(l, [])


in_files = []
for r, d, f in os.walk(dbpath):
    for file in f:
        if 'DBFAM_' in file or 'rfam-' in file or 'silva-' in file:
            if os.path.join(r, file) not in in_files:
                in_files.append(os.path.join(r,file))
keys_ref={}
for f in in_files:
    if "rfam" in f:
        type="RFAM"
    elif "silva" in f:
        type="SILVA"
    elif "chloroplast" in f:
        type="DBFAM_chloroplast"
    elif "mitochondrion" in f:
        type="DBFAM_mitochondrion"
    elif "nucrdna" in f:
        type="DBFAM_nucrdna"
    with open(f,'r',encoding='utf-8', errors='ignore') as ffile:
        for line in ffile:
            if line.startswith(">"):
                id=line.rstrip().replace(">","").split(" ")[0]
                keys_ref[id]=str(type)

for taxa in taxa_dict.keys():
        #chloro=list()
        try:
            with open(os.path.join(str(path),str("bad_"+taxa+".log")), 'r') as out:
                contam=out.read().splitlines()
            chloro_size=0
            mito_size=0
            nucrdna_size=0
            rfam_size=0
            silva_size=0
            for cont in contam:
                line=cont.split("\t")
                seltax=",".join(line[2].split(",")[0:4])
                size=int(line[0].split("length_")[1].split("_")[0])
                contid=line[1]
                if contid in keys_ref:
                    conttype=keys_ref[contid]
                    toprint=list()
                    toprint.append(taxa)
                    toprint.append(conttype)
                    toprint.append(str(size))
                    toprint.append(seltax)
                    with open(os.path.join(str(outdirect),'contaminant_full_statistics.txt'), 'a+') as full:
                        full.write("\t".join(toprint)+'\n')
                        full.close()
                    if "RFAM" in conttype:
                        rfam_size=rfam_size+size
                    elif "SILVA" in conttype:
                        silva_size=silva_size+size
                    elif "DBFAM_nucrdna" in conttype:
                        nucrdna_size=nucrdna_size+size
                    elif "DBFAM_chloroplast" in conttype:
                        chloro_size=chloro_size+size
                    elif "DBFAM_mitochondrion" in conttype:
                        mito_size=mito_size+size
                else:
                    pass

        except:
            chloro_size=0
            mito_size=0
            nucrdna_size=0
            rfam_size=0
            silva_size=0

        toprint=list()
        toprint.append(taxa)
        toprint.append(str(chloro_size))
        toprint.append(str(mito_size))
        toprint.append(str(nucrdna_size))
        toprint.append(str(rfam_size))
        toprint.append(str(silva_size))
        with open(os.path.join(str(outdirect),'contaminant_statistics.txt'), 'a+') as out:
            out.write("\t".join(toprint)+'\n')
            out.close()
