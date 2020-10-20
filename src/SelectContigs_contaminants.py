#!/usr/bin/env python

import glob
import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Identification of contaminant according to taxonomic assignment of contigs into rRNA database. Script was writen by C. Pouchon (2020).')
parser.add_argument("--blast", help="output from blast of contigs into database",
                    type=str)
parser.add_argument('--rank', help='expected taxonomic rank', type=str)
parser.add_argument("--taxo", help="NCBI accessions to taxonomy for rRNA databases IDs",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

blastout=args.blast
expected=args.rank
NCBI_accessions=args.taxo

contigs={}
toremove=list()
cont_to_remove=list()
accessions={}
access_dict={}
bitscore_dict={}
best_hit={}

with open(NCBI_accessions) as acessf:
    for l in acessf:
        tab=l.rstrip().split('\t')
        if 'accession' in tab:
            pass
        else:
            code=tab[0]
            code_txid=tab[2]
            access_dict[code]=str(code_txid)

with open(blastout) as blastf:
    for l in blastf:
        contigid=l.rstrip().split("\t")[0]
        contigs.setdefault(contigid, {})
        refid=l.rstrip().split("\t")[1]
        accessions.setdefault(refid, "NA")
        bitscore=float(l.rstrip().split("\t")[11])

        if contigid in contigs.keys():
            if refid in contigs[contigid].keys():
                old_score=float(contigs[contigid][refid])
                new_score=float(old_score)+float(bitscore)
                contigs[contigid][refid]=new_score
            else:
                contigs[contigid][refid]=bitscore
        else:
            pass

for allcont in contigs.keys():
    for elem in contigs[allcont].keys():
        if allcont in bitscore_dict.keys():
            score=contigs[allcont][elem]
            if float(score) >= float(bitscore_dict[allcont]):
                best_hit[allcont]=elem
                bitscore_dict[allcont]=float(score)
            else:
                pass
        else:
            score=contigs[allcont][elem]
            best_hit[allcont]=elem
            bitscore_dict[allcont]=float(score)

for k in accessions.keys():
    if k in access_dict.keys():
        tax=access_dict[k]
        if expected.lower() in tax.lower():
            accessions[k]="TRUE"
        else:
            accessions[k]="FALSE"

for cont in best_hit.keys():
    if cont in cont_to_remove:
        continue
    else:
        refcont=best_hit[cont]
        if refcont in accessions.keys():
            toprint=list()
            if accessions[refcont]=="FALSE":
                if cont in cont_to_remove:
                    pass
                else:
                    toprint.append(cont)
                    toprint.append(refcont)
                    toprint.append(access_dict[refcont])
                    toremove.append(toprint)
                    cont_to_remove.append(cont)
        else:
            pass

if len(toremove)>0:
    for l in toremove:
        print("\t".join(l))
else:
    print("none")
