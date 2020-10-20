#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
#from Bio.Alphabet import *
from Bio.SeqRecord import *
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from ete3 import NCBITaxa

import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Selection of chloro-nucrdna contigs from blast outputs')
parser.add_argument("--rdna", help="blast out tab for nucrdna matches",
                    type=str)
parser.add_argument("--chloro", help="blast out tab for chloroplast matches",
                    type=str)
parser.add_argument("--out_rdna", help="Out file for selected nucrdna contigs",
                    type=str)
parser.add_argument("--out_chloro", help="Out file for selected chloroplastic contigs",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


def mkdir(path, overwrite=False):
    '''
    function to create a directory for output fasta
    '''
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not overwrite:
                print ("path '%s' already exists" % path)   # overwrite == False and we've hit a directory that exists
        else: raise



rdnatab=args.rdna
chlorotab=args.chloro
out_rdna_cont=args.out_rdna
out_chloro_cont=args.out_chloro

rdna_contigs={}
chloro_contigs={}
all_contigs={}

with open(rdnatab) as f:
    for l in f:
        contigid=l.rstrip().split("\t")[0]
        rdna_contigs.setdefault(contigid, {})
        all_contigs.setdefault(contigid, {})

        refid=l.rstrip().split("\t")[1]
        lenmapp=int(l.rstrip().split("\t")[3])

        if refid in rdna_contigs[contigid].keys():
            oldval=int(rdna_contigs[contigid][refid])
            newval=oldval+lenmapp
            rdna_contigs[contigid][refid]=newval
        else:
            rdna_contigs[contigid][refid]=lenmapp

with open(chlorotab) as f:
    for l in f:
        contigid=l.rstrip().split("\t")[0]
        chloro_contigs.setdefault(contigid, {})
        all_contigs.setdefault(contigid, {})

        refid=l.rstrip().split("\t")[1]
        lenmapp=int(l.rstrip().split("\t")[3])

        if refid in chloro_contigs[contigid].keys():
            oldval=int(chloro_contigs[contigid][refid])
            newval=oldval+lenmapp
            chloro_contigs[contigid][refid]=newval
        else:
            chloro_contigs[contigid][refid]=lenmapp


select_rdna=list()
select_chloro=list()

for cont in all_contigs.keys():
    if cont in chloro_contigs.keys() and cont in rdna_contigs.keys():
        chloro_aln=max(list(chloro_contigs[cont].values()))
        rdna_aln=max(list(rdna_contigs[cont].values()))
        if rdna_aln>chloro_aln:
            select_rdna.append(cont)
        elif rdna_aln<chloro_aln:
            select_chloro.append(cont)
        elif rdna_aln==chloro_aln:
            pass
    elif cont in chloro_contigs.keys() and cont not in rdna_contigs.keys():
        select_chloro.append(cont)
    elif cont not in chloro_contigs.keys() and cont in rdna_contigs.keys():
        select_rdna.append(cont)

if len(select_rdna)>0:
    with open(out_rdna_cont,"a+") as file:
        file.write("\n".join(select_rdna))
        file.close()
else:
    pass
if len(select_chloro)>0:
    with open(out_chloro_cont,"a+") as file:
        file.write("\n".join(select_chloro))
        file.close()
else:
    pass
