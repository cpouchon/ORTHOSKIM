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


parser = argparse.ArgumentParser(description='Selection of organelles contigs from blast outputs')
parser.add_argument("--chloro", help="blast out tab for chloroplast matches",
                    type=str)
parser.add_argument("--mito", help="blast out tab for mitochondrion matches",
                    type=str)
parser.add_argument("--out_chloro", help="Out file for selected chloroplast contigs",
                    type=str)
parser.add_argument("--out_mito", help="Out file for selected mitochondrion contigs",
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



chlorotab=args.chloro
mitotab=args.mito
out_chloro_cont=args.out_chloro
out_mito_cont=args.out_mito

chloro_contigs={}
mito_contigs={}
all_contigs={}

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

with open(mitotab) as f:
    for l in f:
        contigid=l.rstrip().split("\t")[0]
        mito_contigs.setdefault(contigid, {})
        all_contigs.setdefault(contigid, {})

        refid=l.rstrip().split("\t")[1]
        lenmapp=int(l.rstrip().split("\t")[3])

        if refid in mito_contigs[contigid].keys():
            oldval=int(mito_contigs[contigid][refid])
            newval=oldval+lenmapp
            mito_contigs[contigid][refid]=newval
        else:
            mito_contigs[contigid][refid]=lenmapp


select_chloro=list()
select_mito=list()

for cont in all_contigs.keys():
    if cont in mito_contigs.keys() and cont in chloro_contigs.keys():
        mito_aln=max(list(mito_contigs[cont].values()))
        chloro_aln=max(list(chloro_contigs[cont].values()))
        if chloro_aln>mito_aln:
            select_chloro.append(cont)
        elif chloro_aln<mito_aln:
            select_mito.append(cont)
        elif chloro_aln==mito_aln:
            pass
    elif cont in mito_contigs.keys() and cont not in chloro_contigs.keys():
        select_mito.append(cont)
    elif cont not in mito_contigs.keys() and cont in chloro_contigs.keys():
        select_chloro.append(cont)


if len(select_chloro)>0:
    with open(out_chloro_cont,"a+") as file:
        file.write("\n".join(select_chloro))
        file.close()
else:
    pass
if len(select_mito)>0:
    with open(out_mito_cont,"a+") as file:
        file.write("\n".join(select_mito))
        file.close()
else:
    pass
