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


parser = argparse.ArgumentParser(description='Selection of mito-nucrdna contigs from blast outputs')
parser.add_argument("--rdna", help="blast out tab for nucrdna matches",
                    type=str)
parser.add_argument("--mito", help="blast out tab for mitochondrion matches",
                    type=str)
parser.add_argument("--out_rdna", help="Out file for selected nucrdna contigs",
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



rdnatab=args.rdna
mitotab=args.mito
out_rdna_cont=args.out_rdna
out_mito_cont=args.out_mito

rdna_contigs={}
mito_contigs={}
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


select_rdna=list()
select_mito=list()

for cont in all_contigs.keys():
    if cont in mito_contigs.keys() and cont in rdna_contigs.keys():
        mito_aln=max(list(mito_contigs[cont].values()))
        rdna_aln=max(list(rdna_contigs[cont].values()))
        if rdna_aln>mito_aln:
            select_rdna.append(cont)
        elif rdna_aln<mito_aln:
            select_mito.append(cont)
        elif rdna_aln==mito_aln:
            pass
    elif cont in mito_contigs.keys() and cont not in rdna_contigs.keys():
        select_mito.append(cont)
    elif cont not in mito_contigs.keys() and cont in rdna_contigs.keys():
        select_rdna.append(cont)

if len(select_rdna)>0:
    with open(out_rdna_cont,"a+") as file:
        file.write("\n".join(select_rdna))
        file.close()
else:
    pass
if len(select_mito)>0:
    with open(out_mito_cont,"a+") as file:
        file.write("\n".join(select_mito))
        file.close()
else:
    pass
