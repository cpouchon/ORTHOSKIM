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

if args.trim:
    'Search all fasta files files'
    path = args.path
    in_files = []
    for r, d, f in os.walk(path):
        for file in f:
            if '.trim' in file:
                in_files.append(os.path.join(r, file))
else:
        'Search all fasta files files'
        path = args.path
        in_files = []
        for r, d, f in os.walk(path):
            for file in f:
                if '.fa' in file:
                    in_files.append(os.path.join(r, file))




for file in in_files:
    cf = file.rstrip()
    cfile_path = cf
    cfile_name = os.path.basename(cfile_path)
    cindex_of_dot = cfile_name.index('.')
    gene_name = cfile_name[:cindex_of_dot]

    alignment = AlignIO.read(file, "fasta")
    aln_taxa_size=len(alignment)
    miss_taxa=taxa_size-aln_taxa_size
    aln_size=alignment.get_alignment_length()

    info_sites=0
    mp1=0
    mp2=0
    mp3=0
    mp4=0

    #def informative Sites and missing sites
    for s in range(aln_size):
        gap=list(alignment[:,s]).count("-")
        miss=list(alignment[:,s]).count("n")
        amb=int(gap)+int(miss)+int(miss_taxa)
        
        amb_ratio=round(float(amb)/float(taxa_size),3)
        if amb_ratio<=0.1:
            mp1=mp1+1
        elif amb_ratio<=0.3 and amb_ratio>0.1:
            mp2=mp2+1
        elif amb_ratio<=0.6 and amb_ratio>0.3:
            mp3=mp3+1
        elif amb_ratio>0.6:
            mp4=mp4+1

        uniq=set(list(alignment[:,s]))
        cond_info_site=0
        if len(uniq)>0:
            if len(uniq)==2:
                if "n" in uniq or "-" in uniq:
                    pass
                else:
                    for n in uniq:
                        count=list(alignment[:,s]).count(n)
                        if count>=2:
                            cond_info_site=cond_info_site+1
                        else:
                            pass
                    if cond_info_site>=2:
                        info_sites=info_sites+1
                    else:
                        pass
            else:
                for n in uniq:
                    if n=="n" or n=="-":
                        pass
                    else:
                        count=list(alignment[:,s]).count(n)
                        if count>=2:
                            cond_info_site=cond_info_site+1
                        else:
                            pass
                if cond_info_site>=2:
                    info_sites=info_sites+1
                else:
                    pass
        else:
            pass


    print ("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (str(gene_name),int(aln_size),int(mp1),int(mp2),int(mp3),int(mp4),int(info_sites)))
