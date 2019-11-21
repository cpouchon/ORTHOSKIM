#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
import random

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation




parser = argparse.ArgumentParser(description='Change name with OrthoSKim format for chloroplastic DB')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-i","--infile", help="input file list of sequence files, or sequence file",
                    type=str)
parser.add_argument("--pathfind", help="[mode] search all sequences files in given path (-p)",
                    action="store_true")
parser.add_argument("--protein", help="[mode] for the CDS regions",
                    action="store_true")
parser.add_argument("-o","--outfile", help="output file list of ref sequences",
                    type=str)
parser.add_argument("-e","--extension", help="extension of files to search (i.e fa,fna,fasta or custom)",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


def cap_gene(phrase):
    return ' '.join(map(lambda s: s[:-1]+s[-1].upper(), phrase.split()))



file_extension=args.extension
if args.pathfind:
    'Search all fasta files files'
    path = args.path
    in_files = []
    for r, d, f in os.walk(path):
        for file in f:
            if str("."+file_extension) in file:
                in_files.append(os.path.join(r, file))
else:
    'We parse a list of files'
    input_list = args.infile
    with open(input_list) as f:
        in_files = f.readlines()


fname=args.outfile

for file in in_files:
    'we parse files'
    seqs={}
    cur_genome = SeqIO.parse(file, "fasta")
    for record in cur_genome:
        seqID=record.id.split(";")[0]
        sequence=record.seq
        seqs.setdefault(seqID, []).append(sequence)

    genename=os.path.basename(file).replace(str("."+file_extension),"")

    for taxa in seqs:
        if args.protein:
            header=">"+str(cap_gene(genename))+"_"+str(taxa.split("@")[0])
            taxa_sequence=str(seqs[taxa][0]).replace("-","X").upper()
        else:
            'condition for RNA'
            infos=taxa.split(" ")[0].split("_")
            header=">"+str(genename)+"_"+infos[1]+"_"+infos[2]
            taxa_sequence=str(seqs[taxa][0]).replace("-","N").upper()

        if os.path.isfile(fname):
            with open(fname, 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not taxa in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(taxa_sequence)+'\n')
                else:
                    pass
        else :
            with open(fname, 'w') as out:
                out.write(header+'\n')
                out.write(str(taxa_sequence)+'\n')
