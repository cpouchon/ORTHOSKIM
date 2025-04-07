#! /usr/bin/env python

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

import argparse
import sys
import os, errno
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='Get contigs format for megahit output')
parser.add_argument("-in","--infile", help="input reference DB (fasta format)",
                    type=str)
parser.add_argument("-q","--query", help="taxa name",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)

args = parser.parse_args()

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

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

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


input=args.infile
taxname=args.query
outd=args.outdir

fname=taxname+".fa"
cur_genome=to_dict_remove_dups(SeqIO.parse(input, "fasta"))
n=0
for s in list(cur_genome.keys()):
    n=n+1
    infos=cur_genome[s].description.split(" ")
    seq=cur_genome[s].seq
    for i in infos:
        if i.startswith("multi"):
            cov=i.split("=")[1]
        elif i.startswith("len"):
            clen=i.split("=")[1]
        else:
            pass
    nodename="NODE_"+str(n)+"_length_"+str(clen)+"_cov_"+str(cov)
    out_header=">"+str(nodename)
    if os.path.isfile(os.path.join(outd, fname)):
        with open(os.path.join(outd, fname), 'a+') as file:
            old_headers = []
            end_file=file.tell()
            file.seek(0)
            for line in file:
                if line.startswith(">"):
                    old_headers.append(line.rstrip())
            if not out_header in old_headers:
                file.seek(end_file)
                file.write(out_header+'\n')
                file.write(str(seq)+'\n')
            else:
                pass
    else :
        with open(os.path.join(outd, fname), 'a') as out:
            out.write(out_header+'\n')
            out.write(str(seq)+'\n')

print("--- %s seconds ---" % (time.time() - start_time))
