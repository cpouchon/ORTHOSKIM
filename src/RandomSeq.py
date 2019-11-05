#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
import random

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation




parser = argparse.ArgumentParser(description='Extraction of random sequences')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-i","--infile", help="input file list of sequence files, or sequence file",
                    type=str)
parser.add_argument("--pathfind", help="[mode] search all sequences files in given path (-p)",
                    action="store_true")
parser.add_argument("--multiple", help="[mode] extraction from multiple sequences files in given path (-p) otherwise parse a given list (-i)",
                    action="store_true")
parser.add_argument("--single", help="[mode] extraction from a single given fasta file (-i)",
                    action="store_true")
parser.add_argument("-o","--outfile", help="output file list of concatenated sequences",
                    type=str)
parser.add_argument("-c","--count", help="number of sequences targeted",
                    type=int)
parser.add_argument("-e","--extension", help="extension of files to search (i.e fa,fna,fasta or custom)",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


class ProgressBar:
	'''
	Progress bar
	'''
	def __init__ (self, valmax, maxbar, title):
		if valmax == 0:  valmax = 1
		if maxbar > 200: maxbar = 200
		self.valmax = valmax
		self.maxbar = maxbar
		self.title  = title
	def update(self, val):
		import sys
		perc  = round((float(val) / float(self.valmax)) * 100)
		scale = 100.0 / float(self.maxbar)
		bar   = int(perc / scale)
		out = '\r %20s [%s%s] %3d %% ' % (self.title, '.' * bar, ' ' * (self.maxbar - bar), perc)
		sys.stdout.write(out)
		sys.stdout.flush()


file_extension=args.extension
number_of_seqs=args.count

if args.multiple:
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

if args.single:
    in_files = args.infile

#cat *.fa | grep '>' | sort | uniq | perl -pe 's/>//' > list_taxa

fname=args.outfile

'read sequences and stock the id/seqs in dictionnary'
for file in in_files:
    seqs={}
    dicinfo={}
    cur_genome = SeqIO.parse(file, "fasta")
    for record in cur_genome:
        seqID=record.id.split(";")[0]
        sequence=record.seq
        seqs.setdefault(seqID, []).append(sequence)
        dicinfo.setdefault(seqID, []).append(record.description.replace(";",""))

    genename=os.path.basename(file).replace(str("."+file_extension),"")

    if len(seqs) >= number_of_seqs:
        random_keys=random.sample(list(seqs),k=number_of_seqs)
    else:
        print ("WARN: the number of targeted sequences is greater than the number of sequences found into files. By default, all sequences are kept.")
        random_keys=random.sample(list(seqs),k=len(seqs))

    for taxa in random_keys:
        taxa_sequence=str(seqs[taxa][0]).replace("-","n")
        infos = dicinfo[taxa][0].split(" ")[1:len(dicinfo[taxa][0].split(" "))]
        taxinfo = dict(item.split("=") for item in infos)
        header=">"+str(taxinfo['gene'])+"_"+str(taxinfo['taxid'])

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
