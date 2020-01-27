#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation




parser = argparse.ArgumentParser(description='Extraction of given positions start:end of alignment')
parser.add_argument("--infile", help="input file list of sequence files, sequences have to be aligned before",
                    type=str)
parser.add_argument("--outfile", help="output file",
                    type=str)
parser.add_argument("--start", help="start of alignment",
                    type=int)
parser.add_argument("--end", help="end of alignment",
                    type=int)
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


start=args.start
end=args.end
fname=args.outfile

file=args.infile
seqs={}
cur_genome = SeqIO.parse(file, "fasta")
for record in cur_genome:
    seqID=record.id.split(";")[0]
    sequence=record.seq
    seqs.setdefault(seqID, []).append(sequence)

for taxa in seqs:
    s = start-1
    e = end
    header = ">"+str(taxa)
    extractseq=str(seqs[taxa][0][int(s):int(e)])

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
                file.write(str(extractseq)+'\n')
            else:
                pass
    else :
        with open(fname, 'w') as out:
            out.write(header+'\n')
            out.write(str(extractseq)+'\n')
