#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


parser = argparse.ArgumentParser(description='Concatenation of sequenes alignments')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-i","--infile", help="input file list of sequence files, sequences have to be aligned before",
                    type=str)
parser.add_argument("-pfind","--pathfind", help="[mode] search all sequences files in given path (-p) otherwise parse a given list (-i)",
                    action="store_true")
parser.add_argument("-o","--outfile", help="output file list of concatenated sequences",
                    type=str)
parser.add_argument("-t","--taxa", help="list of taxa to include in the final concatenated file",
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


if args.pathfind:
    'Search all fasta files files'
    path = args.path
    in_files = []
    for r, d, f in os.walk(path):
        for file in f:
            if '.fa' in file:
                in_files.append(os.path.join(r, file))
else:
    'We parse a list of files'
    input_list = args.infile
    with open(input_list) as f:
        in_files = f.readlines()

#cat *.fa | grep '>' | sort | uniq | perl -pe 's/>//' > list_taxa

fname=args.outfile

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    #taxa_list.append(l)
    taxa_dict.setdefault(l, [])

end=0
'read sequences and stock the id/seqs in dictionnary'
for file in in_files:
    seqs={}
    cur_genome = SeqIO.parse(file, "fasta")
    for record in cur_genome:
        seqID=record.id
        sequence=record.seq
        seqs.setdefault(seqID, []).append(sequence)

    start=end+1
    end=(start-1)+len(str(seqs[seqs.keys()[0]][0]))
    genename=os.path.basename(file)
    'a partition file is writed to help for phylogenetic inference'
    print "%s\t%s\t%s" % (start,end,genename)

    for taxa in seqs.keys():
        if taxa in taxa_dict.keys():
            taxa_dict[taxa].append(str(seqs[taxa][0]))
        else:
            'we create a null sequence for missing taxa per gene'
            lenseq = len(str(seqs[seqs.keys()[0]][0]))
            nullseq = "N"*lenseq
            taxa_dict[taxa].append(nullseq)

'concatenation of sequences in dictionary and output sequences'
for taxa in taxa_dict:
    header= ">"+str(taxa)
    concatenate_seq="".join(taxa_dict[taxa])

    if os.path.isfile(fname):
        with open(fname, 'a+') as file:
            old_headers = []
            end_file=file.tell()
            file.seek(0)
            for line in file:
                if line.startswith(">"):
                    old_headers.append(line.rstrip())

            if not header in old_headers:
                file.seek(end_file)
                file.write(header+'\n')
                file.write(str(concatenate_seq)+'\n')
            else:
                pass
    else :
        with open(fname, 'w') as out:
            out.write(header+'\n')
            out.write(str(concatenate_seq)+'\n')
