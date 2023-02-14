#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import sys
import os, errno
import argparse
import time
from Bio.Seq import Seq

start_time = time.time()

parser = argparse.ArgumentParser(description='translation of 353 reference sequences')
parser.add_argument("-i","--infile", help="input fasta file",
                    type=str)
parser.add_argument("-o","--outfile", help="translated fasta file",
                    type=str)
parser.add_argument("--genetic_code", help="NCBI genetic code number used for DNA translation",
                    type=int)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """

    remainder = len(sequence) % 3

    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))

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




orfcond = args.orf
genefile = args.infile
outf=args.outfile
genet_code=args.genetic_code


seqs={}

print("parsing input file")
seqs=to_dict_remove_dups(SeqIO.parse(genefile, "fasta"))
print("done")

Bar = ProgressBar(len(seqs), 60, '\t ')
barp=0
print("extracting sequences")
for s in seqs:
    cond_toprint_var=0
    cond_cs=0
    barp=barp+1
    seqname=s
    seq=str(seqs[s].seq)
    tseq=pad_seq(Seq(seq)).translate(table=genet_code)
    if "*" in tseq:
        newseq_rev=str(Seq(seq).reverse_complement())
        seqframes=[pad_seq(Seq(seq)).translate(table=genet_code),pad_seq(Seq(seq[1:len(seq)])).translate(table=genet_code),pad_seq(Seq(seq[2:len(seq)])).translate(table=genet_code),pad_seq(Seq(newseq_rev)).translate(table=genet_code),pad_seq(Seq(newseq_rev[1:len(newseq_rev)])).translate(table=genet_code),pad_seq(Seq(newseq_rev[2:len(newseq_rev)])).translate(table=genet_code)]
        frame_select=0
        len_orf=0
        orf_select=""
        for p in range(len(seqframes)):
            orf_parts=[]
            start=0
            for i in range(len(str(seqframes[p]))):
                if seqframes[p][i]!="*":
                    end=i
                else:
                    end=i
                    orf_parts.append(str(start)+"-"+str(end))
                    start=i+1
            orf_parts.append(str(start)+"-"+str(end+1))
            for i in range(len(orf_parts)):
                part=orf_parts[i]
                s=int(part.split("-")[0])
                e=int(part.split("-")[1])+1
                if len(seqframes[p][s:e])>len_orf and len(seqframes[p][s:e])>=(float(orfcond)*len(seqframes[p])):
                    len_orf=len(seqframes[p][s:e])
                    orf_select=part
                    frame_select=p
        if len(orf_select)>0:
            newseq=str(seqframes[frame_select]).replace("*","")
            geneid=str(seqname)+"_CS"
            cond_cs=cond_cs+1
        else:
            cond_toprint_var=cond_toprint_var+1
            print(seqname)
    else:
        newseq=tseq
        geneid=str(seqname)
    if cond_toprint_var<1:
        header=">"+str(geneid)
        sequence=newseq
        if os.path.isfile(outf):
            with open(outf, 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not header.replace(">","") in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(sequence)+'\n')
                else:
                    pass
        else :
            with open(outf, 'w') as out:
                out.write(header+'\n')
                out.write(str(sequence)+'\n')
    Bar.update(barp)
print("done")
