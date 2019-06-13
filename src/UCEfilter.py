#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import sys
import os, errno
import argparse





parser = argparse.ArgumentParser(description='')
parser.add_argument("-i","--infile", help="Fasta input alignment from MACSE run of amino acids sequences",
                    type=str)
parser.add_argument("-p","--outdir", help="Out directory path",
                    type=str)
parser.add_argument("-r","--allelefreq", help="minimal allele frequency to make consense (0.5 means majority rule consensus)",
                    type=float)
parser.add_argument("-o","--outname", help="Outname file",
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
                print "path '%s' already exists" % path   # overwrite == False and we've hit a directory that exists
        else: raise


class L(list):
    """
    A subclass of list that can accept additional attributes.
    Should be able to be used just like a regular list.

    The problem:
    a = [1, 2, 4, 8]
    a.x = "Hey!" # AttributeError: 'list' object has no attribute 'x'

    The solution:
    a = L(1, 2, 4, 8)
    a.x = "Hey!"
    print a       # [1, 2, 4, 8]
    print a.x     # "Hey!"
    print len(a)  # 4

    You can also do these:
    a = L( 1, 2, 4, 8 , x="Hey!" )                 # [1, 2, 4, 8]
    a = L( 1, 2, 4, 8 )( x="Hey!" )                # [1, 2, 4, 8]
    a = L( [1, 2, 4, 8] , x="Hey!" )               # [1, 2, 4, 8]
    a = L( {1, 2, 4, 8} , x="Hey!" )               # [1, 2, 4, 8]
    a = L( [2 ** b for b in range(4)] , x="Hey!" ) # [1, 2, 4, 8]
    a = L( (2 ** b for b in range(4)) , x="Hey!" ) # [1, 2, 4, 8]
    a = L( 2 ** b for b in range(4) )( x="Hey!" )  # [1, 2, 4, 8]
    a = L( 2 )                                     # [2]
    """
    def __new__(self, *args, **kwargs):
        return super(L, self).__new__(self, args, kwargs)

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and hasattr(args[0], '__iter__'):
            list.__init__(self, args[0])
        else:
            list.__init__(self, args)
        self.__dict__.update(kwargs)

    def __call__(self, **kwargs):
        self.__dict__.update(kwargs)
        return self

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


def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] >= valueToFind:
            listOfKeys.append(item[0])
    return  listOfKeys


def GetConsensusBaseN(dict_nt,val_factor):
    nt=''
    threshold = val_factor*(sum(dict_nt.values()))
    keys_max = getKeysByValue(dict_nt,threshold)
    if len(keys_max)==1:
        nt = getKeysByValue(dict_nt,threshold)[0]
    elif len(keys_max)>1:
        nt = "X"
    return nt

'#########################################################'
'function'
'#########################################################'

input_file = args.infile
outdir = args.outdir
name = args.outname
frq=args.allelefreq

fname = str(name)+".fna"
'we keep the name of alignment (name of UCE) for header'

file_name = os.path.basename(input_file)
index_of_dot = file_name.index('.')
file_name_without_extension = file_name[:index_of_dot]
seqname = file_name_without_extension.replace("_AA","")

mkdir(outdir)

cur_genome = SeqIO.parse(input_file, "fasta")
seqList = L()
for record in cur_genome:
    sequence=str(record.seq)
    seqList.append(sequence.replace("!","X").replace("*","X"))

n = len(seqList[0])

'we create a profile for each Nt at each position'
profile = { 'A':[0]*n,'B':[0]*n ,'C':[0]*n,'D':[0]*n , 'E':[0]*n, 'F':[0]*n, 'G':[0]*n,'H':[0]*n,'I':[0]*n,'K':[0]*n,'L':[0]*n,'M':[0]*n,'N':[0]*n,'P':[0]*n,'Q':[0]*n,'R':[0]*n,'S':[0]*n,'T':[0]*n,'V':[0]*n,'W':[0]*n,'X':[0]*n,'Y':[0]*n,'Z':[0]*n,'-':[0]*n}

for seq in seqList:
    for i, char in enumerate(seq):
        profile[char][i] += 1

'we keep a consensus position according to the profile'
consensus = ""

Bar = ProgressBar(len(seqList[0]), 60, '\t Extraction of genes for %s' % input_file)
barp=0
for i in range(n):
    Bar.update(barp)
    res=[]
    cons_nt ='x'
    for nt in "ABCDEFGHIKLMNPQRSTVWXYZ-":
        if profile[nt][i] > 0:
            tmp=(nt,profile[nt][i])
            res.append(tmp)
        val=dict(res)
    if len(val)>1:
        val2={}
        for key in val:
            if key != "-":
                val2[key] = val[key]
        if len(val2)>1:
            cons_nt=GetConsensusBaseN(val2,frq)
        else:
            cons_nt=val2.keys()[0]
    else:
        cons_nt=val.keys()[0]
    consensus += cons_nt

'output results'
out_seq=consensus
header=">"+str(seqname)+"_UCE"

if os.path.isfile(os.path.join(outdir, fname)):
    with open(os.path.join(outdir, fname), 'a+') as file:
        old_headers = []
        end_file=file.tell()
        file.seek(0)
        for line in file:
            if line.startswith(">"):
                old_headers.append(line.rstrip())

        if not header in old_headers:
            file.seek(end_file)
            file.write(header+'\n')
            file.write(str(out_seq)+'\n')
        else:
            pass
else :
    with open(os.path.join(outdir, fname), 'w') as out:
        out.write(header+'\n')
        out.write(str(out_seq)+'\n')
