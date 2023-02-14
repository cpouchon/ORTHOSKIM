#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import *
from Bio.Seq import *
from Bio.SeqUtils import *
from ete3 import NCBITaxa
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
#from Bio.Alphabet import IUPAC, Gapped

#import numpy as np
#import random
import argparse
import sys
import os, errno
import time
import glob

start_time = time.time()

parser = argparse.ArgumentParser(description='Formate UCE mega353 references to Orthoskim DB. Script was writen by C. Pouchon (2020).')
parser.add_argument("-i","--infile", help="Fasta input contigs file",
                    type=str)
parser.add_argument("-n","--names", help="seq names codes",
                    type=str)
parser.add_argument("-o","--outfile", help="Out directory path",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

try:
    import ast
    import inspect
    import sys
    print("Patching NCBITaxa's base methods.\n")
    code_to_patch = """db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))"""
    patched_code = """db.execute("INSERT OR REPLACE INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))"""

    ncbiquery = sys.modules[NCBITaxa.__module__]
    lines_code = [x.replace(code_to_patch, patched_code)
                  for x in inspect.getsourcelines(ncbiquery.upload_data)[0]]
    # Insert info message to see if patch is really applied
    lines_code.insert(1, "    print('\\nIf this message shown, then the patch is successful!')\n")
    # Insert external import and constants since only this function is patched and recompiled
    lines_code.insert(1, "    import os, sqlite3, sys\n")
    lines_code.insert(1, "    DB_VERSION = 2\n")
    lines_code = "".join(lines_code)

    # Compile and apply the patch
    ast_tree = ast.parse(lines_code)
    patched_function = compile(ast_tree, "<string>", mode="exec")
    mod_dummy = {}
    exec(patched_function, mod_dummy)
    ncbiquery.upload_data = mod_dummy["upload_data"]
except Exception:
    print("Patching failed, current taxonomy data downloaded from FTP may be failed to update with ETE3!")
finally:
    print("Patch finished.")

ncbi = NCBITaxa()

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
		out = '\r %20s |%s%s| %3d %% ' % (self.title, '.' * bar, ' ' * (self.maxbar - bar), perc)
		sys.stdout.write(out)
		sys.stdout.flush()


inf=args.infile
outf=args.outfile
names=args.names

with open(names) as f:
    tab = f.readlines()
dictname=dict()
for g in tab:
    g_tab = g.rstrip().split('\t')
    dictname[g_tab[0]]=g_tab[1]

# importer code IUPAC
ht_basis=["R","Y","S","W","K","M","B","D","H","V","N"]


# et virer sequences avec IUPAC

print("parsing input file")
seqs=to_dict_remove_dups(SeqIO.parse(inf, "fasta"))
print("done")

Bar = ProgressBar(len(seqs), 60, '\t ')
barp=0
badseq=0
print("extracting sequences")
for s in seqs:
    barp=barp+1
    tmp=s.split(" ")[0]
    geneid=tmp.split("-")[len(tmp.split("-"))-1]
    seqname=tmp.split("-")[0].split("_")[0]
    if seqname in list(dictname.keys()):
        spname=dictname[seqname]
        if len(ncbi.get_name_translator([spname]))>0:
            taxid=int(ncbi.get_name_translator([spname])[spname][0])
        else:
            gename=spname.split(" ")[0]
            if len(ncbi.get_name_translator([gename]))>0:
                taxid=int(ncbi.get_name_translator([gename])[gename][0])
            else:
                continue
        header=">"+str(geneid)+"_"+str(taxid)+"_"+str(spname).replace(" ","_")
        sequence=str(seqs[s].seq).upper().replace("N","").replace("-","")
        hb=[x for x in ht_basis if x in sequence]
        if len(hb)>0:
            badseq+=1
            continue
        else:
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
print("total sequences: ",len(seqs))
print("sequences removed: ",str(badseq))
print("done")
