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
from collections import Counter
#from Bio.Alphabet import IUPAC, Gapped

#import numpy as np
#import random
import argparse
import sys
import os, errno
import time
import glob
import random

start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering of errors from alignments. Script was writen by C. Pouchon (2020).')
parser.add_argument("-i","--inpath", help="in path of alignments fasta files",
                    type=str)
parser.add_argument("-e","--extension", help="file extension (eg. trim/fa/fasta/fna/etc)",
                    type=str)
parser.add_argument("-q","--query", help="queried taxonomic level for consensus",
                    type=str, choices=["family","genus"])
parser.add_argument("-o","--outpath", help="outpath for selected fasta",
                    type=str)
parser.add_argument("-gw","--gwindow", help="sliding window size to identify polymorphic sites within the genus (eg. 100bp)",
                    type=int)
parser.add_argument("-gp","--fpsites", help="maxmimal polymorphic sites allowed within sliding window within the genus (eg. 20bp)",
                    type=int)
parser.add_argument("-fw","--fwindow", help="sliding window size to identify polymorphic sites within the genus (eg. 100bp)",
                    type=int)
parser.add_argument("-fp","--fpsites", help="maxmimal polymorphic sites allowed within sliding window within the genus (eg. 20bp)",
                    type=int)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

print("Errors checking and filtering")

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

def mkdir(path, overwrite=False):
    '''
    function to create a directory for output fasta
    '''
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not overwrite:
                pass#print ("path '%s' already exists" % path)   # overwrite == False and we've hit a directory that exists
        else: raise

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


try:
    # Python 2
    xrange
except NameError:
    # Python 3, xrange is now named range
    xrange = range

def get_maj_sites(string):
    from collections import Counter
    counts=Counter(string)
    keys=[]
    #tot=sum(counts.values())
    for key, value in counts.items():
        if value == max(counts.values()):
            keys.append(key)
    return keys

def majority_rule_consense(msa_summary):
    import random
    consensus = ''
    for i in xrange(msa_summary.alignment.get_alignment_length()):
        possibles = get_maj_sites(msa_summary.get_column(i))
        if len(possibles) == 1:
            consensus += possibles[0]
        else:
            consensus += random.choice(possibles)
    return consensus

# load all fasta
inpath=args.inpath
w_size=args.window
w_site=args.psites
outpath=args.outpath
query_level=args.query


in_files = []
path_to_seek=os.walk(os.path.join(inpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(str("."+str(args.extension))):
            in_files.append(os.path.join(r, file))

stored=dict()
alignments=dict()
stored_for_cons=dict()
