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

parser = argparse.ArgumentParser(description='Filtering of paralogs from alignments. Script was writen by C. Pouchon (2020).')
parser.add_argument("-i","--inpath", help="in path of alignments fasta files",
                    type=str)
parser.add_argument("-c","--consensus", help="consensus rule (eg. 0.5 for majority rule consensus)",
                    type=float)
parser.add_argument("-e","--extension", help="file extension (eg. trim/fa/fasta/fna/etc)",
                    type=str)
parser.add_argument("-q","--query", help="queried taxonomic level for consensus",
                    type=str, choices=["order","family","genus"])
parser.add_argument("-o","--outpath", help="outpath for selected fasta",
                    type=str)
parser.add_argument("-w","--window", help="sliding window size to identify polymorphic sites (eg. 100bp)",
                    type=int)
parser.add_argument("-p","--psites", help="maxmimal polymorphic sites allowed within sliding window (eg. 20bp)",
                    type=int)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

print("Paralog checking and filtering")

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




# load all fasta
inpath=args.inpath
#inpath="Alignment/mitochondrion_CDS/trim/"

rule_cons=args.consensus
w_size=args.window
w_site=args.psites
outpath=args.outpath
#outpath="/Users/pouchonc/Lecythidaceae/alignment/filter_paralogs"
query_level=args.query


in_files = []
path_to_seek=os.walk(os.path.join(inpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(str("."+str(args.extension))):
            in_files.append(os.path.join(r, file))

stored=dict()
alignments=dict()


# store all sequences
print("parsing of input files")
Bar = ProgressBar(len(in_files), 60, '\t ')
barp=0
for file in in_files:
    #barp=barp+1
    geneID=os.path.basename(file).replace(".fasta","").replace(str("."+"fa"),"").replace(".trim","")
    seqs=to_dict_remove_dups(SeqIO.parse(file, "fasta"))
    if geneID in list(stored.keys()):
        pass
    else:
        stored[geneID]=dict()
        alignments[geneID]=dict()
        stored[geneID]["other"]=list()
        #alignments[geneID]["other"]=MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
        #alignments[geneID]["other"]=MultipleSeqAlignment([])
    for s in seqs:
        genus=s.split("_")[0]
        if len(ncbi.get_name_translator([genus]))>0:
            taxid=int(ncbi.get_name_translator([genus])[genus][0])
            sub_lineages=ncbi.get_lineage(taxid)
            sub_names=ncbi.get_taxid_translator(sub_lineages)
            sub_ranks=ncbi.get_rank(sub_names.keys())
            list_fam=[key  for (key, value) in sub_ranks.items() if value == query_level.lower()]
            if len(list_fam)==0:
                stored[geneID]["other"].append(seqs[s])
                #alignments[geneID]["other"].add_sequence(s, str(seqs[s].seq))
            else:
                fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                fam_taxid=int(list_fam[0])
                if fam in list(stored[geneID].keys()):
                    stored[geneID][fam].append(seqs[s])
                    #alignments[geneID][fam].add_sequence(s, str(seqs[s].seq))
                else:
                    stored[geneID][fam]=list()
                    stored[geneID][fam].append(seqs[s])
                    #alignments[geneID][fam] = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
                    #alignments[geneID][fam] = MultipleSeqAlignment([])
                    #alignments[geneID][fam].add_sequence(s, str(seqs[s].seq))
    for f in list(stored[geneID].keys()):
        alignments[geneID][f]=MultipleSeqAlignment(stored[geneID][f])
    Bar.update(barp)
print("")
print("done")
print("")

# get consensus seq per family
print("Get consensus sequences")
consense=dict()
Bar = ProgressBar(len(stored), 60, '\t ')
barp=0
for g in stored:
    consense[g]=dict()
    barp=barp+1
    for f in stored[g]:
        if len(stored[g][f])>0:
            align = alignments[g][f]
            summary_align = AlignInfo.SummaryInfo(align)
            cons=str(summary_align.dumb_consensus(rule_cons))
            consense[g][f]=cons
        else:
            continue
    Bar.update(barp)
print("")
print("done")
# stored bad samples

Bar = ProgressBar(len(stored), 60, '\t ')
barp=0
print("checking for paralogs")
bad_seq=dict()
for g in stored:
    barp=barp+1
    bad_seq[g]=list()
    for f in stored[g]:
        if len(stored[g][f])>0:
            #align = alignments[g][f]
            #summary_align = AlignInfo.SummaryInfo(align)
            #cons=str(summary_align.dumb_consensus(rule_cons))
            cons=consense[g][f]
            lenseq=len(cons)
            # cond to check polymorphic sites
            for i in range(lenseq - w_size + 1):
                bmin=i
                bmax=i+w_size
                poly_dict=dict.fromkeys(list(stored[g][f].keys()), 0)
                for j in range(bmin,bmax):
                    #tmpsite=list()
                    for s in stored[g][f]:
                        #tmpsite.append(stored[g][f][s][j])
                        # polymorphic sites - paralogs
                        if stored[g][f][s][j]=="-" or stored[g][f][s][j]=="N" or stored[g][f][s][j]=="n":
                            continue
                        else:
                            if stored[g][f][s][j]!=cons[j]:
                                poly_dict[s]+=1
                            else:
                                continue
                tmp_samp=[k for (k,v) in poly_dict.items() if v >= w_site]
                if len(tmp_samp)>0:
                    for s in tmp_samp:
                        if s in bad_seq[g]:
                            pass
                        else:
                            bad_seq[g].append(s)
            Bar.update(barp)
print("")
print("done")

# start_time = time.time()
# Bar = ProgressBar(len(stored), 60, '\t ')
# barp=0
# print("checking for paralogs")
# bad_seq=dict()
# poly_dict=dict()
# for g in stored:
#     barp=barp+1
#     bad_seq[g]=list()
#     lenseq=len(consense[g][list(consense[g].keys())[0]])
#     poly_dict[g]=dict()
#     for i in range(lenseq - w_size + 1):
#         bmin=i
#         bmax=i+w_size
#         for f in stored[g]:
#             if len(stored[g][f])>0:
#                 poly_dict[g][f]=dict.fromkeys(list(stored[g][f].keys()), 0)
#                 for j in range(bmin,bmax):
#                     for s in stored[g][f]:
#                         if stored[g][f][s][j]=="-":
#                             continue
#                         else:
#                             if stored[g][f][s][j]!=consense[g][f][j]:
#                                 poly_dict[g][f][s]+=1
#                             else:
#                                 continue
#                 tmp_samp=[k for (k,v) in poly_dict[g][f].items() if v >= w_site]
#                 if len(tmp_samp)>0:
#                     for s in tmp_samp:
#                         if s in bad_seq[g]:
#                             pass
#                         else:
#                             bad_seq[g].append(s)
#             else:
#                 continue
#     Bar.update(barp)
# print("")
# print("done")
# print("--- %s seconds ---" % (time.time() - start_time))

# write fasta out without bad SAMPLES
mkdir(str(outpath))

for g in stored:
    fname=str(g)+".fa"
    for f in stored[g]:
        if len(stored[g][f])>0:
            for s in stored[g][f]:
                if s in bad_seq[g]:
                    pass
                else:
                    header=">"+str(s)
                    newseq2=str(stored[g][f][s].seq)
                    if os.path.isfile(os.path.join(outpath, fname)):
                        with open(os.path.join(outpath, fname), 'a+') as file:
                            old_headers = []
                            end_file=file.tell()
                            file.seek(0)
                            for line in file:
                                if line.startswith(">"):
                                    old_headers.append(line.rstrip().replace(">","").split(";")[0])
                            if not str(s) in old_headers:
                                file.seek(end_file)
                                file.write(header+'\n')
                                file.write(str(newseq2)+'\n')
                            else:
                                pass
                    else:
                        with open(os.path.join(outpath, fname), 'w') as out:
                            out.write(header+'\n')
                            out.write(str(newseq2)+'\n')



# print bad
cond_print=0
for g in bad_seq:
    if len(bad_seq[g])>0:
        cond_print=cond_print+1

if cond_print>0:
    print("Paralogs found:")
    for g in bad_seq:
        if len(bad_seq[g])>0:
            print(g,len(bad_seq[g]),",".join(bad_seq[g]))
print("")
print("Time ellapsed:")
print("--- %s seconds ---" % (time.time() - start_time))
