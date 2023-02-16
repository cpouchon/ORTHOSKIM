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
parser.add_argument("-o","--outpath", help="outpath for selected fasta",
                    type=str)
parser.add_argument("-gw","--gwindow", help="sliding window size to identify polymorphic sites within the genus (eg. 100bp)",
                    type=int)
parser.add_argument("-gp","--gpsites", help="maxmimal polymorphic sites allowed within sliding window within the genus (eg. 20bp)",
                    type=int)
parser.add_argument("-fw","--fwindow", help="sliding window size to identify polymorphic sites within the genus (eg. 100bp)",
                    type=int)
parser.add_argument("-fp","--fpsites", help="maxmimal polymorphic sites allowed within sliding window within the genus (eg. 20bp)",
                    type=int)
parser.add_argument("-gt","--gtaxa", help="minimal taxa required for the consensus per genus level",
                    type=int)
parser.add_argument("-ft","--ftaxa", help="minimal taxa required for the consensus per family level",
                    type=int)
parser.add_argument("-t","--taxonomy", help="expected taxonomy rank for the taxid selection (eg. Viridiplantae)",
                    type=str)
parser.add_argument("-co","--cons_outpath", help="outpath for consensus sequences",
                    type=str)
parser.add_argument("--export", help="[mode] export consensus sequences",
                    action="store_true")
# mettre condition pour write outpath
# et ajouter le cons_path pour out les consensus

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

print("Errors checking and filtering")
print("")

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
gw_size=args.gwindow
gw_site=args.gpsites
fw_size=args.fwindow
fw_site=args.fpsites
outpath=args.outpath
expected_taxo=args.taxonomy
gmin=args.gtaxa
fmin=args.ftaxa

in_files = []
path_to_seek=os.walk(os.path.join(inpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(str("."+str(args.extension))):
            in_files.append(os.path.join(r, file))

stored=dict()
alignments=dict()
stored_for_cons=dict()
stored_for_cons_fam=dict()
alignments=dict()
alignments_fam=dict()
link_gen_fam=dict()

# store all sequences
print("parsing of input files")
Bar = ProgressBar(len(in_files), 60, '\t ')
barp=0
for file in in_files:
    barp=barp+1
    geneID=os.path.basename(file).replace(".fasta","").replace(str("."+"fa"),"").replace(".trim","")
    seqs=to_dict_remove_dups(SeqIO.parse(file, "fasta"))
    if geneID in list(stored.keys()):
        pass
    else:
        stored[geneID]=dict()
        stored_for_cons[geneID]=dict()
        stored_for_cons_fam[geneID]=dict()
        alignments[geneID]=dict()
        alignments_fam[geneID]=dict()
        stored[geneID]["other"]=dict()
        stored_for_cons[geneID]["other"]=list()
        stored_for_cons_fam[geneID]["other"]=list()
    for s in seqs:
        genus=s.replace("_R_","").split("_")[0]
        taxids=[]
        for i in s.split("_"):
            if i.isdigit():
                taxids.append(int(i))
        if len(taxids)>0:
            if len(ncbi.get_taxid_translator([taxids[0]]))>0:
                tax=taxids[0]
                sub_lineages=ncbi.get_lineage(int(tax))
                sub_names=ncbi.get_taxid_translator(sub_lineages)
                sub_ranks=ncbi.get_rank(sub_names.keys())
                list_genus=[key  for (key, value) in sub_ranks.items() if value == "genus"]
                list_fam=[key  for (key, value) in sub_ranks.items() if value == "family"]
                if len(list_fam)==0:
                    fam="other"
                    stored_for_cons_fam[geneID]["other"].append(seqs[s])
                else:
                    fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                    if fam in list(stored_for_cons_fam[geneID].keys()):
                        stored_for_cons_fam[geneID][fam].append(seqs[s])
                    else:
                        stored_for_cons_fam[geneID][fam]=list()
                        stored_for_cons_fam[geneID][fam].append(seqs[s])
                if len(list_genus)==0:
                    gen="other"
                    stored[geneID]["other"][s]=seqs[s]
                    stored_for_cons[geneID]["other"].append(seqs[s])
                else:
                    gen=ncbi.get_taxid_translator(list_genus)[list_genus[0]]
                    if gen in list(stored[geneID].keys()):
                        if s in list(stored[geneID][gen].keys()):
                            pass
                        else:
                            stored[geneID][gen][s]=seqs[s]
                            stored_for_cons[geneID][gen].append(seqs[s])
                    else:
                        stored[geneID][gen]=dict()
                        stored[geneID][gen][s]=seqs[s]
                        stored_for_cons[geneID][gen]=list()
                        stored_for_cons[geneID][gen].append(seqs[s])
                link_gen_fam[gen]=fam
            else:
                if len(ncbi.get_name_translator([genus]))>0:
                    taxids=ncbi.get_name_translator([genus])[genus]
                    for tax in taxids:
                        if str(expected_taxo) not in list(ncbi.get_taxid_translator(ncbi.get_lineage(tax)).values()):
                            continue
                        else:
                            sub_lineages=ncbi.get_lineage(int(tax))
                            sub_names=ncbi.get_taxid_translator(sub_lineages)
                            sub_ranks=ncbi.get_rank(sub_names.keys())
                            list_genus=[key  for (key, value) in sub_ranks.items() if value == "genus"]
                            list_fam=[key  for (key, value) in sub_ranks.items() if value == "family"]
                            if len(list_fam)==0:
                                fam="other"
                                stored_for_cons_fam[geneID]["other"].append(seqs[s])
                            else:
                                fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                                if fam in list(stored_for_cons_fam[geneID].keys()):
                                    stored_for_cons_fam[geneID][fam].append(seqs[s])
                                else:
                                    stored_for_cons_fam[geneID][fam]=list()
                                    stored_for_cons_fam[geneID][fam].append(seqs[s])
                            if len(list_genus)==0:
                                gen="other"
                                stored[geneID]["other"][s]=seqs[s]
                                stored_for_cons[geneID]["other"].append(seqs[s])
                            else:
                                gen=ncbi.get_taxid_translator(list_genus)[list_genus[0]]
                                if gen in list(stored[geneID].keys()):
                                    if s in list(stored[geneID][gen].keys()):
                                        pass
                                    else:
                                        stored[geneID][gen][s]=seqs[s]
                                        stored_for_cons[geneID][gen].append(seqs[s])
                                else:
                                    stored[geneID][gen]=dict()
                                    stored[geneID][gen][s]=seqs[s]
                                    stored_for_cons[geneID][gen]=list()
                                    stored_for_cons[geneID][gen].append(seqs[s])
                            link_gen_fam[gen]=fam
        else:
            if len(ncbi.get_name_translator([genus]))>0:
                taxids=ncbi.get_name_translator([genus])[genus]
                for tax in taxids:
                    if str(expected_taxo) not in list(ncbi.get_taxid_translator(ncbi.get_lineage(tax)).values()):
                        continue
                    else:
                        sub_lineages=ncbi.get_lineage(int(tax))
                        sub_names=ncbi.get_taxid_translator(sub_lineages)
                        sub_ranks=ncbi.get_rank(sub_names.keys())
                        list_genus=[key  for (key, value) in sub_ranks.items() if value == "genus"]
                        list_fam=[key  for (key, value) in sub_ranks.items() if value == "family"]
                        if len(list_fam)==0:
                            fam="other"
                            stored_for_cons_fam[geneID]["other"].append(seqs[s])
                        else:
                            fam=ncbi.get_taxid_translator(list_fam)[list_fam[0]]
                            if fam in list(stored_for_cons_fam[geneID].keys()):
                                stored_for_cons_fam[geneID][fam].append(seqs[s])
                            else:
                                stored_for_cons_fam[geneID][fam]=list()
                                stored_for_cons_fam[geneID][fam].append(seqs[s])
                        if len(list_genus)==0:
                            gen="other"
                            stored[geneID]["other"][s]=seqs[s]
                            stored_for_cons[geneID]["other"].append(seqs[s])
                        else:
                            gen=ncbi.get_taxid_translator(list_genus)[list_genus[0]]
                            if gen in list(stored[geneID].keys()):
                                if s in list(stored[geneID][gen].keys()):
                                    pass
                                else:
                                    stored[geneID][gen][s]=seqs[s]
                                    stored_for_cons[geneID][gen].append(seqs[s])
                            else:
                                stored[geneID][gen]=dict()
                                stored[geneID][gen][s]=seqs[s]
                                stored_for_cons[geneID][gen]=list()
                                stored_for_cons[geneID][gen].append(seqs[s])
                        link_gen_fam[gen]=fam
    for g in list(stored[geneID].keys()):
        alignments[geneID][g]=MultipleSeqAlignment(stored_for_cons[geneID][g])
    for f in list(stored_for_cons_fam[geneID].keys()):
        alignments_fam[geneID][f]=MultipleSeqAlignment(stored_for_cons_fam[geneID][f])
    Bar.update(barp)
print("")
print("done")
print("")

# get consensus seq per family
print("Get consensus sequences - genus level")
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
            cons=majority_rule_consense(summary_align)
            consense[g][f]=cons
        else:
            continue
    Bar.update(barp)
print("")
print("done")

# get consensus seq per family
print("Get consensus sequences - family level")
consense_fam=dict()
Bar = ProgressBar(len(stored), 60, '\t ')
barp=0
for g in stored_for_cons_fam:
    consense_fam[g]=dict()
    barp=barp+1
    for f in stored_for_cons_fam[g]:
        if len(stored_for_cons_fam[g][f])>0:
            align = alignments_fam[g][f]
            summary_align = AlignInfo.SummaryInfo(align)
            cons=majority_rule_consense(summary_align)
            consense_fam[g][f]=cons
        else:
            continue
    Bar.update(barp)
print("")
print("done")


print("checking for errors/paralogs")
bad_seq=dict()
Bar = ProgressBar(len(stored), 60, '\t ')
barp=0
for g in stored:
    barp=barp+1
    bad_seq[g]=list()
    for f in stored[g]:
        if len(stored[g][f])>0:
            if len(stored[g][f])<gmin:
                #less than gmin taxa so consensus comp at the family level
                famil=link_gen_fam[f]
                if famil in list(consense_fam[g].keys()):
                    # voir condition fam
                    cons=consense_fam[g][famil]
                    lenseq=len(cons)
                    if lenseq==0:
                        #Bar.update(barp)
                        continue
                    else:
                        # cond for minimal number of taxa at the family level
                        if len(stored_for_cons_fam[g][famil])<fmin:
                            continue
                        else:
                            # cond to check polymorphic sites
                            if lenseq<fw_size:
                                # condition to use the ratio if the gene size is shorter than the window size
                                fratio=float(fw_site)/float(fw_size)
                                for i in range(lenseq+1):
                                    bmin=i
                                    bmax=lenseq+1
                                    poly_dict=dict.fromkeys(list(stored[g][f].keys()), 0)
                                    for j in range(bmin,bmax):
                                        for s in stored[g][f]:
                                            # polymorphic sites - paralogs
                                            if stored[g][f][s][j]=="-" or stored[g][f][s][j]=="N" or stored[g][f][s][j]=="n":
                                                continue
                                            else:
                                                if cons[j]=="-" or cons[j]=="N" or cons[j]=="n":
                                                    continue
                                                else:
                                                    if stored[g][f][s][j]!=cons[j]:
                                                        poly_dict[s]+=1
                                                    else:
                                                        continue
                                    tmp_samp=[k for (k,v) in poly_dict.items() if float(v)/float(lenseq) >= fratio]
                                    if len(tmp_samp)>0:
                                        for s in tmp_samp:
                                            if s in bad_seq[g]:
                                                pass
                                            else:
                                                bad_seq[g].append(s)
                            else:
                                for i in range(lenseq - fw_size + 1):
                                    bmin=i
                                    bmax=i+fw_size
                                    poly_dict=dict.fromkeys(list(stored[g][f].keys()), 0)
                                    for j in range(bmin,bmax):
                                        for s in stored[g][f]:
                                            # polymorphic sites - paralogs
                                            if stored[g][f][s][j]=="-" or stored[g][f][s][j]=="N" or stored[g][f][s][j]=="n":
                                                continue
                                            else:
                                                if cons[j]=="-" or cons[j]=="N" or cons[j]=="n":
                                                    continue
                                                else:
                                                    if stored[g][f][s][j]!=cons[j]:
                                                        poly_dict[s]+=1
                                                    else:
                                                        continue
                                    tmp_samp=[k for (k,v) in poly_dict.items() if v >= fw_site]
                                    if len(tmp_samp)>0:
                                        for s in tmp_samp:
                                            if s in bad_seq[g]:
                                                pass
                                            else:
                                                bad_seq[g].append(s)
            else:
                cons=consense[g][f]
                lenseq=len(cons)
                if lenseq==0:
                    #Bar.update(barp)
                    continue
                else:
                    if lenseq<gw_size:
                        # condition to use the ratio if the gene size is shorter than the window size
                        gratio=float(gw_site)/float(gw_size)
                        for i in range(lenseq+1):
                            bmin=i
                            bmax=lenseq+1
                            poly_dict=dict.fromkeys(list(stored[g][f].keys()), 0)
                            for j in range(bmin,bmax):
                                for s in stored[g][f]:
                                    # polymorphic sites - paralogs
                                    if stored[g][f][s][j]=="-" or stored[g][f][s][j]=="N" or stored[g][f][s][j]=="n":
                                        continue
                                    else:
                                        if cons[j]=="-" or cons[j]=="N" or cons[j]=="n":
                                            continue
                                        else:
                                            if stored[g][f][s][j]!=cons[j]:
                                                poly_dict[s]+=1
                                            else:
                                                continue
                            tmp_samp=[k for (k,v) in poly_dict.items() if float(v)/float(lenseq) >= gratio]
                            if len(tmp_samp)>0:
                                for s in tmp_samp:
                                    if s in bad_seq[g]:
                                        pass
                                    else:
                                        bad_seq[g].append(s)
                    else:
                        # cond to check polymorphic sites
                        for i in range(lenseq - gw_size + 1):
                            bmin=i
                            bmax=i+gw_size
                            poly_dict=dict.fromkeys(list(stored[g][f].keys()), 0)
                            for j in range(bmin,bmax):
                                for s in stored[g][f]:
                                    # polymorphic sites - paralogs
                                    if stored[g][f][s][j]=="-" or stored[g][f][s][j]=="N" or stored[g][f][s][j]=="n":
                                        continue
                                    else:
                                        if cons[j]=="-" or cons[j]=="N" or cons[j]=="n":
                                            continue
                                        else:
                                            if stored[g][f][s][j]!=cons[j]:
                                                poly_dict[s]+=1
                                            else:
                                                continue
                            tmp_samp=[k for (k,v) in poly_dict.items() if v >= gw_site]
                            if len(tmp_samp)>0:
                                for s in tmp_samp:
                                    if s in bad_seq[g]:
                                        pass
                                    else:
                                        bad_seq[g].append(s)
    Bar.update(barp)
print("")
print("done")

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

# write consensus sequence in an output dir by genus and family ()

if args.export:
    print("writing of consensus sequences")
    cons_out=args.cons_outpath
    mkdir(str(cons_out))
    mkdir(str(cons_out)+"/family/")
    mkdir(str(cons_out)+"/genus/")
    print("family level sequences")
    Bar = ProgressBar(len(consense_fam), 60, '\t ')
    barp=0
    for g in consense_fam:
        outpath=str(cons_out)+"/family/"
        barp=barp+1
        for f in consense_fam[g]:
            outname=str(f)+".fa"
            header=">"+str(g)+"_"+str(f)
            newseq2=str(consense_fam[g][f]).replace("-","").replace("N","").replace("n","")
            if os.path.isfile(os.path.join(outpath, outname)):
                with open(os.path.join(outpath, outname), 'a+') as file:
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
                with open(os.path.join(outpath, outname), 'w') as out:
                    out.write(header+'\n')
                    out.write(str(newseq2)+'\n')
            # il faut printer par fam.fa et genus.fa
        Bar.update(barp)
    print("genus level sequences")
    Bar = ProgressBar(len(consense), 60, '\t ')
    barp=0
    for g in consense:
        outpath=str(cons_out)+"/genus/"
        barp=barp+1
        for f in consense[g]:
            outname=str(f)+".fa"
            header=">"+str(g)+"_"+str(f)
            newseq2=str(consense[g][f]).replace("-","").replace("N","").replace("n","")
            if os.path.isfile(os.path.join(outpath, outname)):
                with open(os.path.join(outpath, outname), 'a+') as file:
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
                with open(os.path.join(outpath, outname), 'w') as out:
                    out.write(header+'\n')
                    out.write(str(newseq2)+'\n')
            # il faut printer par fam.fa et genus.fa
        Bar.update(barp)
    print("")
    print("done")

print("")
print("Time ellapsed:")
print("--- %s seconds ---" % (time.time() - start_time))


## Gene-like sequences ? voir comment introduire ça car il faudrait parser le rep d'Extraction (en fonction du mode) donc preciser le input_ext, faire la liste des sequences et parser
