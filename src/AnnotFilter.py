#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
#from Bio.Alphabet import *
from Bio.SeqRecord import *
#from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA;
from Bio.Seq import *
from Bio.SeqUtils import *
from ete3 import NCBITaxa


import numpy as np
import random
import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Filtering of annotations by taxonomy. Script was writen by C. Pouchon (2020).')
parser.add_argument("-i","--infile", help="input genomes annotation file",
                    type=str)
parser.add_argument("-f","--format", help="format of annotated file",choices=["embl", "genbank"],
                    type=str)
parser.add_argument("-l","--lineage", help="required lineage for single lineage (e.g. viridiplantae) or multiple lineages (e.g. asteraceae,helianthae)",
                    type=str)
parser.add_argument("-o","--outfile", help="outfile",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

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
		out = '\r %20s [%s%s] %3d %% ' % (self.title, '.' * bar, ' ' * (self.maxbar - bar), perc)
		sys.stdout.write(out)
		sys.stdout.flush()

outf=args.outfile
query_lineages = [str(item) for item in args.lineage.split(',')]
formatin = args.format

open(outf, 'w').close()


try:
    import ast
    import inspect
    import sys
    print("Patching NCBITaxa's base methods. For reason, see https://github.com/etetoolkit/ete/issues/469.\n")
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
#ncbi.update_taxonomy_database()

stored=list()

f=args.infile
file_path = f
print("Filtering annotations on taxonomy")
print('%s level(s) of taxonomy set: %s' % (len(query_lineages)," - ".join(query_lineages)))

#cur_genome =  SeqIO.index(f,formatin)
cur_genome=to_dict_remove_dups(SeqIO.parse(f,formatin))

Bar = ProgressBar(len(cur_genome), 60, '\t parsing annotations')
barp=0

for file in cur_genome:
    barp=barp+1
    Bar.update(barp)
    record=cur_genome[file]
    break_record=0
    for feat in record.features:
        if feat.type == 'source':
            if "/" in feat.qualifiers['organism'][0]:
                break
            else:
                org=feat.qualifiers['organism'][0].split(" ")
                for parts in feat.qualifiers['db_xref']:
                    if "taxon" in parts:
                        taxid=parts.split(":")[1]

                        if len(ncbi.get_taxid_translator([int(taxid)]))==0:
                            if len(ncbi.get_name_translator([" ".join(org)]))==0:
                                break
                            else:
                                taxid=list(ncbi.get_name_translator([" ".join(org)]).values())[0][0]
                                lineages=ncbi.get_lineage(taxid)
                                names = list(ncbi.get_taxid_translator(lineages).values())
                                lin_names = [name.lower() for name in names]
                        else:
                            lineages=ncbi.get_lineage(taxid)
                            names = list(ncbi.get_taxid_translator(lineages).values())
                            lin_names = [name.lower() for name in names]

                        for req_lin in query_lineages:
                            if str(req_lin).lower() in lin_names:
                                break_record=break_record+1
                            else:
                                pass
                    else:
                        break
        else:
            pass

    if break_record>0:
        stored.append(record)

print("\n%s / %s annotations selected on taxonomy"  % (len(stored),len(cur_genome)))

with open(outf, "a+") as output_handle:
    SeqIO.write(stored, output_handle, formatin)
