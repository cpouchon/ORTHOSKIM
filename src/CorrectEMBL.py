#! /usr/bin/env python
import sys
import os, errno
import argparse

parser = argparse.ArgumentParser(description='Correction of embl files for tRNA features in the end of contigs')
parser.add_argument("-f","--file", help="embl file to modify",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

embl_file=args.file

with open(embl_file) as f:
    lines = f.readlines()

l_to_remove=list()
cond_tRNA=0
for l in range(len(lines)):
    line=lines[l]
    if line.startswith("FT"):
        if line.split(" ")[3]!="":
            if line.split(" ")[3]=="tRNA":
                cond_tRNA=cond_tRNA+1
                l_tRNA_start=l
            else:
                cond_tRNA=0
        else:
            pass
    elif line.startswith("XX"):
        if cond_tRNA>0:
            l_tRNA_end=l-1
            for i in list(range(l_tRNA_start,l_tRNA_end+1)):
                l_to_remove.append(i)
            cond_tRNA=0
        else:
            pass
    else:
        pass

for l in range(len(lines)):
    line=lines[l].rstrip()
    if l not in l_to_remove:
        print(line)
