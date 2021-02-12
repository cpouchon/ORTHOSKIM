#! /usr/bin/env python

import sys
import os, errno
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-t","--taxa", help="taxa list",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path",
                    type=str)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()


def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.
    Args:
        list_of_lengths (list): List of numbers.
    Returns:
        float: N50 value.
    """
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
    return median


path = args.path
outdirect = args.outdir
giventaxa = args.taxa

with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    #taxa_list.append(l)
    taxa_dict.setdefault(l, [])

for taxa in taxa_dict.keys():
    cont_size=list()
    seq_all=list()
    count_GC=0
    try:
        with open(os.path.join(str(path), str(taxa+".fa")), 'r') as out:
            for l in out:
                if l.startswith(">"):
                    cont=l.rstrip().split("_")
                    cont_size.append(int(cont[3]))
                else:
                    count_GC=count_GC+l.rstrip().count("G")+l.rstrip().count("C")
        contnumber=len(cont_size)
        totalsize=sum(cont_size)
        N50=calculate_N50(cont_size)
        L50=sum(i<=sum(cont_size)/2 for i in cont_size)
        GCcontent=round(count_GC/sum(cont_size)*100,2)
    except:
        contnumber=0
        totalsize=0
        N50=0
        L50=0
        GCcontent=0
    toprint=list()
    toprint.append(taxa)
    toprint.append(str(contnumber))
    toprint.append(str(totalsize))
    toprint.append(str(N50))
    toprint.append(str(L50))
    toprint.append(str(GCcontent))
    with open(os.path.join(str(outdirect),'assemblies_statistics.txt'), 'a+') as out:
        out.write("\t".join(toprint)+'\n')
        out.close()
