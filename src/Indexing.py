#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch



parser = argparse.ArgumentParser(description='creation of the sample table for phyloskims users')
parser.add_argument("-i","--infile", help="input file list of assembled annotated chloroplasts",
                    type=str)
parser.add_argument("-p","--path", help="searching path of annotated chloroplasts",
                    type=str)
parser.add_argument("-cfind","--chlorofind", help="[mode] search all annotated chloroplast in given path (-p) otherwise parse a given list (-i)",
                    action="store_true")
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)



args = parser.parse_args()


if args.chlorofind:
    'Search all annotated chloroplasts files'
    path = args.path
    in_files = []
    for r, d, f in os.walk(path):
        for file in f:
            if 'indexing.done' in file:
                if os.path.join(r, file) not in in_files:
                    in_files.append(os.path.join(r,file))
else:
    'We parse a list of files'
    input_list = args.infile
    with open(input_list) as f:
        in_files = f.readlines()


for file in in_files:
    try:
        cf = file.rstrip()
        cfile_path = cf
        cfile_name = os.path.basename(cfile_path)

        cf_parts = cf.split("/")

        code=cf_parts[-3].replace(":","_")
        sequencing = cf_parts[-2].replace(":","_")
        taxainfos=cf_parts[-4]
        ctaxaname=taxainfos.split(":")[0]
        taxid=taxainfos.split(":")[1]
        #ingenus=ctaxaname.split("_")[0]
        samplename=str(ctaxaname+"_"+taxid+"_"+code+"_"+sequencing)
        embl_chloro=".".join([taxainfos,cf_parts[-3],cf_parts[-2],"chloro","embl"])

        if "GWM" in cf:
            r1 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_R1.fastq.gz')
            r1.sort()
            r2 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_R2.fastq.gz')
            r2.sort()
            r1path=os.path.dirname(cf)+"/"+r1[0]
            r2path=os.path.dirname(cf)+"/"+r2[0]
        elif "RSZ" in cf:
            r1 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_1.fq.gz')
            r1.sort()
            r2 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_2.fq.gz')
            r2.sort()
            r1path=os.path.dirname(cf)+"/"+r1[0]
            r2path=os.path.dirname(cf)+"/"+r2[0]
        elif "AXZ" in sequencing.split("_")[0]:
            r1 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_[0-9]_1_*.fastq.*')
            r1.sort()
            r2 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_[0-9]_2_*.fastq.*')
            r2.sort()
            r1path=os.path.dirname(cf)+"/"+r1[0]
            r2path=os.path.dirname(cf)+"/"+r2[0]
        else:
            r1 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_[0-9]_1_*_clean.fastq.gz')
            r1.sort()
            r2 = fnmatch.filter(os.listdir(os.path.dirname(cf)), '*_[0-9]_2_*_clean.fastq.gz')
            r2.sort()
            r1path=os.path.dirname(cf)+"/"+r1[0]
            r2path=os.path.dirname(cf)+"/"+r2[0]

        print ('%s\t%s\t%s\t%s' % (str(samplename),str(r1path),str(r2path),str(cfile_path.replace("indexing.done",embl_chloro))))

    except:
        pass
        with open("./indexing_error.log", 'a+') as file:
            file.write(cfile_path+'\n')
