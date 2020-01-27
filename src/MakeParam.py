#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch



parser = argparse.ArgumentParser(description='Creation of parameters table for assembly path analysis according to annotated chloroplasts previously assembled.\
 Output format: (1) chloroplast file name; (2) Genus of species sample; (3) sample name; (4) forward reads, (5) reverse reads, (6) Output path for assemblies samples subdirectories\
 Script was writen by C. Pouchon (2019).')
parser.add_argument("-i","--infile", help="input file list of assembled annotated chloroplasts",
                    type=str)
parser.add_argument("-p","--path", help="searching path of annotated chloroplasts",
                    type=str)
parser.add_argument("-o","--outdir", help="out directory path to perform assembly for each sample",
                    type=str)
parser.add_argument("-cfind","--chlorofind", help="[mode] search all annotated chloroplast in given path (-p) otherwise parse a given list (-i)",
                    action="store_true")
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)



args = parser.parse_args()

outdir_assembly_path= args.outdir

if args.chlorofind:
    'Search all annotated chloroplasts files'
    path = args.path
    in_files = []
    for r, d, f in os.walk(path):
        for file in f:
            if '.chloro.embl' in file or '.rdnanuc.embl' in file:
                if os.path.join(r, "".join(file.split(".")[0:-2])) not in in_files:
                    in_files.append(os.path.join(r, "".join(file.split(".")[0:-2])))
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
        ingenus=ctaxaname.split("_")[0]
        samplename=str(ctaxaname+"_"+taxid+"_"+code+"_"+sequencing)

        #
        # if "subsp." in cfile_name:
        #     cfile_name_without_extension = cfile_name.split(".chloro.embl")[0]
        #     ctaxaname = cfile_name_without_extension.split(":")[0]
        #     taxid = cfile_name_without_extension.split(":")[1].split(".")[0]
        #     ingenus=ctaxaname.split("_")[0]
        #     sequencing = cfile_name.split(".")[3].split(":")
        #     samplename=str(ctaxaname+"_"+taxid+"_"+str(cfile_name.split(".")[2])+"_"+str(sequencing[0])+"_"+str(sequencing[1]))
        #
        # elif "var." in cfile_name:
        #     cfile_name_without_extension = cfile_name.split(".chloro.embl")[0]
        #     ctaxaname = cfile_name_without_extension.split(":")[0]
        #     taxid = cfile_name_without_extension.split(":")[1].split(".")[0]
        #     ingenus=ctaxaname.split("_")[0]
        #     sequencing = cfile_name.split(".")[3].split(":")
        #     samplename=str(ctaxaname+"_"+taxid+"_"+str(cfile_name.split(".")[2])+"_"+str(sequencing[0])+"_"+str(sequencing[1]))
        #
        # elif "sp." in cfile_name:
        #     cfile_name_without_extension = cfile_name.split(".chloro.embl")[0]
        #     ctaxaname = cfile_name_without_extension.split(":")[0]
        #     taxid = cfile_name_without_extension.split(":")[1].split(".")[0]
        #     ingenus=ctaxaname.split("_")[0]
        #     sequencing = cfile_name.split(".")[3].split(":")
        #     samplename=str(ctaxaname+"_"+taxid+"_"+str(cfile_name.split(".")[2])+"_"+str(sequencing[0])+"_"+str(sequencing[1]))
        #
        # else:
        #     cfile_name_without_extension = cfile_name.split(".chloro.embl")[0]
        #     ctaxaname = cfile_name_without_extension.split(":")[0]
        #     taxid = cfile_name_without_extension.split(":")[1].split(".")[0]
        #     ingenus=ctaxaname.split("_")[0]
        #     sequencing = cfile_name.split(".")[2].split(":")
        #     samplename=str(ctaxaname+"_"+taxid+"_"+str(cfile_name.split(".")[1])+"_"+str(sequencing[0])+"_"+str(sequencing[1]))

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


        print ('%s\t%s\t%s\t%s\t%s\t%s' % (str(str(cfile_path)+str(".chloro.embl")),str(ingenus),str(samplename),str(r1path),str(r2path),str(outdir_assembly_path)))

    except:
        pass
        with open("./indexing_error.log", 'a+') as file:
            file.write(cfile_path+'\n')
