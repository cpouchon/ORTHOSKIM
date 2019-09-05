#!/usr/bin/env python

import glob
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import sys
import os, errno
import argparse



parser = argparse.ArgumentParser(description='Export regions of query database from exonerate alignment. Filtering is made according to rawscore. If there is multiple non overlapping parts of alignment, a \'partial\' tag is add to the header. Script was writen by C. Pouchon (2019).')
parser.add_argument("-i","--infile", help="Fasta input file from exonerate run",
                    type=str)
parser.add_argument("-t","--table", help="out table of diamond run",
                    type=str)
parser.add_argument("-m","--model", help="molecular type position",
                    type=str,choices=["chloroplast", "mitochondrion","nucleus"])
parser.add_argument("-o","--outdir", help="Out directory path",
                    type=str)
parser.add_argument("-n","--namesample", help="sample name used in sequence header for output",
                    type=str)
parser.add_argument("-c","--maxcont", help="maximal threshold of contigs number allowed to match on reference sequence",
                    type=int)
parser.add_argument("-l","--minlength", help="minimal length of alignment into reference to consider output",
                    type=int)
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
                print ("path '%s' already exists" % path)   # overwrite == False and we've hit a directory that exists
        else: raise

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

def len_overlapping(x, y):
    return len(list(set(x) & set(y)))

'#########################################################'
'function'
'#########################################################'

input_file = args.infile
model = args.model
outdir = args.outdir
max_contigs = args.maxcont
nameofsample = args.namesample
min_length = args.minlength
diamond_tab = args.table

with open(diamond_tab) as f:
    tab = f.readlines()
hits = {}
for line in tab:
    l = line.rstrip().split("\t")
    seqid = l[0]
    geneid = l[1].split("_")[0]
    if geneid not in hits.keys():
        hits[geneid] = list()
    if seqid not in hits[geneid]:
        hits[geneid].append(seqid)

mkdir(outdir)

stored = {}
cur_genome = SeqIO.parse(input_file, "fasta")
for record in cur_genome:
    info = record.description.split(" ")[1:len(record.description.split(" "))]
    dicinfo = dict(item.split("=") for item in info)
    dicinfo['seq']=str(record.seq)
    gene_ref = dicinfo['gene'].split("_")

    if gene_ref[0] not in hits.keys():
        'we check if a refgene is mapped in exonerate but not in diamond and keep all seqID mapping the gene'
        stored.setdefault(gene_ref[0], []).append(dicinfo)
    else:
        'if the gene is mapped in diamond AND exonerate, we keep only the seqID mapping in both softs'
        if dicinfo['ref'] in hits[gene_ref[0]]:
            stored.setdefault(gene_ref[0], []).append(dicinfo)
        else:
            pass

Bar = ProgressBar(len(stored.keys()), 60, '\t Extraction of genes for %s' % input_file)
barp=0
for geneid in stored.keys():
    'we seek if there is overlapping position according to the reference in records for the gene'
    indexint = range(len(stored[geneid]))
    groups = {}
    for x in indexint:
        loverlap=list()
        for y in range(x,len(indexint)):
            range_x = range(int(stored[geneid][x]['qbal']),int(stored[geneid][x]['qeal']))
            range_y = range(int(stored[geneid][y]['qbal']),int(stored[geneid][y]['qeal']))
            if len_overlapping(range_x, range_y) != 0:
                loverlap.append(x)
                loverlap.append(y)
            else:
                pass
        'we create/append a dict with overlapping groups as key'
        if len(set(loverlap))>0:
            if len(groups)==0:
                groups[len(groups)+1]= list(set(loverlap))
            else:
                'if the dict exist we check if we retrieve groups'
                cond = 0
                for i in range(len(groups)):
                    if len(set(loverlap)&set(groups[groups.keys()[i]])) > 0:
                        'if we retriewe group we add to the corresponding dict[key] the new record for the overlapping group'
                        temp = [x for x in loverlap if x not in set(groups[groups.keys()[i]])]
                        groups[groups.keys()[i]].extend(temp)
                        cond = cond+1
                    else:
                        pass
                if cond == 0:
                    'if not retrieved we create a new key with a new group on the dict'
                    groups[len(groups)+1]= list(set(loverlap))

    'now we have to check if we have multiple key on the dict --> we have multiple overlapping regions for a same gene'

    'if not subgroup we select sequences based on score (1) and length (2)'

    if len(groups)==1:
        lscore=list()
        llen=list()
        passFilter1=list()
        passFilter2=list()
        tmp_score=list()
        contigs_count=0
        count=0
        contignames=list()
        count_parts=list()
        sbindex=groups[groups.keys()[0]]
        'add a condition to test if there is multiple contigs matching the same region (indication of duplicated area)'
        for i in range(len(sbindex)):
            contignames.append(stored[geneid][sbindex[i]]['ref'])
        contigs_count = contigs_count+len(set(contignames))
        count=count+len(set(contignames))

        'We retrieve all alignment score for each unique gene'
        for i in range(len(stored[geneid])):
            llen.append(int(stored[geneid][i]['alnlen']))
        'we extract dict of our list with max(score)'
        for item in stored[geneid]:
            if int(item['alnlen']) == max(llen):
                passFilter1.append(item)

        'we check if there is multiple sequences with the highest score'
        if len(passFilter1) > 1:
            'if there is multiple highest score we check the length of the alignment'
            for i in range(len(passFilter1)):
                lscore.append(int(passFilter1[i]['score']))
            'if there is multiple highest score with the same length we choose the first record'
            for item in passFilter1:
                if int(item['score']) == max(lscore):
                    tmp_score.append(item)
            passFilter2 = tmp_score[0]
        else:
            passFilter2 = passFilter1[0]

        count_parts.append(count)
        aln_length= int(passFilter2['qlen'])*3
        aln_seq=list("N"*aln_length)
        bestid=passFilter2["ref"]

        for i in range(len(stored[geneid])):
            if stored[geneid][i]['ref'] == bestid:
                if int(stored[geneid][i]['qbal'])<int(stored[geneid][i]['qeal']):
                    sstart = int(stored[geneid][i]['qbal'])
                    send = int(stored[geneid][i]['qeal'])
                else:
                    sstart = int(stored[geneid][i]['qeal'])
                    send = int(stored[geneid][i]['qbal'])
                al_start = (int(sstart+1)*3)-2
                al_end = int(send)*3
                aln_seq[al_start-1:al_end] = list(stored[geneid][i]['seq'])

        out_seq="".join(aln_seq)
        catseq=out_seq.replace("N","")
        newlength = len(catseq)

        if "UCE" in stored[geneid][i]["gene"]:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"_UCE"+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+" match_contigs="+str(contigs_count)
            nuctype="_UCE"
        elif "BUSCO" in stored[geneid][i]["gene"]:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"_BUSCO"+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+" match_contigs="+str(contigs_count)
            nuctype="_BUSCO"
        else:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; "+" match_contigs="+str(contigs_count)
            nuctype=""

        cond_min_length=0
        if newlength>=min_length:
            cond_min_length=cond_min_length+1

    else:
        'we aim to concatenate sequences of genes if we have multiple partial part'
        MultipleParts=list()
        contigs_parts=list()
        count_parts=list()
        for subgp in range(len(groups)):
            sbindex=groups[groups.keys()[subgp]]
            lscore=list()
            llen=list()
            sbgpassFilter1=list()
            sbgpassFilter2=list()
            subtmp_score=list()
            contigs_count=0
            count=0
            contignames=list()
            'add a condition to test if there is multiple contigs matching the same region (indication of duplicated area)'
            for i in range(len(sbindex)):
                contignames.append(stored[geneid][sbindex[i]]['ref'])
            contigs_count = str(contigs_count+len(set(contignames)))
            count=count+len(set(contignames))

            'We retrieve all alignment score for each unique gene'
            for i in range(len(sbindex)):
                llen.append(int(stored[geneid][sbindex[i]]['alnlen']))
            'we extract dict of our list with max(score)'
            for i in range(len(sbindex)):
                if int(stored[geneid][sbindex[i]]['alnlen']) == max(llen):
                    sbgpassFilter1.append(stored[geneid][sbindex[i]])

            # 'We retrieve all alignment score for each unique gene'
            # for i in range(len(sbindex)):
            #     lscore.append(stored[geneid][sbindex[i]]['score'])
            # 'we extract dict of our list with max(score)'
            # for i in range(len(sbindex)):
            #     if stored[geneid][sbindex[i]]['score'] == max(lscore):
            #         sbgpassFilter1.append(stored[geneid][sbindex[i]])

            # 'we check if there is multiple sequences with the highest score'
            # if len(sbgpassFilter1) > 1:
            #     'if there is multiple highest score we check the length of the alignment'
            #     for i in range(len(sbgpassFilter1)):
            #         llen.append(sbgpassFilter1[i]['alnlen'])
            #     'if there is multiple highest score with the same length we choose the first record'
            #     if llen.count(max(llen)) > 1:
            #         sbgpassFilter2 = sbgpassFilter1[0]
            #     else:
            #         sbgpassFilter2=(next(x for x in sbgpassFilter1 if x['alnlen'] == max(llen)))
            # else:
            #     sbgpassFilter2 = sbgpassFilter1[0]

            'we check if there is multiple sequences with the highest score'
            if len(sbgpassFilter1) > 1:
                'if there is multiple highest score we check the length of the alignment'
                for i in range(len(sbgpassFilter1)):
                    lscore.append(int(sbgpassFilter1[i]['score']))
                'if there is multiple highest score with the same length we choose the first record'
                for item in sbgpassFilter1:
                    if int(item['score']) == max(lscore):
                        subtmp_score.append(item)
                sbgpassFilter2 = subtmp_score[0]
            else:
                sbgpassFilter2 = sbgpassFilter1[0]

            MultipleParts.append(sbgpassFilter2)
            contigs_parts.append(contigs_count)
            count_parts.append(count)

        'we check the order parts of gene according to the reference'

        aln_length= int(MultipleParts[0]['qlen'])*3
        aln_seq=list("N"*aln_length)

        for subgp in range(len(groups)):
            bestid=MultipleParts[subgp]["ref"]
            sbindex=groups[groups.keys()[subgp]]

            for i in range(len(sbindex)):
                if stored[geneid][sbindex[i]]['ref'] == bestid:
                    sstart = int(stored[geneid][sbindex[i]]['qbal'])
                    send = int(stored[geneid][sbindex[i]]['qeal'])
                    'we define start/end of the sequence according to the reference in protein'
                    'sstart+1 because in exonerate, alignment starts at +1 (ex: bal=0 or it starts to 1)'
                    al_start = (int(sstart+1)*3)-2
                    al_end = int(send)*3
                    'we put -1 on al_start because of python index / and no for the end because of not inclusion of the last value Y in X:Y'
                    aln_seq[(al_start)-1:(al_end)] = list(stored[geneid][sbindex[i]]['seq'])
                else:
                    pass

        out_seq="".join(aln_seq)
        catseq=out_seq.replace("N","")
        newlength = len(catseq)
        if "UCE" in MultipleParts[subgp]["gene"]:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"_UCE"+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; info=concatenated_parts"+"; "+" match_contigs="+str("-".join(contigs_parts))
            nuctype="_UCE"
        elif "BUSCO" in MultipleParts[subgp]["gene"]:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"_BUSCO"+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; info=concatenated_parts"+"; "+" match_contigs="+str("-".join(contigs_parts))
            nuctype="_BUSCO"
        else:
            header = ">"+str(nameofsample)+"; "+"gene="+str(geneid)+"; "+"type="+str(model)+"; "+"length="+str(newlength)+"; info=concatenated_parts"+"; "+" match_contigs="+str("-".join(contigs_parts))
            nuctype=""

        cond_min_length=0
        if newlength>=min_length:
            cond_min_length=cond_min_length+1


    'Output results: we write/append fasta file with sequence taxa for each gene'
    fname = geneid+str(nuctype)+".fa"

    cond_max=0
    'we check the number of contigs matching each overlapping parts'
    for cont in range(len(count_parts)):
        if count_parts[cont] > max_contigs:
            cond_max=cond_max+1
        else:
            pass

    if cond_max==0 and cond_min_length!=0:
        if os.path.isfile(os.path.join(outdir, fname)):
            with open(os.path.join(outdir, fname), 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.rstrip())

                if not header in old_headers:
                    Bar.update(barp)
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(out_seq)+'\n')
                else:
                    pass
        else :
            with open(os.path.join(outdir, fname), 'w') as out:
                Bar.update(barp)
                out.write(header+'\n')
                out.write(str(out_seq)+'\n')

    else:
        continue
