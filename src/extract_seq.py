#!/usr/bin/env python3 
# -*- coding: utf-8 -*

import sys,os,re, time, csv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC 
from Bio.SeqRecord import SeqRecord

mito=sys.argv[1]
conta=sys.argv[2]
output=sys.argv[3]

nf=open(output,"w")

list_id=[]

ensemble_SeqRecord1=SeqIO.parse(mito,"fasta")
for seq_record1 in ensemble_SeqRecord1:
    if seq_record1.id not in list_id:
        list_id.append(seq_record1.id)
        nf.write(">" + str(seq_record1.description) + "\n")
        nf.write(str(seq_record1.seq) + "\n")

ensemble_SeqRecord2=SeqIO.parse(conta,"fasta")
for seq_record2 in ensemble_SeqRecord2:
    if seq_record2.id not in list_id:
        list_id.append(seq_record2.id)
        nf.write(">" + str(seq_record2.description) + "\n")
        nf.write(str(seq_record2.seq) + "\n")

nf.close()