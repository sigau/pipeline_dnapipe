#!/usr/bin/env python3 
# -*- coding: utf-8 -*

import sys,os,re, time, csv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC 
from Bio.SeqRecord import SeqRecord

f_seq=sys.argv[1]
tab_cores=sys.argv[2]
specie=sys.argv[3]
acc_num=sys.argv[4]
output=sys.argv[5]

####create a dictionary mapped sequence : ref sequence
cores_dict={}
with open(tab_cores,"r") as table:
    for line in table:
        seq=line.split()[0]
        ref=line.split()[1]
        cores_dict[seq]=ref

### rename the sequence 
nf=open(output,"w")
ensemble_SeqRecord=SeqIO.parse(f_seq,'fasta')
for seq in ensemble_SeqRecord:
    orig_desc=seq.description
    new_desc=f"{specie}#{acc_num}_{orig_desc}@{cores_dict[orig_desc]}"
    nf.write(f">{new_desc}\n{seq.seq}\n") 

nf.close()   