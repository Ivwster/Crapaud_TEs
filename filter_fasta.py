## A script to filter a fasta file given a list with sequence names, the output fasta will only contain the listed sequences
## 
## IW 30-10-2023

import sys
import os
from Bio import SeqIO

fasta_file = sys.argv[1]
list_file = sys.argv[2]

ids = []

with open(list_file) as lis:
    for row in lis:
        ids.append(row.strip())


with open(fasta_file) as fas:

    for record in SeqIO.parse(fas, "fasta"):
            
        if record.id in ids:
                
            print('>' + record.id + '\n' + record.seq)
                
                
