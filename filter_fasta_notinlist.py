## This script will filter sequences in a fasta given a list of sequence names.  
## Unlike the "filter_fasta.py" script, this script prints those not in the list
##
## IW 30-10-2023

import sys
import os
import argparse
from Bio import SeqIO

fasta_file = sys.argv[1]
list_file = sys.argv[2]


ids = []

with open(list_file) as lis:
    for row in lis:
        ids.append(row.strip())


with open(fasta_file) as fas:

    for record in SeqIO.parse(fas, "fasta"):
            
        if record.id not in ids:
                
            print('>' + record.id + '\n' + record.seq)
                
                
