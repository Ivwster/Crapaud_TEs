## A simple script to get the GC-content of fasta sequences in a multi-fasta file
#
# IW 30-10-2023

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import os


input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "rU") as file:
    
    with open(output_file, "w") as out:
        for record in SeqIO.parse(file, "fasta"):

            cur_gc = GC(record.seq)

            print(f'{record.id}\t{cur_gc}')
            out.write(f"{record.id}\t{cur_gc}\n")