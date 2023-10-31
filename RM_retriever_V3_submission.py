#!/usr/bin/env python3

### This aim of this script is to retrieve total BPs of different repeats in the genome based on RepeatMasker
### outputs. It splits these outputs into Full copies (Overlapping with full repeat positions)
### and Solo-LTRs+fragments. Overlapping works a bit differently too, for internal variations any overlap is
### counted towards the full element. For LTRs, only the overlap of starting position and ending position is counted
### towards the full element (This allows counting for solo-LTRs insertions inside of the copy). 

import sys
import os
import argparse
from Bio import SeqIO
import pysam

parser = argparse.ArgumentParser(description="Counts repeatmasker hits of crapaud LTR variants")

# Add the arguments to the parser
parser.add_argument("-rf", "--rfile", dest="path_to_rmfiles", required=True,
                    help="Directory with the Repeatmasker output")
parser.add_argument("-cf","--classified_file", dest = "classified_file",required=True,
                    help = "File containing the positions and classifications and solo or full information of the 150 classified LTR")
parser.add_argument("-m","--minsize",default=250, dest = "minsize_count",type=int,
                    help = "minimum size considered as an LTR hit")
parser.add_argument("-s","--stitch",default=False, action ='store_true',
                    help = "If stitching of small hits should be done for the counting")                    
parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output file. Returns the table of LTR counts")
args = parser.parse_args()

## List of the Podospora species
species_list = ["PaWa63p","PcWa139m","CBS124.78p","CBS237.71m","CBS411.78m","CBS415.72m","CBS112042p"]


# Calculates overlap between two ranges and returns the length of the overlap

def overlap(start1, end1, start2, end2):
    y = range(start1,end1)
    x = range(start2, end2) 

    """how much does the range (start1, end1) overlap with (start2, end2)"""
    r = range(max(x[0],y[0]),min(x[-1],y[-1]+1))
    return len(r)
    #return max(max((end2-start1), 0) - max((end2-end1), 0) - max((start2-start1), 0), 0)


## uses the package pysam to retrieve fasta sequences based on chromosome and positions

def retriever(fasta,chrom, start, end):

    seq = pysam.FastaFile(fasta)
    cur_seq = seq.fetch(chrom, start, end)
    
    return cur_seq



full_ids = {}
error_count = {}

## Loads in the full copy fasta files, they have coordinates in the names
# which are used to determine overlaps with rm-hits.  
fulltab = [] 
solotab = []
full_count = 0
solo_count = 0

with open(args.classified_file) as classified:

    for line in classified:
        
        
            
        seq, _, LTR, cl = line.strip().split('\t')

        spec = seq.split('_')[-4]
        chromo1 = seq.split('_')[-3]
        chromo2 = seq.split('_')[-2]
        pos = seq.split('_')[-1]


        start, end = pos.split('-')
            
        if cl == "Full":
            fulltab.append((spec.strip('>'),f"{chromo1}_{chromo2}",int(start),int(end), LTR))
        elif cl == "Solo":
            solotab.append((spec.strip('>'),f"{chromo1}_{chromo2}",int(start),int(end), LTR))



for s in species_list:

    rmtab = []
    filter_rmtab=[]
    stitch_tab =[]
    


    ## Loading in the repeatmasker table

    with open(f"{args.path_to_rmfiles}{s}.nice.fa.ori.out", 'r') as inp:

        for i, line in enumerate(inp):

            rm = line.strip().split()

            rmtab.append((rm[4],int(rm[5]),int(rm[6]),rm[8],rm[9]))

    for hit in rmtab:

        hitsize = max(hit[1],hit[2])-min(hit[1],hit[2])

        if hitsize >= args.minsize_count:

            filter_rmtab.append(hit)
    
        elif args.stitch:

            stitch_tab.append(hit)



    ## Sets the LTR counts as 0 to begin.

    fullBP = {'LTR1':0,'LTR2':0,'LTR3':0,'LTR4':0,'LTR5':0,'LTR6':0,'LTR7':0,'LTR8':0,'LTR9':0,'LTR10':0,'LTR11':0,'LTR12':0,'LTR13':0,'LTR14':0}
    solofragBP = {'LTR1':0,'LTR2':0,'LTR3':0,'LTR4':0,'LTR5':0,'LTR6':0,'LTR7':0,'LTR8':0,'LTR9':0,'LTR10':0,'LTR11':0,'LTR12':0,'LTR13':0,'LTR14':0}
    fullLTR = {'LTR1':0,'LTR2':0,'LTR3':0,'LTR4':0,'LTR5':0,'LTR6':0,'LTR7':0,'LTR8':0,'LTR9':0,'LTR10':0,'LTR11':0,'LTR12':0,'LTR13':0,'LTR14':0}
    soloLTR = {'LTR1':0,'LTR2':0,'LTR3':0,'LTR4':0,'LTR5':0,'LTR6':0,'LTR7':0,'LTR8':0,'LTR9':0,'LTR10':0,'LTR11':0,'LTR12':0,'LTR13':0,'LTR14':0}

    fastatab = []
    passed_ids = []
    lines = 0

    for line in filter_rmtab:

        species, _,_ = line[0].split('_')


        if 'LTR' in str(line[4]):
                
            for rep in fulltab:
                
                if line[0] == f"{rep[0]}_{rep[1]}":
                    
                    
                    if overlap(line[1],line[2], min(rep[2],rep[3])-20,min(rep[2],rep[3])+400) > 300 or overlap(line[1],line[2],max(rep[2],rep[3])-400,max(rep[2],rep[3])+20) > 300:    
                                                
                        if line[4] in fullBP.keys():
                            fullBP[line[4]] = fullBP[line[4]] + abs(line[2]-line[1])
                        else:
                            fullBP[line[4]] = abs(line[2]-line[1])
                        
                        for LTR in fullLTR.keys():

                            if f"{LTR}_" in line[4]:
                                
                                fullLTR[LTR] += 1
                                lines +=1
                                passed_ids.append(line)

                                if f"{rep[0]}_{rep[1]}_{rep[2]}-{rep[3]}" in full_ids:
                                    full_ids[f"{rep[0]}_{rep[1]}_{rep[2]}-{rep[3]}"].append(line)
                                else:
                                    full_ids[f"{rep[0]}_{rep[1]}_{rep[2]}-{rep[3]}"] = [line]

                                    fastatab.append((f">full_{line[0]}_{line[1]}_{line[2]}_{line[4]}", retriever(f"/home/ivar/Crapaud/data/PodosporaGenomes/{species}.nice.fa",line[0], min(line[2],line[1]),max(line[2],line[1]))))
                                
                                ## Checks if the classified LTR matches the RM hit
                                if rep[4] != LTR and rep[4] != "unclassified":
                                    if f"{rep[4]}_{LTR}" in error_count.keys():
                                        error_count[f"{rep[4]}_{LTR}"] += 1
                                    else:
                                        error_count[f"{rep[4]}_{LTR}"] = 1
            # All full copy LTRs has been added to the passed_ids variable, so here it takes all other hits, ie the solo/fragments    
            if line not in passed_ids:

                for LTR in soloLTR.keys():

                    if f"{LTR}_" in line[4]:

                        soloLTR[LTR] += 1
                        lines += 1
                        fastatab.append((f">solo_{line[0]}_{line[1]}_{line[2]}_{line[4]}", retriever(f"/home/ivar/Crapaud/data/PodosporaGenomes/{species}.nice.fa",line[0], min(line[2],line[1]),max(line[2],line[1]))))

                        for cl in solotab:
                            
                            if overlap(line[1],line[2],min(cl[2],cl[3]),max(cl[2],cl[3])) > 0.9*(max(cl[2],cl[3])-min(cl[2],cl[3])):

                                if cl[4] != LTR and cl[4] != "unclassified":
                                    if f"{cl[4]}_{LTR}" in error_count.keys():
                                        error_count[f"{cl[4]}_{LTR}"] += 1
                                    else:
                                        error_count[f"{cl[4]}_{LTR}"] = 1

        #If the library also contains inner region fastas
        else: 
                
            for rep in fulltab:

                if line[0] == f"{rep[0]}_{rep[1]}" and overlap(line[1],line[2], min(rep[2],rep[3]), max(rep[2],rep[3])) != 0:

                    if line[4] in fullBP.keys():
                        fullBP[line[4]] = fullBP[line[4]] + abs(line[2]-line[1])
                    else: 
                        fullBP[line[4]] = abs(line[2]-line[1])
                else:
                    if line[4] in solofragBP.keys():
                        solofragBP[line[4]] = solofragBP[line[4]] + abs(line[2]-line[1])
                    else: 
                        solofragBP[line[4]] = abs(line[2]-line[1])

    
    ## Finds full elements that has three rmhits assigned to it, then if there are two of one variant and
    ## one of another it will reassign the lone one to a solo in the counting
    
    for id in full_ids:
        splitid = id.split('_')
        if splitid[0] == s:
            if len(full_ids[id]) > 2:
                print("yo")
            if len(full_ids[id]) == 3:    
                if full_ids[id][0][3] == full_ids[id][1][3] and full_ids[id][0][4] == full_ids[id][1][4]:
                    if full_ids[id][0][3] != full_ids[id][2][3] or full_ids[id][0][4] != full_ids[id][2][4]:
                        if species == full_ids[id][2][0]:
                            fullLTR[full_ids[id][2][4]] = fullLTR[full_ids[id][2][4]] -1
                            soloLTR[full_ids[id][2][4]] += 1
                elif full_ids[id][1][3] == full_ids[id][2][3] and full_ids[id][1][4] == full_ids[id][2][4]:
                    if full_ids[id][0][3] != full_ids[id][1][3] or full_ids[id][0][4] != full_ids[id][1][4]:
                        if species == full_ids[id][0][0]:
                            fullLTR[full_ids[id][0][4]] = fullLTR[full_ids[id][0][4]] -1
                            soloLTR[full_ids[id][0][4]] += 1
                elif full_ids[id][0][3] == full_ids[id][2][3] and full_ids[id][0][4] == full_ids[id][2][4]:
                        if full_ids[id][0][3] != full_ids[id][1][3] or full_ids[id][0][4] != full_ids[id][1][4]:
                            if species == full_ids[id][1][0]:
                                fullLTR[full_ids[id][1][4]] = fullLTR[full_ids[id][1][4]] -1
                                soloLTR[full_ids[id][1][4]] += 1
                

    


    with open(args.file_out,'a') as o:
        for count in fullLTR:
            o.write(f"{count}\t{species}\tfull\t{fullLTR[count]}\n")
        for count in soloLTR:
            o.write(f"{count}\t{species}\tsolo\t{soloLTR[count]}\n")

    full_count += sum(fullLTR.values())
    solo_count += sum(soloLTR.values())


print(f"Finished running, in total {full_count} full copies and {solo_count} solo/fragments were counted after filtering\nThis is {(full_count+solo_count)-1079} more than previously classified")
print(f"A more detailed count table can be found at: {args.file_out}",'\n')
print('\n',f"The error count compared to the classified full and solo/fragment copies was: {sum(error_count.values())}/1079, or {sum(error_count.values())/1079}",'\n')
print(f"Detailed information on error count can be found in: {args.file_out}.counts\n")

count_perc = sum(error_count.values())/1079

with open(f"{args.file_out}.count","w") as o2:
    o2.write(f"The error count compared to the classified full and solo/fragment copies was: {sum(error_count.values())}/1079, or {count_perc}\n")
    o2.write("\nCount table, structured 'classifiedLTR_RMclassifiedLTR: count'\n")
    for count in error_count:
        o2.write(f"{count}\t{error_count[count]}\n")

