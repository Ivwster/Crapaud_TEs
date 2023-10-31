# Crapaud_TEs




The first  three scripts are modifications on the from  https://github.com/MiguelMSandin/SSNetworks that works better for fragmented hits, such as Transposable Elements:

**1.1_blastn_allAgainstAll.sh**, which does an all-vs-all blast of the sequences in a fasta file. It is modified from Miguel's scrips to use blastn instead of megablast

**1.2_blastClean_modIW_v2.py**, which cleans up the blast output by removing self-hits and reciprocal hits. In my modified version it also deals with fragmented and overlapping hits which you get with blastn of transposable elements.

**2.1_buildNetwork_IWmod**, which builds the network by considering coverage and percent identity. It is modified to calculate it from the output of the modified scrips.


**Network_x.py**, utilizes the package networkx to calculate communities in the network generated bythe three first scripts


**RM_retriever_V3_submission.py**, Calculates number of RepeatMasker hits belonging to different subfamilies of _crapaud_


The three R-scripts wereused for plottingand statistical analysis ofthe data:

**TE_abundance.R**, Used for the calculation and analyses of TE-abundances

**rogue_stats.R**, Used for the analysis of _crapaud_ subfamily copies falling outside of their main clade,i.e rogue copies. 

**RepeatMasker_barplots_SSN-V2.R**, Used for plotting and analysis of the outputfromthe RM_retriever_V3_submission.py script


The last scripts were usedfor various small functions throughout the study:

**filter_fasta.py**, used to filtera fasta file using a list of sequence names, only keeping the sequences in the list

**filter_fasta_notinlist.py**, Similar to the filter_fasta.py script but instead only keeps the sequences not in the list. 

**GC_fullsequence.py**, Calculated the GC-content of sequences in a fasta file
