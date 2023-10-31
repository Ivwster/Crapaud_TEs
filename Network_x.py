## This script runs the greedy modularity communities algorithm from networkx and then outputs what community/cluster the nodes fall in
#
# IW 30-10-2023

#!/usr/bin/env python3

import os
import sys
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities as gmc


in_net = sys.argv[1]



with open(in_net) as in_n:

    graph = nx.read_weighted_edgelist(in_n)

c = gmc(graph)

for i, community in enumerate(c):

    if len(community) > 4: 

        for seq in community:
            print(f"{seq}\tcluster_{i+1}")

        