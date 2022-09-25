#
# -------------------------------------------------------------------------
#                              alpha_com.py
# --------------------------------------------------------------------------
#
# Monte Carlo simulator of a system of spins that can take 3 values 
# corresponding to the vectors (-1,0),(0,alpha) and (1,0) 
# respectively, where alpha is a real parameter. Spins are located in a
# pregenerated network with Newmann  community structure, and energy is
# defined as the negative sum of dot products over linked nodes.
# From a random initial opinion configuration, the program proposes 
# sequencial changes and accepts those that minimize the energy following 
# a Boltzmann distribution at a given range of temperatures usinc a Monte-Carlo
# Metropolis dynamics.
# Averages over 'imctot' final states and variances of the energy and the
# proportions of each possible value are calculated, written to a file 
# and printed on the screen for every temp.
#
from __future__ import division
import numpy as np
from networkx import nx
import time
import matplotlib.pyplot as plt
import pandas as pd
#
# ------------------------------------------------------------------------------
#   DATA
# ------------------------------------------------------------------------------
#
start=time.time()      # Start measuring time
com_file = 'leon.net'
#com_file = 'm3n12.net'
G = nx.read_edgelist(com_file,nodetype=int) # Reads the community from the input file and creates the graph 
number_nodes = G.number_of_nodes()
# -----------------------------------------------------------------------------

#  PLOTTING -------------------------------------------------------------------
#
label_fontsize = 20
tick_fontsize = 20
legend_fontsize = 15
#
fig, ax1 = plt.subplots(1, figsize=(5,5))
# Plot the orientation of the last config
pos=nx.spring_layout(G)
node_s = 2
w = 0.5
edge_w = 0.1

#rednodes = [node for node in G.nodes() if G.nodes[node]['orientation'] == -1]
#greennodes = [node for node in G.nodes() if G.nodes[node]['orientation'] == 0]
#bluenodes = [node for node in G.nodes() if G.nodes[node]['orientation'] == 1]

#for i in range (1,N+1):
#  if G.nodes[i]['orientation'] == -1:   # RED
#    nx.draw_networkx_nodes(G,pos,nodelist=[i],node_color='red',node_size=node_s, ax=ax1 )
#  elif G.nodes[i]['orientation'] == 0:    # GREEN 
#    nx.draw_networkx_nodes(G,pos,nodelist=[i],node_color='#00ff00',node_size=node_s,ax=ax1 )
#  else:                           # BLUE
#    nx.draw_networkx_nodes(G,pos,nodelist=[i],node_color='blue',node_size=node_s,ax=ax1 ) 

ax1.set_aspect("equal")
nx.draw_networkx_nodes(G, pos=pos,node_size=node_s )

nx.draw_networkx_edges(G,pos,width=edge_w)
# -----------------------------------------------------------------
# -----------------------------------------------------------------
plt.tight_layout()
plt.show()
