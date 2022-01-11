
# NETWORKX RUNNER v1.0
#
# What this file do?
# Look for average path length
#
# Code written by Zoltan Dul, PhD (2021)
# Contact me at zoltan dul [at] gmail.com
#
#

# Importing libraries

import csv
import sys
import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Maximalize memory read size
csv.field_size_limit(sys.maxsize)

protein_list = ["P31946", "P62258", "Q04917", "P61981", "P31947", "P27348", "P63104"]
proteins = '%0d'.join(protein_list)
url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=9606'
r = requests.get(url)

lines = r.text.split('\n')  # pull the text from the response object and split based on new lines
data = [l.split('\t') for l in lines]   # split each line into its components based on tabs

# convert to dataframe using the first row as the column names; drop empty, final row
df = pd.DataFrame(data[1:-1], columns = data[0])

# dataframe with the preferred names of the two proteins and the score of the interaction
interactions = df[['preferredName_A', 'preferredName_B', 'score']]


G = nx.Graph(name='Protein Interaction Test Graph')
interactions = np.array(interactions)
for i in range(len(interactions)):
    interaction = interactions[i]
    a = interaction[0] # protein a node
    b = interaction[1] # protein b node
    w = float(interaction[2]) # score as weighted edge where high scores = low weight
    G.add_weighted_edges_from([(a,b,w)]) # add weighted edge to graph


pos = nx.spring_layout(G) # position the nodes using the spring layout
plt.figure(figsize=(12, 12), facecolor=[0.7, 0.7, 0.7, 0.4])
nx.draw_networkx(G)
plt.axis('off')
plt.show()

