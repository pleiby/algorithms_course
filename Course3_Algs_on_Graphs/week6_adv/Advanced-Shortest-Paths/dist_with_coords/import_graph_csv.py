#!/usr/bin/python3

import pandas as pd

graphDataFilepath = "../../../data/san-francisco-california_graph.csv"

def readCSV(csvFilepath):
    try:
        csvdata = pd.read_csv(csvFilepath)
        return csvdata # returns a pandas dataframe
    except:
        print('yikes! Could not read file.')

ca_graph = readCSV(graphDataFilepath)

ca_graph.head()

#ca_graph.columns
ca_graph.describe()
ca_graph.dtypes

# select only the needed columns
gr = ca_graph.loc[:, ['edge_id', 'from_id', 'to_id', 'd_weighted']]
# , 'from_lon', 'from_lat', 'to_lon', 'to_lat']
gr # this is a basic edge list, with weights.

def edgelist_to_adjlist()
# convert to edge list to adjacency list
adj = {} # empty dictionary for nodes, each
w = {}
for i in range(gr.shape[0]):
    fr_id = gr.loc[i, 'from_id']
    to_id = gr.loc[i, 'to_id']
    if fr_id in adj: # known in dict, so add to list
        adj[fr_id].append(to_id)
        w[fr_id].append(gr.loc[i, 'd_weighted'])
    else: # unknown, introduce to dict and initialize a list
        adj[fr_id] = [to_id] # new list
        w[fr_id] = [gr.loc[i, 'd_weighted']]

len(adj) # number of distinct nodes
