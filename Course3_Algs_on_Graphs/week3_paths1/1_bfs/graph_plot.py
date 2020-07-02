#Uses python3

import networkx as nx
# see: https://www.python-course.eu/networkx.php
import matplotlib.pyplot as plt

#=============================================

# "In NetworkX, nodes can be any hashable object e.g., a text string,
#    an image, an XML object, another Graph, a customized node object, etc"
# "Four basic graph properties facilitate reporting: G.nodes, G.edges, G.adj and G.degree."
# "...an edge can be associated with any object x using G.add_edge(n1, n2, object=x)."

def build_nx_Graph(adj):
    """
    Specify nodes and edges in an nx Graph one by one,
    from an adjacency list `adj`.
    Return the nx graph object G
    """
    G = nx.Graph()
    for v in range(len(adj)): # ?Fix again, assumes sequential node numbering
        G.add_node(v) # nx.Graphs likes all nodes to be referenced as string labels

    for v in range(len(adj)): # ?Fix again, assumes sequential node numbering
        for u in adj[v]:
            G.add_edge(v, u)

    return(G)

def build_nx_Graph_from_edge_list(E):
    """
    Simple and direct construction of nx undirected Graph object
    from graph edge list `E`
    Return the nx graph object G
    """
    # Aside: Note Other possibilities:
    # DiGraph: DG = nx.DiGraph()
    # DiGraph from weighted edge list: DG.add_weighted_edges_from([(1, 2, 0.5), (3, 1, 0.75)])

    G = nx.Graph()
    G.add_edges_from(E)

def print_nx_Graph(G):
    print("Nodes of graph: ")
    print(G.nodes())
    print("Edges of graph: ")
    print(G.edges())


def plotGraph(adj):
    """
    Plot a graph bases on its adjacency list `adj`
    """
    G = build_nx_Graph(adj)
    nx.draw(G)

plotGraph(adj)


