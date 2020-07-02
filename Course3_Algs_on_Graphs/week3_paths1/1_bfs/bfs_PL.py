#Uses python3

import sys
import queue

#=============================================

sample_undigraph0 = """
4 5
2 1
4 3
1 4
2 4
3 2
1 3
"""

sample_undigraph1 = """
4 4
1 2
4 1
2 3
3 1
2 4
"""

sample_undigraph2 = """
5 4
5 2
1 3
3 4
1 4
3 5
"""

sample_digraph1 = """
5 8
4 3
1 2
3 1
3 4
2 5
5 1
5 4
5 3
"""

sample_digraph2 = """
5 1
4 3
"""

sample_wdigraph1 = """
3 3
2 3 9
1 3 5
1 2 -2
"""

#=============================================

def DFS(): # adj is the list of vertices in graph G
    """DFS is Depth First Search of (undirected) graph.
    adj is adjacency list for graph. 
    Returns number of distinct connected components."""

    global cc
    global visited

    for v in range(len(adj)):
        visited[v] = False
    cc = 1

    for v in range(len(adj)):
        if not visited[v]:
            explore(v)
            # increment connected component count after each return from explore()
            cc = cc + 1 # only increment for each unvisited node explored here
    return cc

#=============================================

import networkx as nx
# see: https://www.python-course.eu/networkx.php
import matplotlib.pyplot as plt


def build_nx_Graph(adj):
    G = nx.Graph()
    for v in range(len(adj)): # ?Fix again, assumes sequential node numbering
        G.add_node(str(v)) # nx.Graphs likes all nodes to be referenced as string labels

    for v in range(len(adj)): # ?Fix again, assumes sequential node numbering
        for u in adj[v]:
            G.add_edge(str(v),str(u))

    return(G)

def print_nx_Graph(G):
    print("Nodes of graph: ")
    print(G.nodes())
    print("Edges of graph: ")
    print(G.edges())


def plotGraph(adj):
    G = build_nx_Graph(adj)
    nx.draw(G)

# plotGraph(adj)

#=============================================

# Queue in Python can be implemented by the following ways <https://www.geeksforgeeks.org/queue-in-python/>:
# list:
#   q = [] # initialize 
#   with append(), pop(). But this is slow
# collections.deque 
#   from collections import deque # double-ended queue
#   q = deque() # Initializing a queue
#   q.append(x) # Adding elements to a queue
#   print(q.popleft()) # Removing elements from a queue 
# queue.Queue
#   from queue import Queue
#   q = Queue(maxsize = 3) # Initializing a queue and setting maxsize
#   queue.maxsize # size limit for queue
#   queue.empty() # Return True if the queue is empty, False otherwise.
#   queue.full() # Return True if there are maxsize items in the queue. If the queue was initialized with maxsize=0 (the default), then full() never returns True.
#   queue.qsize() # Return the number of items in the queue. If no free slot is immediately available, raise QueueFull.
#   queue.put_nowait(item) # Put an item into the queue without blocking.
#   queue.get_nowait() # Return an item if one is immediately available, else raise QueueEmpty.
#   queue.put(item) # Put an item into the queue. If the queue is full, wait until a free slot is available before adding the item.
#   queue.get() # Remove and return an item from the queue. If queue is empty, wait until an item is available.

def enqueue(Q, x):
    """enqueue x on the queue Q"""
    # Q.append(x)
    Q.put_nowait(x)
    if debug: 
        print("enqueue", x, ":", end=" ")
        show_queue(Q)
    return Q


def dequeue(Q):
    """dequeue x from the queue Q"""
    # x = Q.pop(0) # default is to pop from end (LIFO stack), param 0 indicates FIFO queue
    x = Q.get_nowait() # default is to pop from end (LIFO stack), param 0 indicates FIFO queue
    if debug: 
        print("dequeue   :", end=" ")
        show_queue(Q)
    return(Q, x)

def show_queue(Q):
    """display the entire queue in one printed line
    """
    print("(Size of the queue:", Q.qsize(), ")", end=" ")
    for n in list(Q.queue):
        print(n, end=" ")
    print()

#=============================================

def BFS(adj,s):
    """Breadth First Search of a graph
    adj is node adjacency list for graph G defined by vertices V and edgelist E,
    s is the starting node frome which distances are measured
    returns distances to all nodes from s
    """
    # The running time of breadth-first search is O(|E|+|V|)
    V = range(len(adj)) # sequence of nodes
    # Note!!: this is not entirely general - there is no quarantee that
    #   the graph node list is sequentially numbered from 0 to n-1

    approxInf = 2*len(V) # establish an impossibly far distance, signal not visited

    dist = [approxInf for u in V] # initialize distance to unk for all u∈V
    dist[s] = 0 # zero distance to start node
    Q = queue.Queue(maxsize = len(adj)+1) # initialize a sufficiently large queue
    enqueue(Q,s) # queue containing just s
    while not Q.empty(): # Q !empty
        (Q,u) = dequeue(Q) # ? is there a better way than passing this queue around?
        for v in adj[u]: # all (u,v) ∈ E
            if dist[v] == approxInf: # have not explored yet
                Q = enqueue(Q,v)
                dist[v] = dist[u] + 1 # increment distance, & signal node visited
    
    return(dist)

def distance(adj, s, t):
    #write your code here
    dist = BFS(adj,s)
    if dist[t] < (len(adj)+1): # warning: need to check for t outside range of nodes in dist
        return dist[t]
    else:
        return -1




def parse_input_to_G_s_and_t(inputtext):
    """
    Expect text file/string describing graph consistent with
    Graph Algorithms course standard,
    followed by a line with start and target node numbers
    """
    data = list(map(int, inputtext.split()))
    n, m = data[0:2]
    data = data[2:]
    edges = list(zip(data[0:(2 * m):2], data[1:(2 * m):2]))
    adj = [[] for _ in range(n)]
    for (a, b) in edges:
        adj[a - 1].append(b - 1)
        adj[b - 1].append(a - 1)
    s, t = data[2 * m] - 1, data[2 * m + 1] - 1
    return (adj, s, t)

# Task. Given an undirected graph with 𝑛 vertices and 𝑚 edges
# and two vertices 𝑢 and 𝑣,
# compute the length of a shortest path between 𝑢 and 𝑣
# (that is, the minimum number of edges in a path from 𝑢 to 𝑣).

if __name__ == '__main__':
    debug = True
    (adj, s, t) = parse_input_to_G_s_and_t(sample_undigraph2)
    # (adj, s, t) = parse_input_to_G_s_and_t(sys.stdin.read())
    print(distance(adj, s, t))

