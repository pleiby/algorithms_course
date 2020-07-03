#Uses python3

import sys
import queue

#=============================================

sample_bipartite_test_graph1 = """
4 4
1 2
4 1
2 3
3 1
"""

sample_bipartite_test_graph2 = """
5 4
5 2
4 2
3 4
1 4
"""

#=============================================

def parse_input_to_G_adj_list(inputtext):
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
    return adj

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


#==============================================
# Following is a simple algorithm to find out whether a given graph is Birpartite or not using Breadth First Search (BFS).
# 1. Assign RED color to the source vertex (putting into set U).
# 2. Color all the neighbors with BLUE color (putting into set V).
# 3. Color all neighbor’s neighbor with RED color (putting into set U).
# 4. This way, assign color to all vertices such that it satisfies all
#    the constraints of m way coloring problem where m = 2.
# 5. While assigning colors,
#    if we find a neighbor which is colored with same color as current vertex,
#    then the graph cannot be colored with 2 vertices (or graph is not Bipartite)

def BFS_bipartite_labeling(adj,s):
    """Breadth First Search of a graph
    adj is node adjacency list for graph G defined by vertices V and edgelist E,
    s is the starting node from which BFS is conducted
    returns distances to all nodes from s
    """
    # The running time of breadth-first search is O(|E|+|V|)
    # V = range(len(adj)) # sequence of nodes
    # Note!!: this is not entirely general - there is no guarantee that
    #   the graph node list is sequentially numbered from 0 to n-1
    success = -1 # assume not bipartite

    color = [None for u in range(len(adj))] # initialize distance to unk for all u∈V, signal not visited

    color[s] = True # set color of start node (use bool to encode 2 colors)
    Q = queue.Queue(maxsize = len(adj)+1) # initialize a sufficiently large queue
    enqueue(Q,s) # queue containing just s
    while not Q.empty(): # Q !empty
        (Q,u) = dequeue(Q) # ? is there a better way than passing this queue around?
        for v in adj[u]: # all (u,v) ∈ E
            if color[v] == None: # have not explored yet
                Q = enqueue(Q,v)
                color[v] = not color[u] # toggle color & signal node visited
            else: 
                # have explored and colored - adjacent dest nodes cannot match origin node color
                if color[v] == color[u]:
                    return(color, 0) # return colors (so far) and signal failure

    # if made it through two-coloring process with no adjacent nodes of same color
    return(color, 1)

def bipartite(adj):
    """
    For a graph specified by adjacency list `adj`,
    return 1 if the graph is bipartite and 0 otherwise.
    """
    #write your code here
    s = 1 # do we need to call different starting nodes until everything is explored, if not connected?
    (color, success) = BFS_bipartite_labeling(adj,s)
    return success 

# Task. Given an undirected graph with 𝑛 vertices and 𝑚 edges,
#   check whether it is bipartite.
# An undirected graph is called bipartite if its vertices can be split into
#   two independent sets such that each edge of the graph joins two vertices
#   from different sets.
#   I.e., the nodes of G are sepable into disjoint sets U and V such that
#    for every edge (u, v) of G, u is in U and v in V, or vice versa.
#     No edge that connects vertices of same set.
#   An alternative definition: a graph is bipartite if its vertices
#   can be colored with two colors (say, black and white) such that the
#   endpoints of each edge have different colors.

if __name__ == '__main__':
    # Input Format. A graph is given in the standard format.
    # Constraints. 1 ≤ 𝑛 ≤ 10E5, 0 ≤ 𝑚 ≤ 10E5.
    debug = False
    adj = parse_input_to_G_adj_list(sys.stdin.read())
    # adj = parse_input_to_G_adj_list(sample_bipartite_test_graph2)
    print(bipartite(adj)) # return 1 if the graph is bipartite and 0 otherwise
