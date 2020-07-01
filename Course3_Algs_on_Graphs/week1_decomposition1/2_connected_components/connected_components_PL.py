#Uses python3

import sys

def parseGraphInputToAdjlist(stdGraphInput):
    """Read standard input format for Graph (for this set of Coursera courses)
       This is a line with numVerices, numEdges, 
       followed by lines giving Edge list.
       Parse into an Adjacency list, and return that.
    """
    data = list(map(int, input.split())) # read graph description
    n, m = data[0:2] # n vertices, m edges
    data = data[2:]
    # create edge list (m ordered pairs)
    edges = list(zip(data[0:(2 * m):2], data[1:(2 * m):2]))
    # create adjacency list
    # renumber all vertices from 0 to n-1, easing pythonic indexing
    adj = [[] for _ in range(n)]
    for (a, b) in edges:
        adj[a - 1].append(b - 1)
        adj[b - 1].append(a - 1)
    return adj


def previsit(v):
    global CCnum

    CCnum[v] = cc
    return

def postvisit(v):
    return

def explore(v):
    """explore a connected component of a graph starting from node `v`
    Uses globals: graph adjacency list `adj`
    Returns when all direct and subsidiary adjacent nodes have been visited."""
    global visited

    previsit(v)
    visited[v] = True # only explored nodes are marked as visited
    for w in adj[v]: # for nodes w in adj list for v, i.e., for (v,w) in E:
        if not visited[w]: explore(w)
    postvisit(v)

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


def number_of_components(adj):
    """Task. Given an undirected graph specified as an adjacency list,
    compute the number of connected components in it."""

    result = DFS() - 1 # DFS increments cc one too many times?

    #for vert in adj:
    #    continue
    #write your code here
    return result

# Examples graphs as strings
connected_sample1 = """4 2
1 2
3 2
"""

# Task: Given an undirected graph with n vertices and m edges,
#   compute the number of connected components in it.
# Input Format: 
#    A graph is given in the standard format: First line is nverts, nedges. 
#    Subsequent nedges lines are vertex pairs
# Constraints: 1 ≤ n 10^3, 0 ≤ m ≤ 10^3.
# Output Format: Output the number of connected components.
# Time Limits:
#   language  C  C++ Java Python 
#   time(sec) 1  1   1.5  5
# Memory Limit. 512MB

if __name__ == '__main__':
    input = sys.stdin.read()
    
    # file1 = open("connected_sample1.txt", "r")
    # input = file1.read()

    # input = connected_sample1 # use the locally defined text example

    debug = False

    adj = parseGraphInputToAdjlist(input)
    # adj is an adjacency list, i.e. a list over vertices, of lists of adjacencies for each vertex
    if (debug):
        for edge in adj: # or enumerate(ver):
            print(edge, end =" ")
        print()

    # initialize (globals)
    visited = [False for i in range(len(adj))]
    V = range(len(adj))
    cc = 0 # number of connected components
    CCnum = [None for i in range(len(adj))] # conn comp for each vertex


    if (debug): 
        print("Vertices: ", V)
        print("For graph with adjacencies adj =", adj)

    print(number_of_components(adj))
