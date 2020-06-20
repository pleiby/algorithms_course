#Uses python3

import sys

def parseGraphInputToAdjlist(stdGraphInput):
    """Read standard input format for Graph (for this set of Coursera courses)
       This is a line with numVerices, numEdges, follow by lines giving Edge list.
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
    # CCnum(v) = cc
    return

def postvisit(v):
    return

def explore(v):
    previsit(v)
    visited[v] = True
    for w in adj[v]: # i.e., for (v,w) in E:
        if not visited[w]: explore(w)
    postvisit(v)

def DFS(adjV): # adjV is the list of vertices in graph G
    for v in adj:
        visited[v] = False
    cc = 1
    for v in V:
        if not visited[v]:
            explore(v)
        cc = cc + 1
    return


def number_of_components(adj):
    """Task. Given an undirected graph with n vertices and m edges, 
    compute the number of connected components in it."""

    result = 0

    for vert in adj:
        continue
    #write your code here
    return result

# Task: Given an undirected graph with n vertices and m edges, compute the number of connected components in it.
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

    debug = True

    adj = parseGraphInputToAdjlist(input)
    # adj is an adjacency list, i.e. a list over vertices, of lists of adjacencies for each vertex
    if (debug):
        for edge in adj: # or enumerate(ver):
            print(edge, end =" ")
        print()

    # initialize
    visited = [False for i in range(len(adj))]
    V = range(len(adj))


    if (debug): 
        print("Vertices: ", V)
        print("For graph with adjacencies adj =", adj)

    print(number_of_components(adj))
