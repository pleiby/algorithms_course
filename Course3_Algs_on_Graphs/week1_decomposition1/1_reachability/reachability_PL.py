#Uses python3

import sys

def reach(adj, x, y):
    #write your code here
    #print("Hello!")
    #print("Seek to connect (x,y) = (", x, ",", y,")")
    # adj is a list over vertices, of lists of adjacencies for each vertex
    #print(x, end =": ")
    #for edge in adj[x]: # or enumerate(ver):
    #    print(edge, end =" ")
    #print()
    result = 0
    visited[x] = True
    #print("Visited =", visited)
    for edge in adj[x]: # or enumerate(ver):
        if edge == y:
            result = 1
        else:
            if (not visited[edge]):
                result = reach(adj, edge, y)
        if (result==1):
            break
    return result

"""
Starts by reading a graph.
Graphs are given as follows.
The first line contains non-negative integers 𝑛 and 𝑚 
 — the number of vertices and the number of edges respectively. 
The vertices are always numbered from 1 to 𝑛. 
Each of the following 𝑚 lines defines an edge in the format u v where
1 ≤ 𝑢, 𝑣 ≤ 𝑛 are endpoints of the edge. 
If the problem deals with an undirected graph this defines an undirected edge between 𝑢 and 𝑣. 
In case of a directed graph this defines a directed edge from 𝑢 to 𝑣. 
If the problem deals with a weighted graph 
then each edge is given as u v w where 𝑢 and 𝑣 are vertices and 𝑤 is a weight.
"""
if __name__ == '__main__':
    input = sys.stdin.read()
    data = list(map(int, input.split())) # read graph description
    n, m = data[0:2] # n vertices, m edges
    data = data[2:]
    # create edge list (m ordered pairs)
    edges = list(zip(data[0:(2 * m):2], data[1:(2 * m):2]))
    # read the last ordered pair
    x, y = data[2 * m:]
    # create adjacency list
    adj = [[] for _ in range(n)]
    # renumber all vertices from 0 to n-1, easing pythonic indexing
    x, y = x - 1, y - 1
    for (a, b) in edges:
        adj[a - 1].append(b - 1)
        adj[b - 1].append(a - 1)
    visited = [False for i in range(len(adj))]
    # print("For graph with adjacencies adj =", adj)
    print(reach(adj, x, y))
