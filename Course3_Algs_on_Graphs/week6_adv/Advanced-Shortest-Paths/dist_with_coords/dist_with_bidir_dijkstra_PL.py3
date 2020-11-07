#!/usr/bin/python3

import sys
import queue
import math

#=============================================

sample_wdigraph1 = """
4 4
1 2 1
4 1 2
2 3 2
1 3 5
1 3
"""

sample_wdigraph2 = """
5 9
1 2 4
1 3 2
2 3 2
3 2 1
2 4 2
3 5 4
5 4 1
2 5 3
3 4 4
1 5
"""

sample_wdigraph3 = """
3 3
1 2 7
1 3 5
2 3 2
3 2
"""

#=============================================

def make_queue(V):
    """
    Create a queue of paired vertices as a simple list
    """
    H = []
    for i in V:
        H.append(i)
    return(H)

def extract_min(H, ds):
    """
    extract_min(H, ds)

    extract the node from queue `H` with the minimal upper-bound estimate of distance to source s.
    `ds` distance array is also passed, absent a priority queue implementation.
    For this node v, the bound distance ds(v) will be the actual distance d(s,v).
    Return the min dist node, its distance
    (but not the reduced set H).
    """
    minDist = approxInf
    u = None # min vertex unknown
    i = 0
    for v in H:
        if ds[v] <= minDist:
            minDist = ds[v]
            u = v # note that u is unused (instead returned by pop)
            imin = i
        i += 1
    return(H.pop(imin)) # return [u, d]

# In Dijkstra: we generate an SPT (shortest path tree) for a given source `S` as root.
# The algorithm maintains two sets:
#  The set H of "unknown" vertices includes vertices that that have not been fully evalated,
#    and that are not yet included in the shortest-path tree.
#  The other set/region, R = not(H), contains vertices whose distance to root is correctl known,
#    and which are included in the shortest-path tree.

def dijkstra(adj, cost, s, t):
    """
    dijkstra(adj, cost, s, t)

    From weighted, directed graph G represented as adjacency matrix `adj`,
    and (non-negative) edge weights organized similarly in `cost`,
    return the distance or shortest path from source `s` to terminus `t`.
    Return -1 if no path found, or any edge weight is negative.
    """

    V = range(len(adj)) # set of nodes, sequentially numbered
    # Note!!: this is not entirely general - there is no quarantee that
    #   the graph node list is sequentially numbered from 0 to n-1

    # for all u∈V:
    #   dist[u] ← ∞, prev[u] ← nil
    # dist[v] will be an upper bound on the actual distance from s to v.
    dist = [approxInf for u in V] # initialize dist to completely unknown for all u∈V
    prev = [None for u in V]
    # visited = [False for u in V] # this is represented as dist[u] = infinite

    dist[s] = 0 # zero distance to start node

    # H ← MakeQueue(V ) {dist-values as keys} # this is the Unknown region, not(R)
    # the set of unknown (unvisited, or not fully visited) vertices
    H = make_queue(V) #, dist)

    while len(H) > 0: # H, set of unknown vertices is not empty:
        # On each iteration we take a vertex outside of R (in H) with the minimal dist-value,
        #  add it to R, and relax all its outgoing edges.
        u = extract_min(H, dist) # [u, d] = extract_min(H)
        # Lemma: When a node u is selected via ExtractMin, dist[u] = d(S,u), actual minimum distance.
        # First node to be extracted will be the source s (since dist[s]==0)
        # Should we stop early if min node u == t (t is moved to known set R before unknown H is exhausted)?
        for i in range(len(adj[u])): # for all (u,v) ∈ E: Relax(u,v) # relax all _outgoing_ edges from u
            # edge relaxation procedure for an edge (u,v) just checks whether
            #  going from s to v through u improves the current value of dist[v].
            v = adj[u][i] # v in adj[u]
            if dist[v] > (dist[u] + cost[u][i]): # + w(u,v):
                dist[v] = dist[u] + cost[u][i] # update the distance
                prev[v] = u # update the predecessor node
                # ChangePriority(H , v , dist[v]) # rather than priority queue, update dist and scan array for min dist
        if u == t: 
            break # continue untill have processed terminal/target node `t`
    return dist[t] 

def reverse_graphs(adj, cost):
    """
    reverse_graphs(adj, cost)

    reverse a graph, given its adjacency and cost matrix
    Returns reversed adj and costs
    """
    # see https://stackoverflow.com/questions/16158926/how-to-reverse-a-graph-in-linear-time

    radj = [[] for _ in range(len(adj))] # list of n lists, one for each node 
    rcost = [[] for _ in range(len(cost))] # organize costs like edge list

    for u in range(len(adj)):
        # for u,v in G, (v,u) is in G^R
        for i in range(len(adj[u])): # for all (u,v) ∈ E: Relax(u,v) # relax all _outgoing_ edges from u
            v = adj[u][i] # node in ith position of adjacencies to u
            radj[v].append(u) # v in adj[u]
            rcost[v].append(cost[u][i])

    return(radj, rcost)

def Relax(adj, w, dist, u, v):
    # **edge relaxation algorithm** 
    # Relax(u, v):
    if dist[v] > (dist[u] + cost[u][i]): # + w(u,v): # if new shorter way to v via u
        dist[v] = dist[u] + cost[u][i] # update the distance
        prev[v] = u # record/update the predecessor node in path
    return(dist)

# ShortestPath(s, dist, prev, proc, t, distR , prevR , procR )
# distance ← +∞, ubest ← None
# for u in proc+procR:
#     if dist[u]+distR[u]<distance:
#         ubest ← u
#         distance ← dist[u] + distR [u]
# path ← empty
# last ← ubest
# while last ̸= s:
#     path.Append(last)
#     last ← prev[last]
# path ← Reverse(path)
# last ← ubest
# while last ̸= t:
#     last ← prevR[last]
#     path.Append(last)
# return(distance,path)

def bidir_dijsktra(adj, cost, s, t):
    """
    bidir_dijkstra(adj, cost, s, t)

    From weighted, directed graph G represented as adjacency matrix `adj`,
    and (non-negative) edge weights organized similarly in `cost`,
    return the distance or shortest path from source `s` to terminus `t`.
    Return -1 if no path found, or any edge weight is negative.
    """

    Relax(u, v , dist, prev)
    if dist[v]>dist[u]+w(u,v): dist[v] ← dist[u] + w(u, v) prev[v] ← u

    # Process(u, G , dist, prev, proc)
    # for (u,v)∈E(G):
    #     Relax(u, v , dist, prev)
    # proc.Append(u)

    # BidirectionalDijkstra(G, s, t)
    # GR ← ReverseGraph(G)
    # Fill dist,distR with +∞ for each node
    # dist[s] ← 0, distR [t] ← 0
    # Fill prev,prevR with None for each node
    # proc ← empty, procR ← empty
    # do:
    #     v ← ExtractMin(dist)
    #     Process(v, G, dist, prev, proc)
    #     if v in procR:
    #         return ShortestPath(s, dist, prev, proc, t, ... )
    #     vR ← ExtractMin(distR)
    #     repeat symmetrically for vR as for v
    # while True

    # first construct the reverse graph G^R
    radj, rcost = reverse_graphs(adj, cost)
    if debug:
        print(adj, cost, s, t)
        print(radj, rcost, s, t)

    return(dist)

def distance(adj, cost, s, t):
    """
    distance(adj, cost, s, t)

    given adjacency list `adj` and similarly organized costs list `cost`
    return distance from s to t
    """
    dist_to_t = dijkstra(adj, cost, s, t)
    if dist_to_t < approxInf:
        return dist_to_t
    else:
        return -1

def parse_weighted_digraph_input_to_G_s_and_t(inputtext):
    """
    Expect text file/string describing graph consistent with
    Graph Algorithms course standard,
    followed by a line with start and target node numbers.
    Returns:
        adj - adjacency matrix
        cost - edge costs organized in same way as adjacency matrix
        s - source node (origin)
        t - terminal node (destination)
    """
    data = list(map(int, inputtext.split()))
    n, m = data[0:2] # count of verts and edges
    data = data[2:]
    # pick off every third element and organize into a m x 3 edge list
    edges = list(zip(zip(data[0:(3 * m):3], data[1:(3 * m):3]), data[2:(3 * m):3]))
    data = data[3 * m:] # last line is (s, t)
    # following assumes that n vertices are number sequentially 1 to n
    adj = [[] for _ in range(n)] # list of n lists, one for each node 
    cost = [[] for _ in range(n)] # organize costs like edge list
    for ((a, b), w) in edges:
        adj[a - 1].append(b - 1)
        cost[a - 1].append(w)
    s, t = data[0] - 1, data[1] - 1
    return (adj, cost, s, t)

if __name__ == '__main__':
    debug = True
    readFromStandardInput = False

    if readFromStandardInput: # expect weighted digraph followed by s and t
        (adj, cost, s, t) = parse_weighted_digraph_input_to_G_s_and_t(sys.stdin.read())
    else: # expect a named data structure (list) to read from
        (adj, cost, s, t) = parse_weighted_digraph_input_to_G_s_and_t(sample_wdigraph1)

    approxInf = math.inf # establish an impossibly far distance, signal upper bound

    print(distance(adj, cost, s, t))
