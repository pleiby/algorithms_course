#Uses python3

import sys

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

    return dist[t] 

def Relax(adj, w, dist, u, v):
    # **edge relaxation algorithm** 
    # Relax(u, v):
    if dist[v] > (dist[u] + cost[u][i]): # + w(u,v): # if new shorter way to v via u
        dist[v] = dist[u] + cost[u][i] # update the distance
        prev[v] = u # record/update the predecessor node in path
    return(dist)

def BellmanFord(G_adj):
    # **Bellman–Ford algorithm** for distances from source s in graph G
    # BellmanFord(G, s):
    # {no negative weight cycles in G}
    # for all u∈V:
    #     dist[u] ← ∞
    #     prev[u] ← nil
    # dist[v] will be an upper bound on the actual distance from s to v.

    V = range(len(adj)) # set of nodes, sequentially numbered

    dist = [approxInf for u in V] # initialize dist to completely unknown for all u∈V
    prev = [None for u in V]
    # visited = [False for u in V] # this is represented as dist[u] = infinite

    dist[s] = 0 # zero distance to start node

    # repeat |V|−1 times: # for all u∈V != s
    #     for all (u,v) ∈ E:
    #         Relax(u, v)
    #     while 
    #         at least one dist changes
    
    for i in range(length(V)-1): # repeat |V|−1 times:
        for u in V: # for all u∈V
            for i in range(len(adj[u])): # for all (u,v) ∈ E: Relax(u,v) # relax all _outgoing_ edges from u
                # edge relaxation procedure for an edge (u,v) just checks whether
                #  going from s to v through u improves the current value of dist[v].
                v = adj[u][i] # v in adj[u]
                if dist[v] > (dist[u] + cost[u][i]): # + w(u,v):
                    dist[v] = dist[u] + cost[u][i] # update the distance
                    prev[v] = u # update the predecessor node
                    # ChangePriority(H , v , dist[v]) # rather than priority queue, update dist and scan array for min dist


    u = 0
    return(u)

def negative_cycle(adj, cost):
    """ negative_cycle(adj, cost)

    Return 1 if the graph contains a cycle of negative weight
    and 0 otherwise.
    """

    # **Finding Negative Cycle Algorithm:**
    # Run |V| iterations of Bellman–Ford algorithm, 
    for v in adj:
        u = BellmanFord(adj) # u is the node relaxed by Bellman-Ford
    if (u>0): # non-zero node relaxed after last BellmanFord iteration, cycle found by last iteration
        return(1)
    else: # signal no cycle
        return(0)
    
    # If node v ius relaxed on last iteration, negative cycle exists
    # save node v relaxed on the last iteration
    # v is reachable from a negative cycle
    # To ID the cycler:
    # Start from x ← v (x = Prev(v)??), 
    #    follow the link x ← prev[x] for |V| times # will be definitely on the cycle
    #    Repeat
    #        Save y ← x
    #        go x ← prev[x]
    #         until x = y again
    return 0

# Note:
# if edges are multiplicative weights, for path with maximum product:
#    nlr_ei = (−log(r_ei)) and find shortest path length

def parse_weighted_digraph_input_to_G(inputtext):
    """
    Expect text file/string describing graph consistent with
    Graph Algorithms course standard, for digraph with weights
    Returns:
        adj - adjacency matrix
        cost - edge costs organized in same way as adjacency matrix
    """
    data = list(map(int, inputtext.split()))
    n, m = data[0:2] # count of verts and edges
    data = data[2:]
    # pick off every third element and organize into a m x 3 edge list
    edges = list(zip(zip(data[0:(3 * m):3], data[1:(3 * m):3]), data[2:(3 * m):3]))
    # following assumes that n vertices are number sequentially 1 to n
    adj = [[] for _ in range(n)] # list of n lists, one for each node 
    cost = [[] for _ in range(n)] # organize costs like edge list
    for ((a, b), w) in edges:
        adj[a - 1].append(b - 1)
        cost[a - 1].append(w)

    return (adj, cost)

# Task. Given an directed graph with possibly negative edge weights
#  and with 𝑛 vertices and 𝑚 edges, check whether it contains a cycle
#  of negative weight.
# Input Format: A graph is given in the standard format.
# Constraints: 1 ≤ 𝑛 ≤ 10^3, 0 ≤ 𝑚 ≤ 10^4,
#  edge weights are integers of absolute value at most 103.
# Output Format: Output 1 if the graph contains a cycle of negative weight
#  and 0 otherwise.

if __name__ == '__main__':
    debug = False
    readFromStandardInput = True

    if readFromStandardInput:
        (adj, cost) = parse_weighted_digraph_input_to_G(sys.stdin.read())
    else: # expect a named data structure (list) to read from
        (adj, cost) = parse_weighted_digraph_input_to_G(sample_wdigraph1)

    approxInf = math.inf # establish an impossibly far distance, signal upper bound
    print(negative_cycle(adj, cost))
