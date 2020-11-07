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

def distance(adj, cost, s, t):
    #write your code here
    dist_to_t = dijkstra(adj, cost, s, t)
    if dist_to_t < approxInf:
        return dist_to_t
    else:
        return -1
