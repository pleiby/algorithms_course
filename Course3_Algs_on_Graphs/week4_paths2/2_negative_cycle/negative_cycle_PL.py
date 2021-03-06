#Uses python3

import sys, math
import logging


def BellmanFord_iter(adj, cost, dist, prev):
    # Single iteration through the graph of BellmanFord (relaxing all node edges)
    #     for all (u,v) ∈ E:
    #         Relax(u, v)
    #     while 
    #         at least one dist changes
    
    relaxed = False # true if a node relaxed
    V = range(len(adj)) # set of nodes, sequentially numbered
    for u in V: # for all u∈V
        for i in range(len(adj[u])): # for all (u,v) ∈ E: Relax(u,v) # relax all _outgoing_ edges from u
            # edge relaxation procedure for an edge (u,v) just checks whether
            #  going from s to v through u improves the current value of dist[v].
            v = adj[u][i] # v in adj[u]
            if dist[v] > (dist[u] + cost[u][i]): # + w(u,v): (found a shorter path to v thru u)
                relaxed = True
                dist[v] = dist[u] + cost[u][i] # update the distance
                prev[v] = u # update the predecessor node
                # ChangePriority(H , v , dist[v]) # rather than priority queue, update dist and scan array for min dist

    return relaxed

def BellmanFord(adj, cost, s=0, negcycle_test=False): # default start node
    """
    determine distances from source node using BellmanFord, which allows negative edge weights.

    Args:
        adj (list): node adjacency matrix.
        cost (list): edge weights organized like adj
        s (ind): number of source node for distance calc
        negcycle_test (bool): Truen if extra BF iter wanted to test for neg cycles


    Returns:
        (bool): True if some node was relaxed in last iteration
        (list): updated distances from source node
    """
    # **Bellman–Ford algorithm** for distances from source s in graph G
    # BellmanFord(G, s):
    # {no negative weight cycles in G}
    # for all u∈V:
    #     dist[u] ← ∞
    #     prev[u] ← nil
    # dist[v] will be an upper bound on the actual distance from s to v.

    V = range(len(adj)) # set of nodes, sequentially numbered

    dist = [math.inf for u in V] # initialize dist to completely unknown for all u∈V
    prev = [None for u in V]
    # visited = [False for u in V] # this is represented as dist[u] = infinite

    dist[s] = 0 # zero distance to start node

    # repeat |V|−1 times: # for all u∈V != s
    #     for all (u,v) ∈ E:
    #         Relax(u, v)
    #     while 
    #         at least one dist changes
    
    if negcycle_test:
        iters = len(V)
    else:
        iters = len(V)-1
    
    for x in range(iters): # repeat |V|−1 times, |V| if neg-cycle check: (can stop early if no relaxations)
        # do one iteration of BellmanFord
        relaxed = BellmanFord_iter(adj, cost, dist, prev)
        # logging.info("BF Iter :", x, relaxed, dist) # not clear that logging can display non-string objects
        if debug:
            print("BF Iter :", x, relaxed, dist, "\n")
        if (not relaxed): # no node relaxed this iteration
            break
    return(relaxed, dist) # return number of last relaxed node (-1 if none), and distances


def negative_cycle(adj, cost):
    """
    bool function testing for negative weight cycle in graph

    Args:
        adj (list): node adjacency matrix.
        cost (list): edge weights organized like adj.

    Returns:
        int: 1 if graph contains cycle of negative weight; 0 otherwise.
    """

    # **Finding Negative Cycle Algorithm:**
    # Run |V| iterations of Bellman–Ford algorithm, 
    # First call `BellmanFord()``, which is |V|-1 BellmanFord iterations,
    #  gets distances providing no negative cycles
    # negcycle_test=True does one extra iteration and see if any node is "relaxed" 
    did_a_relaxation, d = BellmanFord(adj, cost, negcycle_test=True) # u is the last node relaxed by Bellman-Ford
   
    if (did_a_relaxation): # non-negative node relaxed after last BellmanFord iteration, cycle found by last iteration
        return(1) # required return by the problem set
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

# Note:
# if edges are multiplicative weights, for path with maximum product:
#    nlr_ei = (−log(r_ei)) and find shortest path length

def parse_weighted_digraph_input_to_G(inputtext):
    """
    parse text file/string describing weighted digraph consistent with
    Graph Algorithms course standard

    Returns:
        adj (list): adjacency matrix
        cost (list): edge costs organized same as adjacency matrix
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

# note that the graph parsing fn will ignore trailing line with source and dest note s and t
sample_wdigraph0 = """
4 4
1 2 -5
4 1 2
2 3 2
3 1 1
"""

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
4 4
1 2 -1
4 1 2
2 3 2
3 1 1
"""


sample_wdigraph4 = """
5 9
1 2 4
1 3 2
2 3 2
3 2 1
2 4 2
3 5 -5
5 4 1
2 5 3
4 1 1
1 5
""" # has neg cycle

sample_wdigraph4a = """
5 9
1 2 4
1 3 2
2 3 2
3 2 1
2 4 2
3 5 -4
5 4 1
2 5 3
4 1 1
1 5
""" # zero distance but not neg cycle

# Task. Given an directed graph with possibly negative edge weights
#  and with 𝑛 vertices and 𝑚 edges, check whether it contains a cycle
#  of negative weight.
# Input Format: A graph is given in the standard format.
# Constraints: 1 ≤ 𝑛 ≤ 10^3, 0 ≤ 𝑚 ≤ 10^4,
#  edge weights are integers of absolute value at most 103.
# Output Format: Output 1 if the graph contains a cycle of negative weight
#  and 0 otherwise.

def main():
    # logging.basicConfig(filename='myapp.log', level=logging.INFO) # an alternative to debug print statements
    # logging.info('Started')

    debug = True
    readFromStandardInput = False

    if readFromStandardInput:
        (adj, cost) = parse_weighted_digraph_input_to_G(sys.stdin.read())
    else: # expect a named data structure (list) to read from
        (adj, cost) = parse_weighted_digraph_input_to_G(sample_wdigraph4a)

    print(negative_cycle(adj, cost))
    # logging.info('Finished')


if __name__ == '__main__':
    main()

