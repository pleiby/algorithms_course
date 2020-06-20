Notes on Graph Algorithms - from UC San Diego Course
====================================================

Week 1 Graphs and Decomposition
-------------------------------
- The generality and wide applicability of graphs is surprising, and intriguing
- a graph can describe a general system of discrete states (vertices/nodes) and allowble sequential connections or transitions (edges/arcs)
    - e.g. for a puzzle or a mechanical system, they can describe alternative, feasible, discret configurations, and moves between them.

### Representation
#### Summary
- Different operations are faster in different representations.
- And depending on the density (degree of connectedness, number of edges per vertex)

Operation: | Is Edge? | List Edge | List Nbrs.
-----------|----------|-----------|-------------
Adj. Matrix| Θ(1)     | Θ(|V|^2)  | Θ(|V|)
Edge List  | Θ(|E|)   | Θ(|E|)    | Θ(|E|)
Adj. List  | Θ(deg)   | Θ(|E|)    | Θ(deg)

- "deg" = degree of connectedness
- For many problems, want adjacency list.
- sparse graphs |E| ~ |V|, deg ~ 1
- dense graphs |E| ~ |V|^2

### Exploration Algorithms
- Three parts
    - list of visited vertices/nodes
    - list of vertices that still require more processing (e.g. have edges that have not been explored) - not explicity
    - search order - embodied in a function "explore()", that is called on unvisited vertices
        - which direction to start in (and proceed in after visiting/processing a vertex)
            - for Depth First Search (DFS), proceed via any edge of this vertex to a successive vertex
                - Note: For this, and adjacency list is efficient
        - what to do when reach node with particular quality (e.g. has been already visited)
            - for Depth First Search (DFS), backtrack to first vertex with edge to unexplored node
        - what to do when reach "dead end", last node in connected sequence (has no unvisited neighbors)
            - Backtrack one level
            - Backtracking can be achieved by popping recursive stack of calls to "explore()"
    - DFG(G)
        - calls explore(v)
- Graph algorithm runtime
    - consider how many times something (some computational work) must be done for each vertex (e.g. a call explore())), some function of number of vertices |V|
    - consider how many times each edge from a vertex must be processed (traversed) in some way, implies some function of total number of edges |E|
    - DFS is O(|V|) + O(|E|) = O(|V| + |E|) (linear time)

### Connected Components in Graph
- Thm: The vertices of a (_undirected_) graph G can be partitioned into _Connected Components_ so that v is reachable from w if and only if they are in the same connected component
- reachability R is an _equivalence relation_
    - reflexive: v R v
    - symmetric: v R w => w R v
    - transitive: u R v and v R w => u R w

- Modifying DFS with Previsit and Postvisit operations, to maintain and return useful info about search and its order
    - explore(v )
        - visited(v) ← true previsit(v)
        - for (v,w) ∈ E:
            - if not visited(w): explore(w)
        - postvisit(v)
- previsit() and postvisit() augment the explore() function, to keep track of exploration results
- E.g. to record sequence of visits, use a clock
    - Clock ticks (increments) at each pre-/post- visit.
    - Include data pre(v) and post(v) for each node v
    - Record previsit and postvisit times (clock count) for each v, as data pre(v) and post(v)

Directed Graphs
----------------
- imposes order on the structure, or dependency relations
- Q; can we linearly order notes to indicate depenencies?
    - A: Not always, e.g. in any case where graph has cycles
    - If G contains a cycle, it cannot be linearly ordered. (prove by contraction)
- TO be linearly orderable, must be a **Directed Acyclic Graph (DAG)**
    - Theorem: If G contains a cycle, it cannot be linearly ordered.
    - Is converse true? Yes
    - Theorem: Any DAG can be linearly ordered (topologically sorted)
        - proof by an algorithm for doing that
- Algorithm to Topologically Sort DAGs
    - conside which vertex comes last: no outgoing edges
        - source: has no incoming edges
        - sink: has no outgoing edges
    - to order, start by finding a sink, put at end of the order
        - remove sink and repeat (recording order of removed nodes) until no nodes left
    - to find a sink
        - start anywhere, and just keep following path down until either:
            - cannot extend, dead end (found a sink, vertex with now outgoing edges)
            - repeat a vertex (have a cycle: not a DAG)
- LinearOrder(G)
    - while G non-empty:
        - Follow a path until cannot extend 
        - Find sink v
        - Put v at end of order [put at front of list of removed nodes?] 
        - Remove v from G
    - Runtime
        - O(|V|) paths, (i.e. starting at the beginning, walk a path for each sink found)
        - Each takes O(|V|) time. (may have to walk down all the nodes)
        - Runtime O(|V|^2). (not great)
    - Improved algorithm
        - Rather than Retrace same path from start vertex every time, 
            - Instead only back up as far as necessary (one step from sink).
            - Then go forward again
        - This is same sequence as DFS
            - whenever finish post-visit block at a vertex, "put it at the end of our ordering" [?]
            - sorting vertices based on the postorder!
- Topological Sort(G)
    - DFS(G)
    - sort vertices by reverse post-visit order
- correctness:
    - Theorem:
        - If G is a DAG, with an edge u to v, (u comes before v, so) post(u) > post(v)

### Strongly Connected Graphs in Directed Graphs (Digraphs)
- in Undirected Graphs, have connected components, and unconnected
- for Directed Graphs, vertices may be "connected", but that does not tell you which vertices can be reached from which others
- alternative notions of connectivity in Digraphs
    - connected by edges in any direction
        - (can't be separated from one another, but doesn't say anything about reachability)
    - one vertex reachable from the other (awkward)
    - both vertices reachable from each other (stronger connectivity notion, mutual reachability)
- this last is the sense in which Digraphs are said to be (strongly) connected
- **Strongly Connected Components** in a Digraph
    - Definition:
        - Two vertices v,w in a directed graph are **connected** if you can reach v from w and also can reach w from v.
    - Theorem: Digraphs partitionable into SCCs
        - A directed graph can (always) be partitioned into **strongly connected components (SCC)** where two vertices are connected if and only if they are in the same component.
            - i.e. two vertices are in the same strongly connected) component iff they are connected to each other (mutually reachable)
    - Proof: 
        - follows connectivity partitioning proof for undirected graph
        - need to show that connectivity is an equivalence relation (reflexive, symmetric and transitive)
- metagraph: shows connections among graph SCCs
    - Theorem: The metagraph of a graph G is always a DAG.
        - Proof: by contraction (if a cycle, then components (SCCs) involved are actually connected, i.e. part of the same SCC, by transitivity of connectivity)
- Summary
    - Can partition vertices into strongly connected components.
    - Metagraph describes how strongly connected components connect to each other.
    - Metagraph always a DAG.

### Computing Strongly Connected Components
- goal: Effciently compute SCC(G), the strongly connected components of a directed graph.
- Easy Algorithm
    - EasySCC(G)
        - for each vertex v:
            - run `explore(v)` to determine vertices reachable from v 
        - for each vertex v:
            - find the u reachable from v that can also reach v
            - these are the SCCs
    - A bit slow: Runtime O(|V|^2 + |V||E|). Want faster.
    - key concept: look for "sink Strongly Connected Component" - an SCC with no exit links in the Metagraph
        - that will avoid exploring other nodes in the `explore` step that are not strongly connected (that are part of another SCC)
        - That is, if v is in a sink SCC, `explore(v)` finds vertices "reachable" from v, i.e. finds exactly the SCC of v (a good start)
        - How to find a vertex v in a sink SCC?
    - Useful Theorem
        - If C and C′ are two strongly connected components with an edge _from_ some vertex of C _to_ some vertex of C′, then largest post in C bigger than largest post in C′.
            - Proof: when explore C` cannot reach C, because that would imply a cycle in the metagraph (and the Metagraph(G) is a DAG for all Digraphs G), so contradioction.
    - Conclusion from this Thm: 
        - The vertex with the largest postorder number is in a source SCC
        - Problem: We wanted a sink component.
    - Trick: Explore a *Reverse Graph*: 
        - G^R = G with all edges reversed in direction (easy to computee)
        - Assertion/obvious: 
            - GR and G have same SCCs. (SCC(G^R) = SCC(G))
            - Source components of GR are sink components of G.
        - Can find vertex in a sink components of G, and sink components of G, by
            - running DFS on G^R
            - taking the vertex with the largest postorder as the source (sink of G)
                - "The vertex with largest postorder in Reverse(G) is in a sink SCC of G."
- Algoritm for finding SCCs
    - SCCs_v0(G) # first draft
        - run DFS(G^R)
        - let v have largest post number 
        - run Explore(v)
        - ID vertices found as first SCC # (because v was in a sink component, with no exiting edges)
        - Remove sink SCC found from G and
        - repeat (until empty of vertices)
    - Improvement: Don't need to rerun DFS on G^R each time
        - recall on Digraph
            - upstream vertex has the largest post-order
            - for reverse graph the downstream one has the largest postorder
        - Vertex with largest remaining post number (from DFS on Reverse(G)) comes from sink component.
            - so after removing sink component found, look to remaining vertex with highest postorder, 
                - which will be a sink node of the new graph
                - i.e., which will be in a sink component of new graph
    - SCCs(G)
        - Run DFS(G^R)
            - record reverse postorder
        - for v ∈ V in reverse postorder:
            - if not visited(v): 
                - Explore(v )
                - mark visited vertices as new SCC
    - Runtime SSCs(G)
        - Essentially two DFSs: DFS on GR and then DFS on G. 
        Runtime O(|V|+|E|).
- Decompostion algs tell to find ways to get from one vertex to another
- Next time: efficient/shortest paths

- Questions: doesnt Kane mispeak when he says, to order Digraph, find sink nodes and sequentially add each to the "end" of the ordered list of nodes?
- Questions: "09_graph_decomposition_9_computing-sccs.pdf, Slide 22: how is node numbered 11/18 reached in DFS(G^R)?
    - A: This involves a restart of the DFS at a new start node (11/18), since cannot get from 1 to 11

Shortest Paths on Graphs (Week 3)
---------------------------------

### Paths and Distances
- Motivation: interested in Minimum switches or minimum segments between 2 vertices, on Unirected-graphs or Digraphs
- **Most Direct Route** problem: minimum number of segments (edges) from vertex A to B
    - e.g. flights from Hamburg to Moscow
- Define L(P) Length of Path = Number of edges in path
    - (equivalent to equal-weighted arcs)
    - applies to Directed and Undirected Graphs
- Define Distance(A,B) = Length of _Shortest_ Path(A,B)
- Insight: finding shortest path from A to B "is not simpler than finding shortest paths from A to all other nodes in the graph, or at least almost [?] all the other nodes in the graph."
- Idea of Distance Layers
    - based on a chosen starting vertex S
    - each vertex v is assigned to a number layer based on its distance (L(shortest path) from S
        - implies certain conclusions about the existence (or non-existence) of edges from one node to another
        - different conclusions depending on whether graphs is directed
        - "the general property is that there cannot be any edge from a layer to another layer which is farther from S by at least two."
            - there can be an edge from a layer to the next layer; 
            - There can be an edge within the layer's edges;
            - There can be an edge from a layer to any of the previous layers.
    - use Distance layers to compute distance from starting node S:
        - discuss an efficient algorithm to traverse [from S] the graph layer by layer,
            - so that in the end every node v is assigned to some layer, and
            - we know the distance (from S) to this node v as the number of the layer to which it was assigned.
### BFS
- Idea of Breadth First Search (BFS)
    - BFS is an efficient alg to find distances from origin node S to all other nodes in graph

- **Breadth-first search Algorithm**
    - BFS(G,S) # graph G and origin node S
        - for all u∈V: 
            - dist[u] ← ∞ # initialize all to undiscovered, and potentially unreachable
        - dist[S] ← 0 # init distance to self. "dist" is array or map/dict
        - Q ← {S} {queue containing just S} # all discovered, unprocessed nodes are in Q
        - while Q is not empty:
            - u ← Dequeue(Q) # process as FIFO queue vs LIFO stack, which would be DFS
            - for all (u,v) ∈ E: # use adjacency list E(u)
                - if dist[v] = ∞: # not visited/discovered
                    - Enqueue(Q, v) 
                    - dist[v] ← dist[u] + 1 # changing dist from inf assures it will note be put in Q again

- Note: inf ~ Num(nodes)+1 or Num(edges)+1
- Running time of BFS - Lemma
    - The running time of breadth-first search is O(|E|+|V|).
    - Proof
        - Each vertex is enqueued at most once (then marked as visited/dist < inf)
        - Each edge is examined either once (for directed graphs) or twice (for undirected graphs)

### Proof of Correctness

Question: Discuss: BFS step u ← Dequeue(Q) # process as FIFO queue vs LIFO stack, which would be DFS
        
- Define **Reachability**
    - Node u is reachable from node S if there is a path from S to u
- *Reachability Lemma*
    - Reachable nodes are discovered at some point in BFS(G,S), so they get a finite distance estimate from the source S. 
    - Unreachable nodes are not discovered at any point, and the distance to them stays infinite.
    - Proof: ...
- *Order Lemma*
    - By the time node u at distance d from S is dequeued, all the nodes at distance at most d have already been discovered (enqueued).
    - Proof: ...
- *Correct distances Lemma*
    - When node u is discovered (enqueued), dist[u] is assigned exactly d(S, u).
    - Proof: by induction
    
### Shortest Path Tree
- Shortest path tree usedful, 
    - contains all the distances (and shortest paths) from the origin to other reachable nodes
    - since traversing edges changes distance by exactly 1, so if distance D, can reach D in exactly D steps
- Shortest-path tree Lemma 
    - Shortest-path tree is indeed a tree, i.e. it doesn't contain cycles (it is a connected component by construction).

- to construct shortest-path tree, add concept (attribute) of previous level node to each node
- *Algorithm: Constructing shortest-path tree by BFS*
    - BFS(G,S)
        - for all u∈V: # initialize
            - dist[u] ← ∞, prev[u] ← nil
        - dist[S] ← 0
        - Q ← {S} {queue containing just S} 
        - while Q is not empty:
            - u ← Dequeue(Q)
            - for all (u,v)∈E:
                - if dist[v] = ∞:
                    - Enqueue(Q , v )
                    - dist[v] ← dist[u] + 1, prev[v] ← u
- *Algorithm Reconstructing Shortest Path*
    - ReconstructPath(S, u, prev) # from origin node S to node u, given `prev` 
        - result ← empty
        - while u != S: # start with destination and work back
            - result.append(u)
            - u ← prev[u]
        - return Reverse(result)
- Conclusion
    - Can find the minimum number of segments to get from one node to another
    - Can reconstruct the optimal path
    - Can build the tree of shortest paths from one origin
    - Works in O(|E|+|V|)  

### Bipartite Graph Problem
- [Check whether a given graph is Bipartite or not](https://www.geeksforgeeks.org/bipartite-graph/)
