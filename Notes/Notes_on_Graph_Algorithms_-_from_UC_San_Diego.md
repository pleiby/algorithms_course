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

## Week 4 Shortest Path Dijkstra and Bellman-Ford Algorithms

Week 4 Paths in Graphs: Fastest Route
-------------------------------------
- Concept: **Optimal substructure for optimal path**
    - "Any subpath of an optimal path is also optimal."
    - Corollary: The distance along the shortest path is the sum of the distances along any partiction of the shortest path into subsections.
        - If S → . . . → u → t is a shortest path from S to t, then
            - d(S, t) = d(S, u) + w(u, t)
    - This optimal substructure is the foundation of Bellman Dynamic Programming (for optimal control theory)
        - breaks a sequential optimization problem into a sequence of simpler subproblems
    - 'Bellman equation' usually refers to the dynamic programming equation associated with discrete-time optimization problems, but it is more generally applicable to any discrete decision stages (such as choosing paths through graphs)
    - "Principle of Optimality: An optimal policy has the property that whatever the initial state and initial decision are, the remaining decisions must constitute an optimal policy with regard to the state resulting from the first decision. (See Bellman, 1957, Chap. III.3.)" <https://en.wikipedia.org/wiki/Bellman_equation#Bellman's_principle_of_optimality>
    - a problem that can be broken apart like this is said to have optimal substructure
    - rewrite the problem as a recursive definition of the value function:
        - $$V(x_{0}) = \max_{a_{0}} \{ F(x_{0},a_{0})+ \beta V(x_{1}) \}$$
        - subject to the constraints on feasibility of choice $a_0$ and the transition equation describing how choice $a_i$ alters state $x_i$ to $x_{i+1}$: 
            $$a_{0} \in \Gamma (x_{0})$$
            $$x_{1} = T(x_{0},a_{0})$$

- What is the interpretation of the name "edge relaxation?" How is the edge being relaxed?
    - https://stackoverflow.com/questions/10746727/why-do-we-call-it-relaxing-an-edge#:~:text=In%20general%20mathematically%2C%20relaxation%20is,reducing%20the%20number%20of%20constraints.

    - for estimated distance "dist[v]" from source S to v:
    - _provided_ the estimate dist[v] is an _upper bound_ on the actual distance d(S,v)

        **Relax((u,v) ∈ E)**
        if dist[v] > dist[u] + w(u,v): 
            dist[v] ← dist[u] + w(u, v)
            prev[v] ← u   # record closest know predecessor node as `prev`
    
    - Intuition: 
        - edge (u,v) is "relaxed" in the context of distances from source node S
        - edge (u,v) is relaxed in the sense that we are looking for a shortcut to v (from S) via u
            - "relax edge (u,v)" means "attempt to shorten v's distance to s using edge (u,v)"
            - in OR, relaxation means relaxing a constraint
            - Naive and in Dijkstra relaxation means relaxation of the constraint x_uv = 0?
                - removing (u,v) from set of edges excluded from the shortest path??
                - allowing more non-zero variables

        **NaiveDistance(G , S)**
        for all u∈V: 
            dist[u] ← ∞ 
            prev[u] ← nil
        dist[S] ← 0 
        do:
            relax all the edges # in what order?
        while 
            at least one dist changes

    - notation:
        - dist[v] is our estimate of distance from s to v
        - d(s,v) is the actual (shortest-path) distance from s to v

#### Dijkstra's Algorithm
- Main ideas of Dijkstra’s Algorithm
    - maintain a set R of vertices for which dist is already set correctly ("known region").
        - for w in R: dist[w] = d(s,w)
    - Best to explore (relax edges from) that set of known-distance nodes, 
        - since any node w outside R with new minimum estimate dist[w] will also be its shortest path (correct distance d(s,w))
    - On each iteration we 
        - take the vertex outside of R with the minimal estimated dist-value, 
        - add it to R, and 
        - relax all its outgoing edges [if digraph, o.w. all its edges?].
    - Start with including only source vertex S in R (since dist[S] == 0)

- Dijkstra Algorithm Pseudocode

    **Dijkstra(G, S)**
    for all u∈V:
        dist[u] ← ∞, 
        prev[u] ← nil
    dist[S] ← 0
    H ← MakeQueue(V ) {dist-values as keys} # this is the Unknown region, H = not(R)
    while H is not empty:
        u ← ExtractMin(H)
        # Lemma: When a node u is selected via ExtractMin, dist[u] = d(S,u).
        for all (u,v)∈ E: # relax _outgoing_ edges from u
            if dist[v] > dist[u]+w(u,v):
                dist[v] ← dist[u] + w(u, v)
                prev[v] ← u
                ChangePriority(H , v , dist [v])

#### Dijstra Running Time
- Total running time:
    - T (MakeQueue) 
    - + |V | · T (ExtractMin) 
    - + |E | · T (ChangePriority)
- Actual running time will depend on the efficiency of these operations
    - Priority queue implementations:
        - array:
            - O(|V|+|V|2 +|E|)=O(|V|2)
        - binary heap:
            - O(|V|+|V|log|V|+|E|log|V|) = O((|V|+|E|)log|V|)
    - Conclude DA works in time quadratic in V, or |E| log |V|, depending on implementation
        - choice of implementation depends on density, i.e. E/V (bi-heap better for sparse graph)
- Visualization of DA:
    - [Visualization of Dijkstra'a algorithm](http://www.cs.usfca.edu/~galles/visualization/Dijkstra.html) by David Galles
- References
    - Sections 4.3 and 4.4 in [DPV]
        - [DPV] Sanjoy Dasgupta, Christos Papadimitriou, and Umesh Vazirani. Algorithms (1st Edition). McGraw-Hill Higher Education. 2008.
- Algorithm works for case of non-negative edge weights

#### Variant Graphs and Shortest Distance Problems
- if shortest path is the _product_ of edge weights, can convert that to the sum by taking the log of edge weights
    - relying on $\prod_i x_i = \sum_i ln(x_i)$, provided $x_i > 0$
    - (note $ln(x) < 0$ for $x<1.0$, so may have negative weights)
- if the goal is to maximize path length, rather than min, can take negative of edge weights
    - again we have the issue of negative weights
- Dijsktra algorithm does not work with negative weights, 
    - in part because edge relaxation does not work, i.e. cannot be certain that paths made shorter by relaxation necessarily dominate partial paths that are longer.
    - Dijkstra algorithm relies on fact that shortest path from s to t goes only through vertices that are closer to s than t.
- Assertion: "All the problems in graphs with negative weights come from negative weight cycles"
    - in which case distances along the path in negative weight cycles all got to negative infinity
    - so a key idea in finding distances in graphs with negative weights is to identify negative weight cycles

#### Bellman-Ford Algorithm
- algorithm for finding shortest path in graphs where some edges have negative weight (and cannot use Dijkstra Alborithm)
- like Naive Algorithm
    - relax edges where the distance (estimate via upper bound) changes
- Works with negative weights, but _only_ if there are no negative weight cycles in G

    **Bellman–Ford algorithm** for distances from source S in graph G
    BellmanFord(G, S)
    {no negative weight cycles in G}
    for all u∈V:
        dist[u] ← ∞
        prev[u] ← nil
    dist[S] ← 0
    repeat |V|−1 times: # for all u∈V != S
        for all (u,v) ∈ E:
            Relax(u, v )
        while 
            at least one dist changes

- Running time of Bellman-Ford
    - O(|V||E|)

#### Detecting and Dealing with Negative Cycles in a Graph
- Lemma on Existence of Negative Weight Cycles
    - A graph G contains a negative weight cycle if and only if |V|-th (one additional) iteration of BellmanFord(G,S) updates some dist-value.
        - i.e. if "one edge is relaxed"
- Algorithm for finding negative cycles

    **Finding Negative Cycle Algorithm:**
    Run |V| iterations of Bellman–Ford algorithm, 
    save node v relaxed on the last iteration
    v is reachable from a negative cycle
    Start from x ← v (x = Prev(v)), 
        follow the link x ← prev[x] for |V| times # will be definitely on the cycle
        Repeat
            Save y ← x
            go x ← prev[x]
             until x = y again

Week 5 Spanning Trees: Efficient Algorithms
---------------------------

### Minimum Spanning Trees
- idea is from an (Connected, Undirected) Graph to build a network of selected edges that connects all nodes (is "spanning" over all nodes) that is efficient in the sense of having "minimum" total edge weight.
- so it is minimum and spanning
- **Minimum spanning tree (MST)**
    - Input: A connected, undirected graph G = (V,E) with positive edge
weights.
    - Output: A subset of edges E' ⊆ E of **minimum total weight** such that
the graph G(V, E') is **connected**.
    - Remark: The set E' always forms a **tree**.
- Properties of a tree
    - tree = undirected, connected, acyclic graph
    - tree with N vertices has N-1 edges
    - any connected, acyclic graph with N-1 edges is a tree
    - An undirected graph is a tree iff there is a unique path between any pair of its vertices 
        - (i.e. it is connected (at least one path) and acyclic (at most one path)).
- Note: A Min Spanning Graph must be acyclic (since for any cycle could remove one edge and still be connected with lower total weight), so it must be a tree

### Greedy Algorithms for Minimum Spanning Tree Problem
- two efficient Greedy Alg for MST problem: Kruskal and Prim
    - **Kruskal**: repeatedly add the next edge that is the lightest edge from G that doesn’t produce a cycle
        - (edges need not be to an connected node yet)
    - **Prim**: (select a root vertex;) repeatedly select/add the next edge that is the lightest edge that attaches a new vertex to the current tree
- Both Kruskal and Prim Greedy Algs for MST produce optimal soln

#### Cut property (efficiency of adding shortest edge between two partial portions of a MST)
- For a graph V(V,E), 
    - Let 
        - X ⊆ E be a part of a MST of G(V,E),
        - S ⊆ V be (a partition) such that no edge of X crosses between S and V − S, and 
        - e ∈ E be a lightest edge across this partition.
    - Then X + {e} is a part of some MST.

#### Kruskal's Algorithm for MST

    - **Kruskal’s Algorithm**
        - Algorithm: repeatedly add to X the next lightest edge e that doesn’t produce a cycle
        - At any point of time, the set X is a forest, that is, a collection of trees
        - The next edge e connects two different trees—say, T1 and T2
        - The edge e is the lightest between T1 and V − T1, hence adding e is safe

- Kruskal Algorithm for MST Implementation Details
    - use disjoint sets data structure
    - initially, each vertex lies in a separate set
    - each set is the set of vertices of a connected component
    - when checking about adding a new (minimum weight) edge, to check whether the current edge {u,v} produces a cycle, we check whether u and v belong to the same set

- Kruskal Algorithm for MST - Pseudocode

    **Kruskal(G)**

        for all u∈V:
            MakeSet(v ) # separate set for each node
        X ← empty set # no edges in the forest of partial MSTs
        sort the edges E by weight
        for all {u,v} ∈ E in non-decreasing weight order:
            if Find(u) <> Find(v): # confirm that they are not in the same connected component (same set)
                add edge {u,v} to X
                Union(u, v) # indicate that all nodes of u and v are in same CC
        return X

- Kruskal Running Time
    - Running Time
    - Sorting edges:
        - O(|E|log|E|) = O(|E|log|V|2) = O(2|E|log|V|) = O(|E|log|V|) 
    - Processing edges:
        - 2|E | · T (Find) + |V | · T (Union) = O((|E|+|V|)log|V|) = O(|E|log|V|)
    - Total running time: 
        - O(|E|log|V|)
        - less if edges come pre-sorted, then O(|E|log*|V|

#### Prim’s Algorithm
X is always a subtree, grows by one edge at each iteration
we add a lightest edge between a vertex of the tree and a vertex not in the tree
very similar to Dijkstra’s algorithm

- Prim's Algoritym Pseudocode

    **Prim’s Algorithm**
    - Prim(G )
    - for all u∈V:
        cost[u] ← ∞, 
        parent[u] ← nil
    pick any initial vertex u0
    cost[u0] ← 0
    PrioQ ← MakeQueue(V ) {priority is cost} 
    while PrioQ is not empty:
        v ← ExtractMin(PrioQ) for all {v,z} ∈ E:
        if z ∈ PrioQ and cost[z] > w(v, z): 
            cost[z] ← w(v, z), parent[z] ← v 
            ChangePriority(PrioQ , z , cost [z ])

- Prim's Algorithm for MST - Running Time
- Running Time
    - the running time is
    - |V|·T(ExtractMin)+|E|·T(ChangePriority) 
    - for array-based implementation, the running time is 
        - O(|V |2)
    - for binary heap-based implementation, the running time is
        -O((|V| + |E|) log|V|) = O(|E| log|V|)


#### Summary
- Kruskal MST Greedy Alg: 
    - repeatedly add the next lightest edge if this doesn’t produce a cycle; 
    - use disjoint sets to check whether the current edge joins two vertices from different components
- Prim's MST Greedy Alg: 
    - repeatedly attach a new vertex to the current tree by a lightest edge; 
    - use priority queue to quickly find the next lightest edge


Data Structure and Algorithm Visualizations
--------------------------------------------
- [Data Structure and Algorithm Visualizations, Galles, U. San Fran.](https://www.cs.usfca.edu/~galles/visualization/Algorithms.html)

Week 6 Advanced Shortest Paths
--------------------------------

### Bidirectional Search
- Dijkstra goes in “circles”
- Bidirectional search idea can reduce the search space
- Roughly 2x speedup for road networks
   - roughly sqrt(N)
-  Meet-in-the-middle —
- 1000 times faster for social networks

### Bidirectional Dijkstra algorithm
- Review Dijkstra
    - Steps
        - Initialize dist[s] to 0, all other distances to ∞
        - ExtractMin — choose unprocessed u with the smallest dist[u]
        - Process u — Relax the edges outgoing from u
        - Repeat until t is processed
    - ExtractMin step assures 
        - that estimated distance to node extracted at that point is indeed actual min
        - that nodes are generally explored in order of their distance from the source S
- Notion of Reverse Graph
- For graph G, reversed graph G^R has:
    - same vertices V, and 
    - reversed edges E^R: (u,v) in E iff (v, u) in E^R
- Bidirectional Dijkstra
    - Given G, build G^R
    - Start Dijsktra (forward) from s in G (forward)
    - Start Dijkstra (forward) from t in G^R (but effectively backward for original problem)
    - Alternate between Dijkstra steps in G and G^R [is strict alternation important?]
    - Continue until find vertex v which is processed in both G and G^R
    - Compute the shortest path between s and t (but note this need not be the path through v)

- Computing Distance Lemma
    - Let dist[u] be the distance _estimate_ in the forward Dijkstra from s in G and distR[u] — the same in the backward Dijkstra from t in GR. [algorithm is such that all nodes have a distance estimate, but if it has not been processed]
    - After some node v is processed both in G and GR, some shortest path from s to t passes through some node u which is processed either in G, in GR, or both, and d(s,t) = dist[u] + distR[u].

- Bidirectional Dijkstra Pseudocode (long)

    - **BidirectionalDijkstra(G , s , t )**
GR ← ReverseGraph(G)
Fill dist,distR with +∞ for each node dist[s] ← 0, distR [t] ← 0
Fill prev,prevR with None for each node proc ← empty, procR ← empty
do:
    v ← ExtractMin(dist)
    Process(v , G , dist, prev, proc)
    if v in procR:
        return ShortestPath(s, dist, prev, proc, t, . . . ) 
        vR ← ExtractMin(distR)
    repeat symmetrically for vR as for v
while True

**Relax(u, v , dist, prev)**
if dist[v]>dist[u]+w(u,v):
    dist[v] ← dist[u] + w(u, v)
    prev[v] ← u

**Process(u, G , dist, prev, proc)**
for (u,v)∈E(G): 
    Relax(u, v , dist, prev)
proc.Append(u)

**ShortestPath(s, dist, prev, proc, t, distR , prevR , procR )**
distance ← +∞, ubest ← None 
for u in proc+procR:
    if dist[u]+distR[u] < distance:
        ubest ← u
        distance ← dist[u] + distR [u] 
path ← empty
last ← ubest 
while last ̸= s:
    path.Append(last)
    last ← prev[last] 
path ← Reverse(path) 
last ← ubest
while last ̸= t:
    last ← prevR[last]
    path.Append(last) 
return (distance,path)


- Conclusion on Bidirectional Dijkstra
    - Worst-case running time of Bidirectional Dijkstra is the same as for Dijkstra
    - Speedup in practice depends on the graph
    - Memory consumption is 2x to store G and GR
    - You’ll see the speedup on social network graph in the Programming Assignment


- Visualizations
  - [Path Finding Algorithms (A*, Dijkstra, Bi-Directional BFS)](https://www.youtube.com/watch?v=DINCL5cd_w0)

### Advanced Shortest Paths: A-star Algorithm (A*)
- for more, see [A* search algorithm](https://en.wikipedia.org/wiki/A*_search_algorithm)

#### Seeking a Directed Search
- utilizing bounds on path lengths
- Potential Function

- A* Algorithm is Dijkstra with Potentials
- Potential are lower bounds on distance, hence can guide search _direction_ ("Directed Seargh"

- Dijkstra wi th Potentials
Take some potential function 𝜋 Launch Dijkstra algorithm with edge
weights l𝜋
The resulting shortest path is also a
shortest path initially Does any 𝜋 fit us?
For any edge (u,v), the new length l𝜋(u,v) must be non-negative — such 𝜋 is called feasible

- Intuition
𝜋(v) is an estimation of d(v,t) — “how far is it from here to t?”
If we have such estimation, we can often avoid going wrong direction — directed search
Typically 𝜋(v) is a lower bound on d(v,t) I.e.,onarealmapapathfromv tot cannot be shorter than the straight line segment from v to t

- A* ≡ Dijkstra 
    - On each step, pick the vertex v minimizing dist[v] − 𝜋(s) + 𝜋(v)
    - 𝜋(s) is the same for all v, so v minimizes dist[v ] + 𝜋(v ) — the most promising vertex
    - 𝜋(v) is an estimate of d(v, t)
    - Pick the vertex v with the minimum current estimate of d(s,v)+ d(v,t) 
    - Thus the search is "directed"
- Idea is that while this alg could still explore whole graph, it is more likely to reach the target node t earlier, by exploring in more promising directions.

- Subtract the potential (lower bound fn) from the XXX

### Bidirectional A* Algorithm for MSP
- need 2 potential
- need to make them consistent in some way, i.e. to imply the same edge weights (lengths) in both directions (so that bidirectional seach will work)

### Lower Bounds (for Potential Functions)
- feasibility $l_\pi(x,y) \ge 0$
- l_/pi(P) = 
    - l_𝜋(u,v) = l(u, v) − 𝜋(u) + 𝜋(v)
    - 0 ≤ l𝜋(P) = l(P) − 𝜋(v) + 𝜋(t) ≤ l(P) − 𝜋(v) 
        - ⇒ 𝜋(v) ≤ l(P) = d(v,t)
- Euclidean Distance Potential
    - feasible and \pi(t) = 0
- Landmark Distance Potential
    - Lemma
        - Fix some vertex A ∈ V , we will call it a landmark. 
        - Then the "Landmark" potential 𝜋(v) = d(A, t) − d(A, v) 
            - is feasible, and 𝜋(t) = 0.
    - Landmark selection?
        - On borders of map
    - Landmarks require preprocessing: need distance from landmark(s) to all other vertices
        - worthwhile if expect many queries on distances between other nodes
    
- A* with Landmarks

- Conclusion
    - Directed search can scan fewer vertices
    - A* is a directed search algorithm based on Dijkstra and potential functions
    - A* can also be bidirectional 
    - Euclidean distance is a potential for aplane (road networks)
    - Landmarks can be used for good potential function, but we need preprocessing to use them