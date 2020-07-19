Questions_on_Graph_Algorithms_to_Discuss.md
=================================================

Week 3, DFS, Node distances and Bipartite Graph
-----------------------------------------------
- if use (ordered) list to contain node attributes (visited, distance, color),
    - do we need to be careful about node numberings that are not numbered sequentially from 1?
- Does Source node matter for BFS approach to alternate coloring?
    - it does for numbering, in terms of numbering level
    - but maybe not for relative distance, i.e. if number is even or odd?
- HOw does bipartite test deal with non-connected graphs?
- Graph plotting/visualization
    - blip with note labelling in matplotlib
- What about shortest path tree algorith: how do we use that?
    - see "10_shortest_paths_in_graphs_1_bfs.pdf" slide 196
    - the "prev" concept is needed for shortest (weighted) path algs

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

    - for measured distance "dist[v]" from source S to v:

        **Relax((u,v) ∈ E)**
        if dist[v] > dist[u] + w(u,v): 
            dist[v] ← dist[u] + w(u, v)
            prev[v] ← u   # record closest know predecessor node as `prev`
    
    - Intuition: 
        - edge (u,v) is "relaxed" in the context of distances from source node S
        - edge (u,v) is relaxed in the sense that we are looking for a shortcut to v (from S) via u
            - "relax edge (u,v)" means "attempt to shorten v's distance to s using u"
            - in OR, relaxation means relaxing a constraint
            - Naive and in Dijkstra relaxation means relaxation of the constraint x_uv = 0?
                - removing (u,v) from set of edges excluded from the shortest path??

        **NaiveDist(G , S)**
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
    - Best to explore (relax edges from) that set of known-distance nodes, 
        - since any node w outside R with new minimum estimate dist[w] will also be its shortest path (correct distance d(s,w))
    - On each iteration we 
        - take a vertex outside of R with the minimal dist-value, 
        - add it to R, and 
        - relax all its outgoing edges [if digraph, o.w. all its edges?].
    - Start with including only source vertex S in R (since dist[S] = 0)
- Dijkstra Algorithm Pseudocode

        **Dijkstra(G, S)**
        for all u∈V:
            dist[u] ← ∞, prev[u] ← nil
        dist[S] ← 0
        H ← MakeQueue(V ) {dist-values as keys} # this is the Unknown region, not(R)
        while H is not empty:
            u ← ExtractMin(H)
            # Lemma: When a node u is selected via ExtractMin, dist[u] = d(S,u).
            for all (u,v)∈ E: # relax _outgoing_ edges from u
                if dist[v] > dist[u]+w(u,v):
                    dist[v] ← dist[u] + w(u, v)
                    prev[v] ← u
                    ChangePriority(H , v , dist [v])

#### Topics for Discussion 20200717
- Relation of LAP to bipartite graph
- relation of MIP to graph shortest path problems
- Performance of MIP LAP solution (see R code results)
- Running time of various implementations for Dijkstra steps?
    - priority queue?

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