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

Week 4, Dijkstra and Bellman-Ford Algorithm

#### Topics for Discussion 20200717
- Relation of LAP to bipartite graph
- relation of MIP to graph shortest path problems
- Performance of MIP LAP solution (see R code results)
- Running time of various implementations for Dijkstra steps?
    - priority queue?


Wee 4 continued:
- Bellman-Ford proof and intuition
- Julia HPC approaches
- Julia

### Week 6: Advanced Shortest Paths

#### Bidirectional Dijkstra

#### Advanced Shortest Paths: A-star Algorithm (A*)
    - bidirectional search vs. directional search
    - Dijkstra with Potentials
        - Take some potential function 𝜋 over nodes
        - Launch Dijkstra algorithm with edge weights l𝜋
            - l𝜋(u, v) = l(u,v) - 𝜋_u + 𝜋_v for adjacent vertices
            - l𝜋(s, t) = l(s, t) - 𝜋_s + 𝜋_t
    - The resulting shortest path is also a shortest path initially 
    - Does any 𝜋 fit us?
    - For any edge (u,v), the new length l𝜋(u,v) must be non-negative — such 𝜋 is called feasible
    - Intuition
        - 𝜋(v) is an estimation of d(v,t) — “how far is it from here to t?”
        - If we have such estimation, we can often avoid going wrong direction — directed search
        - Typically 𝜋(v) is a lower bound on d(v,t) (remainig distance to t)

#### Bidirectional A*
    - Same as Bidirectional Dijkstra, but with
potentials
    - Needs two potential functions: 
        - 𝜋_f(v) estimates d(v,t), (l.b. est on forward dist to terminal node)
        - 𝜋_r(v) estimates d(s,v) (l.b. est on distance back to source)

#### Lower Bounds

#### Landmarks
- Landmarks
    - Select several landmarks and precompute their distances to all other vertices
- For any landmark A, use triangle inequality to establish distance bounds
    - triangle inequality
        - d(A,v) + d(v,t) ≥ d(A,t), 
        - d(v,t) + d(t,A) ≥ d(v,A) 
    - so
        - d(v,t) ≥ d(A,t) − d(A,v), 
        - d(v,t) ≥ d(v,A) − d(t,A)

#### Conclusion
- Directed search can scan fewer vertices
- A* is a directed search algorithm based on Dijkstra and potential functions
- A* can also be bidirectional Euclidean distance is a potential for a
plane (road networks)
- Landmarks can be used for good potential function, but we need preprocessing to use them

#### Advanced Shortest Path Algorithms
- interesting that Dijstra's Algorithm is essentially guaranteed to search outward circularly (in regions of monotonically increasing distance from S)
    - because ExtractMin step that chooses node to explore next assures that nodes are generally explored in order of their distance from the source S

- Q: if bidirectional search covers only 1/2 the number of nodes as unidirectional Dijkstra, and the 6-degrees o separation hypothesis is validwhy would be confident that we can make any social network connection between any pair of people with onl 3 degrees (friend-of-friends-of-friend) in each direction?
    - slides 19, slide 67, Six andshakes video, minute 3:00
- "Meet-in-the-Middle strategy for any optimization or shortest path
    - solves in sqrt(n) rather than n    
  
- Computing Distance Lemma - 
    - why is distance approximation dist^R[u] both the upper and lower bound of actual distance d(u, t)? I.e., why must the approximate distance for unprocessed nodes be exact if those nodes are processed in the other direction, and on the shortest path?
    - Lecture, 8:15, Finding Shortest path after meertoign in the middle

- Intuition on Landmark selection?
    - On borders of map

Points to Discuss on Advanced Shortest Path Algs
------------------------------------------------
- Large networks with non-sequential node ID numbers and sparse connections (e.g. streets)

study group notes 20200912: objectives for next time
------------------------------------------------
- "capstone" for Graphs
    - implement bidirectional dijkstra (if possible, understand early stopping condition!)
    - (hopefully) A* with Euclidean distance for potential
    - Apply either of these to the SF (500K) or Austin (1M edges) data sets
- Think about next course to start (~after September, 1-2 more graphs meetings), e.g.
https://computationalthinking.mit.edu/Fall20/
