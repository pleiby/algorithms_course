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
