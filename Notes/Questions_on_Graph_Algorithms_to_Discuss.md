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
