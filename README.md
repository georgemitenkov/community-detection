# community-detection
Louvain's community detection algorithm for big graphs

## The algorithm

This is an implementation of Louvain's algorithm based on the [Fast unfolding of communities in large networks](https://arxiv.org/abs/0803.0476) paper.
For optimizitaion, a metric Q is used.
```
Q = modularity(C) + regularization(C), where C = # of communities
```
When testing, a value of Q = 2.637893 was achieved.

### What is done and what could have been done
- [x] **Speed**
      
      Since all nodes in the graph are in the range [0, V), where V = # vertices, arrays/vectors are used for indexing.
      This gives a significant sppedup when finding neighbouring communities. The overall commplexity is ...

- [x] **Randomization and maximizing**
      
      In order to obtain better results, the vertices selection is randomized. This gives a significant increase in Q
      (around 0.2-0.3 for some graphs). Also, the algorithm was executed a fixed number of times. Therefore, we are
      able to pick the best result from the randomly formed partitions.

- [ ] **Style**
  
      Classes can be factored into separate files for readability. Data structure information printing can be allowed
      in debug mode.

- [ ] **Variation**

      After experimenting with different graphs it is clear that there is no algorithm that would suit all cases.
      Therefore, a better approach may be:
      
          1. Detect graph “type” - whether it’s sparse, dense, with clear communities, etc. This may be done via
          adjacency matrix construction and analyzing weight density.
          
          2. For each type choose the best suitable metric/algorithm. An alternative to greedy Girvan–Newman or
          Louvain would be linear algebra or statistical inference based algorithms.
      
      Also, speaking of the implemented algorithms, instead of pure randomization it may be better to randomize only
      a certain amount of time, and other proportion cluster nodes greedily.
      
- [ ] **Start point and heuristics**

      The algorithm starts with every node being in a single community. Whereas this is a feasible starting point,
      it is clear that some inspection with heuristics in the beginning may allow to identify and cluster together
      potential communities/“centroid” vertices.

      An example of such a heuristic could be density analysis for a chosen N connected vertices. Another heuristic
      that can be used is that if we have clear communities, the edges between them are rare. Therefore we DO NOT
      want to cluster nodes together if their communities have only a small number of edges connecting each other.

## Execution
Input is the file of edges:
```
src dst
```
Output is a partition produced by the algorithm (saved to "partion.graph").
For example:
```
g++ -std=c++11 -O3 -Wall main.cpp -o main
main tests/karate.graph
...
```

