# community-detection
Louvain's community detection algorithm for big graphs

## The algorithm

This is a slight modification of Louvain's algorithm based on the [Fast unfolding of communities in large networks](https://arxiv.org/abs/0803.0476) paper.
For optimizitaion, a metric Q is used.
```
Q = modularity(C) + regularization(C), where C = # of communities
```
The algorithm has the following structure:
```
> Generate `REPS` partitions, for each:
>   while there is an improvement of modularity:
>     Phase1: optimize modularity greedily (with a probability of 1/3 that only c changes are made).
>     Phase2: reconstruct the graph as described in the paper.
>     Save the value of Q and the partition.
>   Pick the partition for the best Q.
> Pick the partition that gives the best Q.
```
The idea behind the use of `MAX_STEPS` for modularity optimization, is that combining  too many communities
is penalized by regularization. This is a particular case when clear communities are presented in the graph.
In order to prevent us from "overoptimizing" modularity (at the expense of regularization) and "overcombining" communities,
less steps are taken for the phase 1.

When testing and tuning `REPS` and `MAX_STEPS`, a value of Q = 2.722822 was achieved.

### What is done and what could have been done
- [x] **Speed**
      
      Since all nodes in the graph are in the range [0, V), where V = # vertices, arrays/vectors
      are used for indexing. This gives a significant sppedup when finding neighbouring communities.
      The overall time commplexity is:
            O(M^2 * V)
            V = # vertices.
            M has an upper bound of (log V) and depends whether only `MAX_STEPS` where made.
      The overall space complexity is:
            O(V + E) (storing graph)
            E = # edges.

- [x] **Randomization and maximizing**
      
      In order to obtain better results, the vertices selection is randomized. This gives a
      significant increasу in Q (around 0.2-0.3 for some graphs). Also, the algorithm is executed
      a fixed number of times. Therefore, we are able to pick the best result from the randomly
      formed partitions.

- [ ] **Style**
  
      Classes can be factored out into separate files for readability. Data structure information
      printing can be allowed in debug mode only and more comments can be added, describing the
      functions/algorithms.

- [ ] **Variation**

      After experimenting with different graphs it is clear that there is no algorithm that
      would suit all cases. Therefore, a better approach may be:
      
          1. Detect graph “type” - whether it’s sparse, dense, with clear communities, etc.
          This may be done via adjacency matrix construction and analyzing weight density.
          
          2. For each type choose the best suitable metric/algorithm. An alternative to
          greedy Girvan–Newman or Louvain would be linear algebra or statistical inference
          based algorithms.
      
      Also, speaking of the implemented algorithms, instead of pure randomization it may be
      better to randomize only a certain amount of time, and other proportion cluster nodes
      greedily.
      
- [ ] **Start point and heuristics**

      The algorithm starts with every node being in a single community. Whereas this is a
      reasonable starting point, it is clear that some inspection with heuristics in the
      beginning may allow to identify and cluster together potential communities/“centroid”
      vertices prior to the main execution.

      An example of such a heuristic could be density analysis for a chosen N connected vertices.
      Another heuristic that can be used is that if we have clear communities, the edges between
      them are rare. Therefore we DO NOT want to cluster nodes together if their communities have
      only a small number of edges connecting each other.

## Execution
Input is the file of edges:
```
src dst
```
Output is a partition produced by the algorithm (saved to "partion.graph") and the value of Q printed to stdout.
For example:
```
g++ -std=c++11 -O3 -Wall main.cpp -o main
main tests/1.graph
1.23532
```

