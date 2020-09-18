#include <algorithm>
#include <iostream>
#include <iterator>
#include <ctime>
#include <map>
#include <vector>

using namespace std;

static const int REPS = 15;
static const int MAX_STEPS = 10000;

//============================================================================//
// Initial graph G(V, E)
//============================================================================//

class InitGraph {

public:
  int V;
  int E;
  vector<vector<int>> adj;

  InitGraph(int V, int E, vector<pair<int, int>> &edges) {
    this->V = V;
    this->E = 2 * E;
    adj.resize(V);
    for (auto e : edges) {
      adj[e.first].push_back(e.second);
      adj[e.second].push_back(e.first);
    }

    for (int i = 0; i < V; i++) {
      vector<int> copy = adj[i];
      sort(copy.begin(), copy.end());
      adj[i] = copy;
    }
  }

  // DEBUG
  void print() {
    cout << "Initial graph G(V: " << V << ", E: "
         << E << ") with adjacency list:\n";
    for (int i = 0; i < V; i++) {
      cout << "("<< i << ") ";
      for (auto v : adj[i]) {
        cout << " -> " << v;
      }
      cout << "\n";
    }
    cout << "\n";
  }

};

//============================================================================//
// Graph network implementation G(V, E, W).
//============================================================================//

class Graph {
public:
  int V;                    // Number of vertices.
  int E;                    // Number of edges.
  int W;                    // Total weight.

  vector<int> degrees;      // A vector of cumulitive degrees.
  vector<int> edges;        // A vector of all `dst` from adjacency list.
  vector<int> weights;      // A vector of weights.

  Graph(){
    V = 0;
    E = 0;
    W = 0;
  };
  Graph(InitGraph &g) {
    V = g.V;
    E = g.E;
    W = E;

    degrees.resize(V);
    int sum = 0;
    for (unsigned i = 0; i < V; i++) {
      sum += g.adj[i].size();
      degrees[i] = sum;
    }

    int edgeCounter = 0;
    for (unsigned int i = 0; i < V; i++) {
      for (int j = 0; j < g.adj[i].size(); j++) {
        edges.push_back(g.adj[i][j]);
        edgeCounter++;
      }
    }

    weights.resize(edgeCounter);
    fill(weights.begin(), weights.end(), 1);
  }

  int degree(int v) {
    if (v == 0)
      return degrees[0];
    return degrees[v] - degrees[v - 1];
  }

  pair<vector<int>::iterator, vector<int>::iterator>
  adjacent(int v) {
  if (v == 0)
    return make_pair(edges.begin(), weights.begin());
  else if (weights.size() != 0)
    return make_pair(edges.begin() + degrees[v - 1],
                     weights.begin() + degrees[v - 1]);
  else
    return make_pair(edges.begin() + degrees[v - 1], weights.begin());
  }

  int selfLoops(int v) {
    auto pointer = adjacent(v);
    for (int i = 0; i < degree(v); i++) {
      if (*(pointer.first + i) == v) {
        if (weights.size() != 0)
	        return *(pointer.second + i);
        else 
	        return 1;
      }
    }
    return 0;
  }

  int weightedDegree(int v) {
    if (weights.size() == 0) {
      return degree(v);
    } else {
      auto pointer = adjacent(v);
      double sum = 0;
      for (int i = 0; i < degree(v); i++) {
        sum += *(pointer.second + i);
      }
      return sum;
    }
    
  }

  // DEBUG
  void print() {
    cout << "Graph G(V: " << V << ", E: " << E << ", W: " << W << ")\ndegrees: ";
    for (int i = 0; i < V; i++) {
      cout << degrees[i] << " ";
    }
    cout << "\nedges: ";
    for (int i = 0; i < edges.size(); i++){
      cout << edges[i] << " ";
    }
    cout << "\n";
    cout << "weights: ";
    for (int i = 0; i < weights.size(); i++) {
      cout << weights[i] << " ";
    }
    cout << "\n\n";
  }
};

//============================================================================//
// Community
//============================================================================//

class Community {
public:

  int originalV; // Original number of vertices in the graph.

  vector<int> adjWeights;
  vector<int> adjPositions;
  int adjLast;

  Graph graph;
  int networkSize;
  vector<int> vertexToCommunityMap;
  vector<int> ein, eout, c;

  Community(Graph g, int v, vector<int> &communitySizes) {
    originalV = v;
    graph = g;
    networkSize = g.V;

    adjWeights.resize(networkSize, -1);
    adjPositions.resize(networkSize);
    adjLast = 0;

    vertexToCommunityMap.resize(networkSize);
    ein.resize(networkSize);
    eout.resize(networkSize);
    c.resize(networkSize);

    for (int i = 0; i < networkSize; i++) {
      vertexToCommunityMap[i] = i;
      ein[i]  = g.selfLoops(i);
      eout[i] = g.weightedDegree(i);
      c[i] = communitySizes[i];
    }
  }

  void remove(int vertex, int community, double edgesToCommunity) {
    eout[community] -= graph.weightedDegree(vertex);
    ein[community] -= 2 * edgesToCommunity + graph.selfLoops(vertex);
    vertexToCommunityMap[vertex] = -1;
  }

  void insert(int vertex, int community, int edgesToCommunity) {
    eout[community] += graph.weightedDegree(vertex);
    ein[community] += 2 * edgesToCommunity + graph.selfLoops(vertex);
    vertexToCommunityMap[vertex] = community;
  }

  double modularity() {
    double mod = 0.;
    for (int i = 0; i < networkSize; i++) {
      if (eout[i] > 0) {
        mod += (double)ein[i] / (double)graph.W
            - ((double)eout[i] / (double)graph.W)
            * ((double)eout[i] / (double)graph.W);
      }
    }
    return mod;
  }

  double regularization() {
    int densitySum = 0;
    for (int i = 0; i < networkSize; i++) {
      if (c[i] == 1) {
        densitySum++;
      } else {
        double den = 0.5 * (double)c[i] * (c[i] - 1);
        densitySum += (double)ein[i] / den;
      }
    }
    double n = static_cast<double>(networkSize);
    double V = static_cast<double>(originalV);
    return 0.5 * (static_cast<double>(densitySum) / n - (n / V));
  }

  // numEdges from v to community.
  double modularityGain(int v, int newCommunity, int numEdges, int degree) {
    double eouts = eout[newCommunity];
    double deg = degree;
    double m = graph.W;
    
    return ((double)numEdges - eouts * deg / m);
  }

  void adjCommunities(int vertex) {
    for (int i = 0; i < adjLast ; i++)
      adjWeights[adjPositions[i]] = -1;
    adjLast = 0;

    auto p = graph.adjacent(vertex);

    int deg = graph.degree(vertex);

    adjPositions[0] = vertexToCommunityMap[vertex];
    adjWeights[adjPositions[0]] = 0;
    adjLast = 1;

    for (int i = 0; i < deg; i++) {
      int neighbour = *(p.first + i);
      int neighbourCommunity = vertexToCommunityMap[neighbour];
      double neighbourWeight = (graph.weights.size() == 0) ? 1 : *(p.second + i);
      
      if (neighbour != vertex) {
        if (adjWeights[neighbourCommunity] == -1) {
          adjWeights[neighbourCommunity] = 0;
          adjPositions[adjLast] = neighbourCommunity;
          adjLast++;
        }
        adjWeights[neighbourCommunity] += neighbourWeight;
      }
    }
  }

  bool phase1(bool dontLimit) {
    bool hasChanged = false;
    int moves;
    double newMod = modularity();
    double currentMod = newMod;

    vector<int> randomOrder(networkSize);
    for (int i = 0; i < networkSize; i++)
      randomOrder[i] = i;
    for (int i = 0; i < networkSize; i++) {
      int randPosition = rand() % (networkSize - i) + i;
      int tmp = randomOrder[i];
      randomOrder[i] = randomOrder[randPosition];
      randomOrder[randPosition] = tmp;
    }

    int counter = 0;

    do {
      currentMod = newMod;
      moves = 0;

      for (int i = 0; i < networkSize; i++) {
        int vertex = randomOrder[i];
        int vertexCommunity = vertexToCommunityMap[vertex];
        double vertexDegree = graph.weightedDegree(vertex);

        adjCommunities(vertex);
        remove(vertex, vertexCommunity, adjWeights[vertexCommunity]);

        int bestCommunity = vertexCommunity;
        int bestNumEdges = 0;
        double bestIncrease = 0.0;
        for (int j = 0 ; j < adjLast; j++) {
          double increase = modularityGain(vertex, adjPositions[j], adjWeights[adjPositions[j]], vertexDegree);
          if (increase > bestIncrease) {
            bestCommunity = adjPositions[j];
            bestNumEdges = adjWeights[adjPositions[j]];
            bestIncrease = increase;
          }
        }

        insert(vertex, bestCommunity, bestNumEdges);
      
        if (bestCommunity != vertexCommunity)
          moves++;

        counter++;

      }

      newMod = modularity();
      if (moves > 0)
        hasChanged = true;

    } while ((dontLimit || counter < MAX_STEPS) && moves > 0 && newMod > currentMod);

    return hasChanged;
  }

  Graph phase2(vector<int> &communitySizes) {
    vector<int> renumber(networkSize, -1);
    for (int v = 0; v < networkSize; v++) {
      renumber[vertexToCommunityMap[v]]++;
    }

    int final = 0;
    for (int i = 0; i < networkSize; i++) {
      if (renumber[i] != -1) {
        renumber[i] = final;
        final++;
      }
    }

    vector<vector<int>> communityVertices(final);
    for (int v = 0; v < networkSize; v++) {
      communityVertices[renumber[vertexToCommunityMap[v]]].push_back(v);
    }

    Graph g;
    int numCommunities = communityVertices.size();
    communitySizes.resize(numCommunities);
    g.V = numCommunities;
    g.degrees.resize(numCommunities);

    for (int i = 0; i < numCommunities; i++) {
      int sum = 0;
      for (int j = 0; j < communityVertices[i].size(); j++) {
        sum += this->c[communityVertices[i][j]];
      }
      communitySizes[i] = sum;
    }


    for (int c = 0; c < numCommunities; c++) {
      map<int, int> m;
      map<int, int>::iterator it;

      int communitySize = communityVertices[c].size();
      for (int v = 0; v < communitySize; v++) {
        auto p = graph.adjacent(communityVertices[c][v]);
        int deg = graph.degree(communityVertices[c][v]);
        for (int i = 0; i < deg; i++) {
          int neighbour = *(p.first + i);
          int neighbourCommunity = renumber[vertexToCommunityMap[neighbour]];
          double neighbourWeight = (graph.weights.size() == 0) ? 1 : *(p.second + i);

          it = m.find(neighbourCommunity);
          if (it == m.end())
            m.insert(make_pair(neighbourCommunity, neighbourWeight));
          else
            it->second += neighbourWeight;
        }
      }

      g.degrees[c] = (c == 0) ? m.size() : g.degrees[c - 1] + m.size();
      g.E += m.size();

      for (it = m.begin(); it != m.end(); it++) {
        g.W += it->second;
        g.edges.push_back(it->first);
        g.weights.push_back(it->second);
      }
    }
    return g;
  }

  void partition(vector<int> &result) {
    vector<int> renumber(networkSize, -1);
    for (int v = 0; v < networkSize; v++) {
      renumber[vertexToCommunityMap[v]]++;
    }

    int final = 0;
    for (int i = 0; i < networkSize; i++) {
      if (renumber[i] != -1) {
        renumber[i] = final;
        final++;
      }
    }

    for (int i = 0; i < networkSize; i++)
      result[i] = renumber[vertexToCommunityMap[i]];
  }

  // DEBUG
  void print() {
    cout << "Community:\n";
    cout << "Network size = " << networkSize << "\n";
    cout << "vertexToCommunityMap: ";
    for (int i = 0; i < networkSize; i++) {
      cout << vertexToCommunityMap[i] << " ";
    }
    cout << "\nein: ";
    for (int i = 0; i < networkSize; i++) {
      cout << ein[i] << " ";
    }
    cout << "\neout: ";
    for (int i = 0; i < networkSize; i++) {
      cout << eout[i] << " ";
    }
    cout << "\nc: ";
    for (int i = 0; i < networkSize; i++) {
      cout << c[i] << " ";
    }
    cout << "\n";
    cout << "modularity: " << modularity() << "\n";
    cout << "reg: " << regularization() << "\n";
  }

};

//============================================================================//
// I/O.
//============================================================================//

pair<int, int> readGraphData(char *filename, vector<pair<int, int>> &edges) {
  freopen(filename, "r", stdin);
  int E = 0, V = 0;

  int i, j;
  while (cin >> i >> j) {
    E++;
    V = max(V, max(i, j));
    edges.push_back({i, j});
  }
  V++;
  fclose(stdin);
  return make_pair(V, E);
}

//============================================================================//
// Main.
//============================================================================//

int main(int argc, char *argv[]) {
  vector<pair<int, int>> edges;
  auto dim = readGraphData(argv[1], edges);
  const int V = dim.first;
  int E = dim.second;
  InitGraph g(V, E, edges);

  double predQ = 0;
  map<int, vector<int>> best;
  for (int rep = 0; rep < REPS; rep++) {
    Graph graph(g);
    vector<int> cs(V, 1);
    Community c(graph, V, cs);

    bool hasChanged = true;
    double mod = c.modularity(), newMod;
    double reg = 0.0, newReg;
    vector<double> Qs;
    vector<vector<int>> levels;
    srand((unsigned) time(0));

    bool dontLimit = rand() % 3;
    do {
      
      // Phase 1: Modularity optimization.
      hasChanged = c.phase1(dontLimit);
      newMod = c.modularity();

      // Phase 2: Community aggregation.
      vector<int> result(c.networkSize);
      c.partition(result);
      levels.push_back(result);
      vector<int> communitySizes;
      graph = c.phase2(communitySizes);

      c = Community(graph, V, communitySizes);
      newReg = c.regularization();
      Qs.push_back((newMod + newReg));

      cout << "m: " << newMod << ", r = " << newReg << "\n";
    
      mod = newMod;
      reg = newReg;
    } while (hasChanged);

    pair<int,double> maxQ = {0, 0.0};
    for (int q = 0; q < Qs.size(); q++) {
      if (Qs[q] > maxQ.second)
        maxQ = {q, Qs[q]};
    }

    cout << "Q max is: " << maxQ.second << "\n";

    if (maxQ.second > predQ) {
      vector<int> vertexToCommunity(V);
      map<int, vector<int> > communities;

      for (int i = 0 ; i < V; i++)
        vertexToCommunity[i]=i;
      
      for (int l = 0 ; l < maxQ.first; l++) {
        for (int v = 0 ; v < V; v++)
          vertexToCommunity[v] = levels[l][vertexToCommunity[v]];
      }
      
      for (int v = 0; v < V; v++) {
        communities[vertexToCommunity[v]].push_back(v);
      }

      predQ = maxQ.second;
      best = communities;
    }
  }

  cout << predQ << "\n";

  freopen("partition.graph", "w", stdout);
  for (auto c : best) {
    for (auto it = c.second.begin(); it != c.second.end(); it++) {
      if (it != c.second.begin())
        cout << " ";
      cout << *it;
    }
    cout << "\n";
  }
  fclose(stdout);
}