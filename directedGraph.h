#include <stdlib.h>
#include <vector>
#include <stack>
#include <unordered_map>
using namespace std;

typedef unordered_map<int, vector<int>> intmap;

// Class to represent a graph
class Graph {
    // Pointer to an array containing adjacency listsList
  intmap adj;

  // A function used by topologicalSort
  void topologicalSortUtil(int v, bool visited[], stack<int> &Stack);

public:
  Graph();   // Constructor

  int V;    // No. of vertices'
   // function to add an edge to graph
  void addEdge(int v, int w);

  // prints a Topological Sort of the complete graph
  stack<int> topologicalSort();

	vector<stack<int>> topologicalSearch(vector<int> sources);

};

unordered_map<int, vector<int>> wavefront(vector<stack<int>> stacks);
