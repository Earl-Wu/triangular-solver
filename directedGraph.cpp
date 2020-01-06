// This file contains a simple directed graph data structure and related helper
// functions.

// The topological sort algorithm is based on https://www.geeksforgeeks.org/topological-sorting/

#include<iostream>
#include <stdlib.h>
#include "directedGraph.h"
using namespace std;

Graph::Graph() {
  V = 0;
}

void Graph::addEdge(int v, int w)
{
  if (adj.find(w) == adj.end()) {
    vector<int> empty_entry = {};
    adj[w] = empty_entry;
    this->V += 1; // Add w to v’s list.
  }

  if (adj.find(v) != adj.end()) {
    adj[v].push_back(w); // Add w to v’s list.
  }

  else {
    vector<int> new_entry = {w};
    adj[v] = new_entry;
    this->V += 1;
  }
}

// A recursive function used by topologicalSearch
void Graph::topologicalSortUtil(int v, bool visited[], stack<int> &Stack)
{

  // Mark the current node as visited.
  visited[v] = true;

  // Recur for all the vertices adjacent to this vertex
  vector<int>::iterator i;
  for (i = adj[v].begin(); i != adj[v].end(); ++i) {
      if (!visited[*i]) {
          topologicalSortUtil(*i, visited, Stack);
      }
  }

  // Push current vertex to stack which stores result
  Stack.push(v);
}

/* The main topologicalSearch function
 * Adapted from https://www.geeksforgeeks.org/topological-sorting/
 * @parameter: A vector containing the starting points of the search
 * @return: A vector containing the reachability status from each starting point */
vector<stack<int>> Graph::topologicalSearch(vector<int> sources)
{

  // Mark all the vertices as not visited
  bool *visited = new bool[V];
  for (int i = 0; i < V; i++) {
      visited[i] = false;
  }

  // Call the recursive helper function to store topological
  // sort results starting from each initial source
  vector<stack<int>> stacks;
  for (int i = 0; i < (int) sources.size(); i++) {
    if (visited[sources[i]] == false) {
      stack<int> cur_stack;
      topologicalSortUtil(sources[i], visited, cur_stack);
      stacks.push_back(cur_stack);
    }
  }

  return stacks;
}

// A helper function to print vector<int>
void print_vector(vector<int> v) {
  for (int i = 0; i < (int) v.size(); i++) {
    cout << v[i] << " ";
  }
}

// A helper function to print map (assuming the value is of type vector<int>)
template<typename K, typename V>
void print_map(std::unordered_map<K,V> const &m)
{
    for (auto const& pair: m) {
        std::cout << "{" << pair.first << ": ";
        print_vector(pair.second);
        cout << "}\n";
    }
}

/* A function that will construct a dictionary containing dependencies information
 * that's ready for wavefront parallelization
 * @parameter: A vector containing dependencies paths for each initial nonzero values
 * @return: A dictionary containing sublevels such that each entry in the sublevel
 * 1) co-indepedent 2) has all preliminary values needed for its computation ready */
unordered_map<int, vector<int>> wavefront(vector<stack<int>> stacks) {
  unordered_map<int, vector<int>> wavefront_vector;
  int counter;
  unordered_map<int, vector<int>> rec;
  for (int i = 0; i < (int) stacks.size(); i++) {
    counter = 0;
    while (stacks[i].empty() == false) {
      int cur_node = stacks[i].top();
      if (wavefront_vector.find(counter) != wavefront_vector.end()) {
        wavefront_vector[counter].push_back(cur_node);
      }
      else {
        vector<int> temp = {cur_node};
        wavefront_vector[counter] = temp;
      }
      counter++;
      stacks[i].pop();
    }
  }

  // DEBUG:
  // for (int i = 0; i < (int) wavefront_vector.size(); i++) {
  //   vector<int> cur = wavefront_vector[i];
  //   cout << i << ": ";
  //   for (int j = 0; j < (int) cur.size(); j++) {
  //     cout << cur[j] << " ";
  //   }
  //   cout << "\n";
  // }
  // print_map(wavefront_vector);

  return wavefront_vector;
}
