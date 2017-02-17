// Copyright 2017 <Biagio Festa>
#include <vector>
#include <set>
#include "TypeUtility.hpp"
#include "ProblemDatas.hpp"

bool Solution::compute_feasibility(
    const for_ch::ProblemDatas& problem) const noexcept {
    // TODO(biagio): incomplete function
  const unsigned num_nodes = problem.m_numNodes;
  const unsigned num_nodes_minusone = num_nodes - 1;

  std::set<VertexIndex> inner_vertices;
  std::vector<std::set<VertexIndex>> nexts(num_nodes);
  std::vector<std::set<VertexIndex>> prevs(num_nodes);
  std::vector<RealNumber> walked_distance(num_nodes, 0);
  bool has_outgoing_edge;
  EdgeIndex index_e;
  VertexIndex target;
  RealNumber delta_d;

  // For each node except school
  for (VertexIndex n = 1; n < num_nodes; ++n) {
    // Compute starting edge index
    index_e = num_nodes_minusone * (n - 1);

    has_outgoing_edge = false;

    // For each edge of that node
    for (unsigned e = 0; e < num_nodes_minusone; ++e, ++index_e) {
      assert(index_e < m_active_edges.size());
      if (m_active_edges[index_e] == true) {
        if (has_outgoing_edge == true) {
          return false;
        }
        has_outgoing_edge = true;

        // Since this is the first outgoing edge the walked distance is 0
        assert(walked_distance[n] == 0);

        // Get the target of that edge
        target = problem.m_mapEdge_index2link.at(index_e).second;

        // Insert target as inner
        inner_vertices.insert(target);

        // Compute the delta distance
        delta_d = problem.m_distance_matrix[n][target] +
            walked_distance[target];

        // Set the upper side of the graph
        nexts[n].insert(target);
        walked_distance[n] += delta_d;
        AddSet(nexts[target], &nexts[n]);
        for (const VertexIndex upper : prevs[n]) {
          AddSet(nexts[n], &nexts[upper]);
          walked_distance[upper] += delta_d;
        }

        // Set the own side of the graph
        prevs[target].insert(n);
        AddSet(prevs[n], &prevs[target]);
        for (const VertexIndex downer : nexts[target]) {
          AddSet(prevs[target], &prevs[downer]);
        }
      }  // if edge is active
    }  // for each edge of that node
  }  // for each node except school

  return true;
}
