// Copyright 2017 <Biagio Festa>
#include <memory>
#include <utility>
#include <vector>
#include <algorithm>
#include "HESolver.hpp"

namespace for_ch {

HESolver::HESolver(const std::shared_ptr<ProblemDatas>& problem) noexcept :
    mp_problem(problem) {
}

bool HESolver::run(std::vector<bool>* active_edges) {
  assert(active_edges != nullptr);
  active_edges->resize(mp_problem->m_numEdges, false);

  bool found_solution = launchProblemSolver<HESolver::Heuristic>();

  if (found_solution) {
    for (const auto& path : m_foundPaths) {
      const unsigned len_path_minusone = path.m_pathVertices.size() - 1;
      for (unsigned i = 0; i < len_path_minusone; ++i) {
        const VertexIndex source = path.m_pathVertices[i];
        const VertexIndex target = path.m_pathVertices[i + 1];
        (*active_edges)[mp_problem->m_mapEdge_link2index.at(
            std::make_pair(source, target))] = true;
      }
    }
  }

  return found_solution;
}

bool HESolver::checkVertexIsInAPath(
    const VertexIndex& v, EdgeLink* outIndices) const {
  assert(outIndices != nullptr);

  // First of all check if v is linked
  if (m_linkedVertices.find(v) != m_linkedVertices.cend()) {
    for (unsigned i = 0; i < m_foundPaths.size(); ++i) {
      const Path& p = m_foundPaths[i];
      for (unsigned j = 0; j < p.m_pathVertices.size(); ++j) {
        if (p.m_pathVertices[j] == v) {
          *outIndices = std::make_pair(i, j);
          return true;
        }
      }
    }
  }
  return false;
}

void HESolver::apply_next_node_to_path(const VertexIndex& v,
                                       Path* path) {
  assert(path != nullptr);

  // Assert the node you want to add is not already in the path
  assert(std::find(path->m_pathVertices.cbegin(),
                   path->m_pathVertices.cend(),
                   v) == path->m_pathVertices.cend());

  // Takes the distance from the last in the current path to v
  const RealNumber& distance_last_to_i =
      mp_problem->m_distance_matrix[path->m_pathVertices.back()][v];

  // Add v as last vertices in the path
  path->m_pathVertices.push_back(v);

  // Now you have to add the walk distance to all nodes in the
  // path
  for (unsigned i = 0;
       i < path->m_routeLenPerLevel.size();
       ++i) {
    path->m_routeLenPerLevel[i] += distance_last_to_i;

    // Assertion the new insertion does not compromize the
    // integrity of the path for no nodes
    assert(path->m_routeLenPerLevel[i] <=
           mp_problem->m_max_distances[path->m_pathVertices[i]]);
  }

  // For the last added in the path (v) the walk distance is zeros
  path->m_routeLenPerLevel.push_back(0);
}

void HESolver::find_possible_next_nodes(const Path& current_path,
                                        std::vector<VertexIndex>* nexts) {
  /// This algorithm give a node should return all possibile next
  /// nodes.
  assert(nexts != nullptr);

  const unsigned num_nodes = mp_problem->m_numNodes;

  // Clean the output
  nexts->clear();

  const unsigned LEN_PATH = current_path.m_pathVertices.size();

  // If current_path is empty, then you can select all except the school
  if (LEN_PATH == 0) {
    for (VertexIndex i = 0; i < num_nodes; ++i) {
      if (i != SCHOOL_INDEX) {
        nexts->push_back(i);
      }
    }
  } else {
    /// In that case there is path

    // Get the last of the path
    const auto& last_node = current_path.m_pathVertices.back();

    for (VertexIndex i = 0; i < num_nodes; ++i) {
      // The i node cannot be already in the path
      if (std::find(current_path.m_pathVertices.cbegin(),
                    current_path.m_pathVertices.cend(),
                    i) == current_path.m_pathVertices.cend()) {
        const RealNumber& distance_last_to_i =
            mp_problem->m_distance_matrix[last_node][i];

        // Check if adding the i node no vertices is not consistent
        bool consistent = true;
        for (unsigned n = 0; n < LEN_PATH && consistent == true; ++n) {
          const auto& level_node = current_path.m_pathVertices[n];
          if (current_path.m_routeLenPerLevel[n] +
              distance_last_to_i > mp_problem->m_max_distances[level_node]) {
            consistent = false;
          }
        }
        if (consistent == true) {
          nexts->push_back(i);
        }
      }
    }
  }
}

}  // namespace for_ch
