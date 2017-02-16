// Copyright 2017 <Biagio Festa> <Alessandro Erba>
#include <vector>
#include <algorithm>
#include <utility>
#include <random>
#include <chrono>
#include "ASolver.hpp"

namespace for_ch {
ASolver::ASolver(const ProblemDatas* problem) :
    mp_problem(problem),
    m_freeNodes(CompareDistanceWithSchool(problem)) {
}

unsigned ASolver::run(std::vector<bool>* active_edges) {
  assert(active_edges != nullptr);

  // Clean state
  m_freeNodes.clear();
  m_found_paths.clear();

  // Set seed
  m_rnd_engine.seed(
      std::chrono::system_clock::now().time_since_epoch().count());

  const unsigned NUM_NODES = mp_problem->m_numNodes;
  for (VertexIndex i = 0; i < NUM_NODES; ++i) {
    if (i != SCHOOL_INDEX) {
      m_freeNodes.insert(i);
    }
  }

  while (m_freeNodes.empty() == false) {
    buildPath();
  }

  active_edges->resize(mp_problem->m_numEdges, false);
  for (const auto& path : m_found_paths) {
    const std::vector<VertexIndex> nodes = path.m_nodes;
    const auto it_last_no_school = nodes.crend() - 1;
    for (auto it = nodes.crbegin(); it != it_last_no_school; ++it) {
      VertexIndex source = *it;
      VertexIndex target = *(it + 1);
      assert(mp_problem->m_mapEdge_link2index.find(
          std::make_pair(source, target)) !=
             mp_problem->m_mapEdge_link2index.cend());
      EdgeIndex edge_index =
          mp_problem->m_mapEdge_link2index.at(std::make_pair(source, target));
      (*active_edges)[edge_index] = true;
    }
  }

  return m_found_paths.size();
}

void ASolver::buildPath() {
  // Get the first node to link
  VertexIndex node_to_link = *m_freeNodes.cbegin();

  // Sort the path in accordance with the distance
  std::sort(m_found_paths.begin(), m_found_paths.end(),
            [&node_to_link, this] (const Path& p1, const Path& p2) -> bool {
              VertexIndex l1 = p1.m_nodes.back();
              VertexIndex l2 = p2.m_nodes.back();
              float d1 = mp_problem->m_distance_matrix.at(node_to_link).at(l1);
              float d2 = mp_problem->m_distance_matrix.at(node_to_link).at(l2);
              return d1 < d2;
            });

  const unsigned N_PATH_FOUND = m_found_paths.size();

  // Add randomess perturbations
  if (N_PATH_FOUND > 1) {
    RealNumber prob_swap = m_inital_prob_swap;
    for (unsigned i = 0; i < N_PATH_FOUND - 1; ++i) {
      const auto x = m_rnd_variable(m_rnd_engine);
      if (x <= prob_swap) {
        std::swap(m_found_paths[i], m_found_paths[i+1]);
      }
      prob_swap *= m_prob_scale;
    }
  }

  // Now we try to link the node with the other path
  bool found = false;
  for (unsigned i = 0; i < N_PATH_FOUND && found == false; ++i) {
    // Get the i-th path and try to attach the node to link
    found = m_found_paths[i].add_node(node_to_link);
  }

  if (found == false) {
    Path temp_new(mp_problem);
    bool added = temp_new.add_node(node_to_link);
    assert(added == true);

    m_found_paths.push_back(std::move(temp_new));
  }

  m_freeNodes.erase(node_to_link);
}

ASolver::Path::Path(const ProblemDatas* problem) :
    m_nodes{SCHOOL_INDEX},
    m_walked_distances{0},
    mp_problem(problem) {
}

bool ASolver::Path::add_node(VertexIndex node) {
  RealNumber maxWalkableDistance = mp_problem->m_max_distances[node];
  RealNumber lastWalkedDistance = m_walked_distances.back();
  VertexIndex lastNode = m_nodes.back();
  RealNumber newWalkedDistance = lastWalkedDistance +
      mp_problem->m_distance_matrix[lastNode][node];
  if (newWalkedDistance <= maxWalkableDistance) {
    m_nodes.push_back(node);
    m_walked_distances.push_back(newWalkedDistance);
    return true;
  }
  return false;
}

ASolver::CompareDistanceWithSchool::CompareDistanceWithSchool(
    const ProblemDatas* problem): mp_problem(problem) {
}

}  // namespace for_ch
