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

void ASolver::run(Solution* out_solution) {
  assert(out_solution != nullptr);

  // Clean state
  while (m_freeNodes.empty() == false) {
    m_freeNodes.pop();
  }
  m_found_paths.clear();

  // Set seed
  m_rnd_engine.seed(
      std::chrono::system_clock::now().time_since_epoch().count());

  const unsigned NUM_NODES = mp_problem->m_numNodes;
  for (VertexIndex i = 0; i < NUM_NODES; ++i) {
    if (i != SCHOOL_INDEX) {
      m_freeNodes.push(i);
      assert(m_freeNodes.size() == i);
    }
  }
  assert(m_freeNodes.size() == NUM_NODES - 1);

  while (m_freeNodes.empty() == false) {
    buildPath();
  }

  out_solution->m_active_edges.clear();
  out_solution->m_active_edges.resize(mp_problem->m_numEdges, false);
  for (const auto& path : m_found_paths) {
    const std::vector<VertexIndex>& nodes = path.m_nodes;
    assert(nodes.size() > 1);
    const auto it_last_no_school = nodes.crend() - 1;
    for (auto it = nodes.crbegin(); it != it_last_no_school; ++it) {
      VertexIndex source = *it;
      VertexIndex target = *(it + 1);
      assert(mp_problem->m_mapEdge_link2index.find(
          std::make_pair(source, target)) !=
             mp_problem->m_mapEdge_link2index.cend());
      EdgeIndex edge_index =
          mp_problem->m_mapEdge_link2index.at(std::make_pair(source, target));
      assert(edge_index < out_solution->m_active_edges.size());
      out_solution->m_active_edges[edge_index] = true;
    }
  }

  out_solution->m_num_leaves = m_found_paths.size();
  out_solution->m_danger = -1;
  assert(out_solution->compute_feasibility(*mp_problem) == true);
}

void ASolver::buildPath() {
  // Compute the node to link (randomically)
  m_nodes_discarded.clear();
  const unsigned num_free_minone = m_freeNodes.size() - 1;
  VertexIndex node_to_link = SCHOOL_INDEX;
  RealNumber prob_swap = m_inital_prob_swap;
  for (unsigned i = 0; i < num_free_minone &&
           node_to_link == SCHOOL_INDEX; ++i) {
    const auto x = m_rnd_variable(m_rnd_engine);
    if (x <= prob_swap) {
      m_nodes_discarded.push_back(m_freeNodes.top());
      prob_swap *= m_prob_scale;
    } else {
      node_to_link = m_freeNodes.top();
    }
    m_freeNodes.pop();
  }
  if (node_to_link == SCHOOL_INDEX) {
    node_to_link = m_freeNodes.top();
    m_freeNodes.pop();
  }

  // Reinsert the discarded nodes
  for (const auto& n : m_nodes_discarded) {
    m_freeNodes.push(n);
  }

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
    prob_swap = m_inital_prob_swap;
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

  // If there is no other path to link, create new path
  if (found == false) {
    Path temp_new(mp_problem);
    bool added = temp_new.add_node(node_to_link);
    assert(added == true);

    m_found_paths.push_back(std::move(temp_new));
  }
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
