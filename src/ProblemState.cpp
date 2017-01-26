// Copyright 2017 <Biagio Festa>
#include "ProblemState.hpp"
#include <set>
#include <stack>
#include <vector>
#include <utility>

ProblemState::ProblemState(const ProblemDatas& problem) noexcept :
    mp_problem(&problem),
    m_active_edges(problem.m_numEdges, false),
    m_walk_distance_per_vertex(problem.m_numNodes, 0),
    m_out_edges_active_per_vertex(problem.m_numNodes),
    m_in_edges_active_per_vertex(problem.m_numNodes),
    m_reachableVertices_per_vertex(problem.m_numNodes) {
  for (unsigned i = 0; i < problem.m_numNodes; ++i) {  // O(|V|)
    m_out_edges_active_per_vertex[i].reserve(problem.m_numEdges);
    m_in_edges_active_per_vertex[i].reserve(problem.m_numEdges);
  }
}

bool ProblemState::compute_distances_starting_v(
    const VertexIndex& v,
    const EdgeIndex* added_edge,
    std::vector<std::pair<VertexIndex, RealNumber>>* newDistances)
    const noexcept {
  // TODO(biagio): this function could be recursive
  // maybe check for performances
  assert(mp_problem != nullptr);
  assert(m_out_edges_active_per_vertex.size() == mp_problem->m_numNodes);

  if (newDistances != nullptr) {
    newDistances->clear();
    newDistances->reserve(mp_problem->m_numNodes);
  }

  const auto& outedges_active = m_out_edges_active_per_vertex;

  // Added additional edge if any
  const VertexIndex* source;
  const VertexIndex* target;
  if (added_edge != nullptr) {
    source = &(mp_problem->m_mapEdge_index2link.at(*added_edge).first);
    target = &(mp_problem->m_mapEdge_index2link.at(*added_edge).second);
    assert(m_out_edges_active_per_vertex[*source].size() == 0);
  }


  std::stack<VertexIndex> path;
  std::stack<RealNumber> distances;
  path.push(v);
  distances.push(0);

  // Fill the stack. Roll the path
  const std::vector<EdgeIndex>* outedges_top =
      &(outedges_active[path.top()]);

  // Until the end of the path
  while (outedges_top->size() != 0 ||
         (added_edge != nullptr ? path.top() == *source : false)) {
    assert(outedges_top->size() == 1 ||
           (added_edge != nullptr &&
            path.top() == *source ?
            true : false));  // check there is only one out edge

    // Get the next step in the path
    const VertexIndex& next_node =
        ((added_edge != nullptr) && (path.top() == *source) ?
         *target :
         mp_problem->m_mapEdge_index2link.at((*outedges_top)[0]).second);

    // Push the internal distance
    distances.push(mp_problem->m_distance_matrix[path.top()][next_node]);

    // Push in the stack the next node in the path
    path.push(next_node);

    // get the out edges of new top of the stack
    outedges_top = &(outedges_active[path.top()]);
  }

  // Unroll the stack
  RealNumber distance_path = 0;
  bool out_of_distance = false;
  while (path.empty() == false) {
    const VertexIndex& step = path.top();

    if (newDistances != nullptr) {
      // Insert the new distance for node step
      newDistances->push_back(std::make_pair(step, distance_path));
    }

    if (distance_path != 0 &&
        distance_path >= mp_problem->m_max_distances[step]) {
      out_of_distance = true;
      if (newDistances == nullptr) {
        return true;  // If you don't care new distance return error
      }
    }

    // Increase distance path
    distance_path += distances.top();

    // Pop the stacks
    distances.pop();
    path.pop();
  }

  return out_of_distance;
}

void ProblemState::update_walk_distance_from_vertex(
    const VertexIndex& v) {
  assert(mp_problem != nullptr);
  assert(m_out_edges_active_per_vertex.size() == mp_problem->m_numNodes);

  std::vector<std::pair<VertexIndex, RealNumber>> newDistances;
  bool out_distance = compute_distances_starting_v(v, nullptr, &newDistances);
  assert(out_distance == false);

  for (const auto& d : newDistances) {
    m_walk_distance_per_vertex[d.first] = d.second;
  }

  /*
  std::stack<VertexIndex> path;
  std::stack<RealNumber> distances;
  path.push(v);
  distances.push(0);

  // Fill the stack. Roll the path
  const std::vector<EdgeIndex>* outedges_top =
      &(m_out_edges_active_per_vertex[path.top()]);

  // Until the end of the path
  while (outedges_top->size() != 0) {
    assert(outedges_top->size() == 1);  // check there is only one out edge

    // Get the next step in the path
    const VertexIndex& next_node =
        mp_problem->m_mapEdge_index2link.at((*outedges_top)[0]).second;

    // Push the internal distance
    distances.push(mp_problem->m_distance_matrix[path.top()][next_node]);

    // Push in the stack the next node in the path
    path.push(next_node);

    // get the out edges of new top of the stack
    outedges_top = &(m_out_edges_active_per_vertex[path.top()]);
  }

  // Unroll the stack
  RealNumber distance_path = 0;
  while (path.empty() == false) {
    const VertexIndex& step = path.top();
    assert(step < m_walk_distance_per_vertex.size());
    m_walk_distance_per_vertex[step] = distance_path;
    assert(m_walk_distance_per_vertex[step] <=
           mp_problem->m_max_distances[step]);
    distance_path += distances.top();

    // Pop the stacks
    distances.pop();
    path.pop();
  }
  */
}

bool ProblemState::add_admissible_edge(const EdgeIndex& e) {
  assert(mp_problem != nullptr);
  assert(m_active_edges.size() == mp_problem->m_numEdges);
  assert(e < mp_problem->m_numEdges);

  // If the edge is already active return false
  if (m_active_edges[e] == true) {
    return false;
  }

  const EdgeLink& link = mp_problem->m_mapEdge_index2link.at(e);
  const VertexIndex& source = link.first;
  const VertexIndex& target = link.second;

  // Active edge
  m_active_edges[e] = true;

  // Add vertices in linked list
  m_linked_vertices.insert(source);
  m_linked_vertices.insert(target);

  // Add out and in edge to source and target
  assert(source < m_out_edges_active_per_vertex.size());
  assert(target < m_in_edges_active_per_vertex.size());
  m_out_edges_active_per_vertex[source].push_back(e);
  m_in_edges_active_per_vertex[target].push_back(e);

  // Add reachable list of source
  std::set<VertexIndex>& reachable_vertices =
      m_reachableVertices_per_vertex[source];
  reachable_vertices.insert(target);
  // transitive property
  for (const auto& v : m_reachableVertices_per_vertex[target]) {
    assert(v != source);
    reachable_vertices.insert(v);
  }

  // Get roots of path
  std::vector<VertexIndex> roots_for_source;
  find_roots_path_with(source, &roots_for_source);

  // For all root update distances path's
  for (const auto& r : roots_for_source) {
    update_walk_distance_from_vertex(r);
  }

  assert(is_admissible() == true);

  return true;
}

bool ProblemState::check_admissible_adding_edge(
    const EdgeIndex& e) const noexcept {
  assert(mp_problem != nullptr);
  assert(m_active_edges[e] == false);  // don't want to consider already talken

  const VertexIndex& source = mp_problem->m_mapEdge_index2link.at(e).first;
  const VertexIndex& target = mp_problem->m_mapEdge_index2link.at(e).second;

  // Check the source does not already have a out edge
  assert(source < m_out_edges_active_per_vertex.size());
  if (m_out_edges_active_per_vertex[source].size() >= 1) {
    return false;
  }

  // Check if edge make loops
  assert(target < m_reachableVertices_per_vertex.size());
  if (m_reachableVertices_per_vertex[target].find(source) !=
      m_reachableVertices_per_vertex[target].cend()) {
    return false;
  }

  // Get roots of path
  std::vector<VertexIndex> roots_for_source;
  find_roots_path_with(source, &roots_for_source);

  for (const auto& r : roots_for_source) {
    if (compute_distances_starting_v(r, &e, nullptr) == true) {
      return false;
    }
  }

  return true;
}
