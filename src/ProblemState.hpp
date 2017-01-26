// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH_PROBLEM_STATE__HPP
#define __FOR_CH_PROBLEM_STATE__HPP
#include <iterator>
#include "ProblemDatas.hpp"

class ProblemState {
 public:
  /// @brief Consructor
  /// @complexity O(|V|)
  explicit ProblemState(const ProblemDatas& problem) noexcept;

  /// @brief Compute all possibile next states from the current state
  /// @param [out] out_states  A container of state
  /// @note The container will not be cleaned! New states will be added
  /// to the container
  template<typename StateContainer>
  void compute_all_possibile_admissible_states(
      StateContainer* out_states) const noexcept;

  /// @brief Compute next state with the outgoing arc from v
  /// starting from the current state
  /// @param [in] v            A vertex you want to explode
  /// @param [out] out_states  A container of state
  /// @note The container will not be cleaned! New states will be added
  /// to the container
  template<typename StateContainer>
  void compute_admissible_states_with_source(
      const VertexIndex& v, StateContainer* out_states) const noexcept;

  template<typename EdgesContainer>
  void compute_all_possibile_admissible_edges(
      EdgesContainer* out_edges) const noexcept;

  template<typename EdgesContainer>
  void compute_admissible_edges_with_source(
      const VertexIndex& v, EdgesContainer* out_edges) const noexcept;

  /// @param [in] e    The edge you want to add to the state
  /// @return 'true' if the edge has been correcly added. 'false' in case that
  /// edge was already active
  bool add_admissible_edge(const EdgeIndex& e);

  /// @complexity O(1)
  /// @note complete does not implies admissible
  inline bool is_complete() const noexcept;

  /// @complexity O(|V|)
  /// @note admissibile not implies complete
  inline bool is_admissible() const noexcept;

 public:
  const ProblemDatas* mp_problem;

  /// @note size of this vector should always be num_edges of problem
  std::vector<bool> m_active_edges;

  /// The set of vertices which have been linked in a path
  std::set<VertexIndex> m_linked_vertices;
  //void compute_linked_vertices_set();

  /// @note size of this vector should always be num_nodes of problem
  std::vector<RealNumber> m_walk_distance_per_vertex;

  /// @note size of this vector should always be num_nodes of problem
  std::vector<std::vector<EdgeIndex>> m_out_edges_active_per_vertex;

  /// @note size of this vector should always be num_nodes of problem
  std::vector<std::vector<EdgeIndex>> m_in_edges_active_per_vertex;

  /// @brief give a vertex the set of all other states it can reach
  /// @note size of this vector should always be num_nodes of problem
  std::vector<std::set<VertexIndex>> m_reachableVertices_per_vertex;

  /// @return true if adding e the state remains admissible
  bool check_admissible_adding_edge(const EdgeIndex& e) const noexcept;

  /// @brief this function go down from v to the end of path
  /// starting from v and update distances.
  ///
  /// @note v will be considered as root. Moreover if going down
  /// the path of v the path is fork there only the path of
  /// v will be updated.
  void update_walk_distance_from_vertex(const VertexIndex& v);

  inline void find_roots_path_with(const VertexIndex& v,
                                   std::vector<VertexIndex>* roots) const noexcept;

  bool compute_distances_starting_v(
      const VertexIndex& v,
      const EdgeIndex* added_edge,
      std::vector<std::pair<VertexIndex, RealNumber>>* newDistances) const noexcept;
};

inline bool ProblemState::is_complete() const noexcept {
  assert(mp_problem != nullptr);
  // State is complete if all vertices are linked
  return m_linked_vertices.size() == mp_problem->m_numNodes;
}

inline bool ProblemState::is_admissible() const noexcept {
  assert(mp_problem != nullptr);
  const auto& numNodes = mp_problem->m_numNodes;
  assert(m_walk_distance_per_vertex.size() == numNodes);
  assert(mp_problem->m_max_distances.size() == numNodes);

  // If a vertex is walking more than its max then the solution is
  // not admissible or if there is more than one outedge from the same node
  // or there is a cycle
  for (unsigned i = 0; i < numNodes; ++i) {  // O(|V|)
    const std::set<VertexIndex>& reach_vertices =
        m_reachableVertices_per_vertex[i];
    if ((m_walk_distance_per_vertex[i] > mp_problem->m_max_distances[i]) ||
        (m_out_edges_active_per_vertex[i].size() > 1) ||
        (reach_vertices.find(i) != reach_vertices.cend())) {
      return false;
    }
  }
  return true;
}

inline void ProblemState::find_roots_path_with(
    const VertexIndex& v, std::vector<VertexIndex>* roots) const noexcept {
  assert(roots != nullptr);
  assert(v < m_out_edges_active_per_vertex.size());

  const std::vector<EdgeIndex>& out_edges =
      m_in_edges_active_per_vertex[v];

  if (out_edges.size() == 0) {
    roots->push_back(v);
  } else {
    for (const auto& edge : out_edges) {
      const VertexIndex& source = mp_problem->m_mapEdge_index2link.at(edge).first;
      find_roots_path_with(source, roots);
    }
  }
}

template<typename StateContainer>
void ProblemState::compute_all_possibile_admissible_states(
    StateContainer* out_states) const noexcept {
  assert(out_states != nullptr);

  // Compute all admissible edges from this state
  std::vector<EdgeIndex> admissible_edges;
  compute_all_possibile_admissible_edges(&admissible_edges);

  auto insert = std::inserter(*out_states, out_states->end());
  ProblemState new_state(*mp_problem);

  // Insert all admissible edges as new states
  for (const auto& e : admissible_edges) {
    new_state = *this;
    bool added = new_state.add_admissible_edge(e);
    assert(added == true);
    insert = std::move(new_state);
  }
}

template<typename StateContainer>
void ProblemState::compute_admissible_states_with_source(
    const VertexIndex& v, StateContainer* out_states) const noexcept {
    assert(out_states != nullptr);

  // Compute all admissible edges from this state
  std::vector<EdgeIndex> admissible_edges;
  compute_admissible_edges_with_source(v, &admissible_edges);

  auto insert = std::inserter(*out_states, out_states->end());
  ProblemState new_state(*mp_problem);

  // Insert all admissible edges as new states
  for (const auto& e : admissible_edges) {
    new_state = *this;
    bool added = new_state.add_admissible_edge(e);
    assert(added == true);
    insert = std::move(new_state);
  }
}

template<typename EdgesContainer>
void ProblemState::compute_all_possibile_admissible_edges(
    EdgesContainer* out_edges) const noexcept {
  assert(out_edges != nullptr);
  out_edges->clear();

  auto insert = std::inserter(*out_edges, out_edges->end());
  for (unsigned i = 0; i < m_active_edges.size(); ++i) {
    if (m_active_edges[i] == false) {
      if (check_admissible_adding_edge(i) == true) {
        insert = i;
      }
    }
  }
}

template<typename EdgesContainer>
void ProblemState::compute_admissible_edges_with_source(
    const VertexIndex& v, EdgesContainer* out_edges) const noexcept {
  assert(out_edges != nullptr);
  assert(v < m_out_edges_active_per_vertex.size());
  out_edges->clear();

  if (v != 0  && m_out_edges_active_per_vertex[v].size() == 0) {
    auto insert = std::inserter(*out_edges, out_edges->end());
    const EdgeIndex first_edge_of_v = (mp_problem->m_numNodes - 1) *
        (v - 1);
    const EdgeIndex last_edge_of_v = first_edge_of_v + mp_problem->m_numNodes - 1;

    for (auto i = first_edge_of_v; i < last_edge_of_v; ++i) {
      if (m_active_edges[i] == false) {
        if (check_admissible_adding_edge(i) == true) {
          insert = i;
        }
      }
    }
  }
}

#endif  // __FOR_CH_PROBLEM_STATE__HPP
