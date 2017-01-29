// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH__GA_SOLVER__HPP
#define __FOR_CH__GA_SOLVER__HPP
#undef GALIB_USE_BORLAND_INST
#include <ga/ga.h>
#include <memory>
#include "ProblemDatas.hpp"

namespace for_ch {

class GASolver {
 public:
  explicit GASolver(const std::shared_ptr<ProblemDatas>& problem) noexcept;

  void run(int argc, char** argv);
 private:
  using Genome = GA1DBinaryStringGenome;
  enum FeasibilityStatus {
    FEASIBLE,
    NOTF_CYCLES,
    NOTF_DOUBLE_EDGE,
    NOTF_DISTANCES_CONST,
    NOTF_NOT_REACH_SCHOOL
  };

  std::shared_ptr<ProblemDatas> mp_problem;
  static const GASolver* mps_running_solver;
  bool m_displayInfo = false;

  void print_current_ga_state(const GAGeneticAlgorithm& ga,
                              std::ostream* os) const noexcept;

  template<typename BinaryGenome>
  void init_genome_w_trivial_solution(BinaryGenome* genome) const noexcept;

  template<typename BinaryGenome>
  bool compute_walked_distances_from_end(
      const BinaryGenome& genome,
      const VertexIndex end_path,
      std::vector<RealNumber>* out_walked_distances) const noexcept;

  template<typename BinaryGenome>
  FeasibilityStatus is_feasible(
      const BinaryGenome& genome,
      unsigned* num_leaves,
      std::vector<RealNumber>* out_walked_distances) const noexcept;

  template<typename BinaryGenome>
  int mutator_genome(BinaryGenome* genome, float pmut) const noexcept;

  template<typename BinaryGenome>
  float evaluator_genome(const BinaryGenome& genome) const noexcept;

  static void ga_genome_init(GAGenome& g) noexcept;
  static float ga_genome_fitness(GAGenome& g) noexcept;
  static int ga_genome_mutator(GAGenome& g, float mp) noexcept;
  static GABoolean ga_algorithm_terminator(GAGeneticAlgorithm & ga) noexcept;
};  // class GASolver

template<typename BinaryGenome>
void GASolver::init_genome_w_trivial_solution(BinaryGenome* genome) const noexcept {
  assert(genome != nullptr);
  assert(mp_problem != nullptr);
  assert(SCHOOL_INDEX == 0);

  const unsigned num_nodes = mp_problem->m_numNodes;
  const unsigned num_nodes_minusone = num_nodes - 1;

  EdgeIndex index_e;

  // For each node except school
  for (unsigned n = 1; n < num_nodes; ++n) {
    // Compute starting edge index
    index_e = num_nodes_minusone * (n - 1);

    // For each edge of that node
    for (unsigned e = 0; e < num_nodes_minusone; ++e, ++index_e) {
      assert(static_cast<int>(index_e) < genome->size());

      if (e == 0) {  // e == 0 is the edge from node to scool
        genome->gene(index_e, 1);
      } else {
        genome->gene(index_e, 0);
      }
    }  // for each edge of that node
  }  // for each node except school
}

template<typename BinaryGenome>
bool GASolver::compute_walked_distances_from_end(
    const BinaryGenome& genome,
    const VertexIndex end_path,
    std::vector<RealNumber>* out_walked_distances) const noexcept {
  assert(out_walked_distances != nullptr);
  assert(mp_problem != nullptr);
  assert(SCHOOL_INDEX == 0);
  assert(out_walked_distances->size() == mp_problem->m_numNodes);

  const unsigned num_nodes = mp_problem->m_numNodes;
  const unsigned num_nodes_minusone = num_nodes - 1;

  EdgeIndex index_e;
  VertexIndex target;

  // For each node except school
  for (unsigned n = 1; n < num_nodes; ++n) {
    // Compute starting edge index
    index_e = num_nodes_minusone * (n - 1);

    // For each edge of that node
    for (unsigned e = 0; e < num_nodes_minusone; ++e, ++index_e) {
      assert(static_cast<int>(index_e) < genome.size());

      if (genome.gene(index_e) == 1) {
        // Get target of that edge
        target = mp_problem->m_mapEdge_index2link.at(index_e).second;

        if (target == end_path) {
          RealNumber& distance_n = (*out_walked_distances)[n];
          const RealNumber& distance_t = (*out_walked_distances)[target];

          distance_n = mp_problem->m_distance_matrix[n][target] +
              distance_t;

          if (distance_n > mp_problem->m_max_distances[n]) {
            return false;
          }

          if (compute_walked_distances_from_end(
                  genome, n, out_walked_distances) == false) {
            return false;
          }
        }
      }
    }  // for each edge of that node
  }  // for each node except school

  return true;
}

template<typename BinaryGenome>
GASolver::FeasibilityStatus GASolver::is_feasible(
    const BinaryGenome& genome,
    unsigned* num_leaves,
    std::vector<RealNumber>* out_walked_distances) const noexcept {
  assert(num_leaves != nullptr);
  assert(out_walked_distances != nullptr);
  assert(mp_problem != nullptr);
  assert(SCHOOL_INDEX == 0);
  assert(out_walked_distances->size() == mp_problem->m_numNodes);

  const unsigned num_nodes = mp_problem->m_numNodes;
  const unsigned num_nodes_minusone = num_nodes - 1;

  std::set<VertexIndex> inner_vertices;
  std::vector<std::set<VertexIndex>> reachable_from_vertex(num_nodes);
  bool has_outgoing_edge;
  EdgeIndex index_e;
  VertexIndex target;

  // For each node except school
  for (VertexIndex n = 1; n < num_nodes; ++n) {
    // Compute starting edge index
    index_e = num_nodes_minusone * (n - 1);

    has_outgoing_edge = false;

    // For each edge of that node
    for (unsigned e = 0; e < num_nodes_minusone; ++e, ++index_e) {
      assert(static_cast<int>(index_e) < genome.size());
      if (genome.gene(index_e) == 1) {
        if (has_outgoing_edge == true) {
          return FeasibilityStatus::NOTF_DOUBLE_EDGE;
        }

        has_outgoing_edge = true;

        // Get the target of that edge
        target = mp_problem->m_mapEdge_index2link.at(index_e).second;

        // Set the target vertex as inner
        inner_vertices.insert(target);

        // insert target as reachable from here
        reachable_from_vertex[n].insert(target);

        // Insert all reachable from target, reachable from here also
        std::for_each(reachable_from_vertex[target].cbegin(),
                      reachable_from_vertex[target].cend(),
                      [&reachable_from_vertex, &n] (const VertexIndex& v) {
                        reachable_from_vertex[n].insert(v);
                      });

        // Update all reachable list of those vertices which reach n
        for (VertexIndex v2 = 1; v2 < num_nodes; ++v2) {
          if (v2 != n) {
            std::set<VertexIndex>& reach_from_v2 = reachable_from_vertex[v2];
            if (reach_from_v2.find(n) != reach_from_v2.cend()) {
              std::for_each(reachable_from_vertex[n].cbegin(),
                            reachable_from_vertex[n].cend(),
                            [&reach_from_v2] (const VertexIndex& v) {
                              reach_from_v2.insert(v);
                            });
            }
          }
        }
      }  // if edge is active
    }  // for each edge of that node
  }  // for each node except school

  // Compute the number of leaves
  *num_leaves = num_nodes - inner_vertices.size();

  for (VertexIndex n = 1; n < num_nodes; ++n) {
    const std::set<VertexIndex>& reachset = reachable_from_vertex[n];

    // Check for cycles
    if (reachset.find(n) != reachset.cend()) {
      return FeasibilityStatus::NOTF_CYCLES;
    }

    // Check reachability to school
    if (reachset.find(SCHOOL_INDEX) == reachset.cend()) {
      return FeasibilityStatus::NOTF_NOT_REACH_SCHOOL;
    }
  }

  // Compute walk distance
  bool distance_constraint =
      compute_walked_distances_from_end(genome,
                                        SCHOOL_INDEX,
                                        out_walked_distances);

  // Check distances constraint
  if (distance_constraint == false) {
    return FeasibilityStatus::NOTF_DISTANCES_CONST;
  }

  return FeasibilityStatus::FEASIBLE;
}

template<typename BinaryGenome>
int GASolver::mutator_genome(BinaryGenome* genome, float pmut) const noexcept {
  assert(genome != nullptr);
  assert(mp_problem != nullptr);
  assert(SCHOOL_INDEX == 0);

  int num_mutation = 0;

  if (pmut > 0.f) {
    const unsigned num_nodes = mp_problem->m_numNodes;
    const unsigned num_nodes_minusone = num_nodes - 1;

    EdgeIndex index_e;
    EdgeIndex edge_to_activate;

    // For each node except school
    for (VertexIndex n = 1; n < num_nodes; ++n) {
      // Mutation active one random edge
      if (GAFlipCoin(pmut)) {
         // Compute starting edge index
        index_e = num_nodes_minusone * (n - 1);
        edge_to_activate = GARandomInt(0, num_nodes_minusone);

        // Now get off all edge except edge_to_active
        for (unsigned e = 0; e < num_nodes_minusone; ++e, ++index_e) {
          assert(static_cast<int>(index_e) < genome->size());
          const auto gene = genome->gene(index_e);

          if (gene == 0 && e == edge_to_activate) {
            // The gene is 0 but we want to activate it
            ++num_mutation;
            genome->gene(index_e, 1);
          } else if (gene == 1 && e != edge_to_activate) {
            // The gene is 1 but we want to deactivate it
            ++num_mutation;
            genome->gene(index_e, 0);
          }
        }  // for each edge of that node
      }  // if mutation to activate
    }  // for each node except school
  }  // if probability mutation > 0

  return num_mutation;
}

template<typename BinaryGenome>
float GASolver::evaluator_genome(const BinaryGenome& genome) const noexcept {
  // Get the number of nodes
  const unsigned num_nodes = mp_problem->m_numNodes;

  std::vector<RealNumber> walked_distance(num_nodes);
  unsigned num_leaves;
  float fitness;

  FeasibilityStatus status = is_feasible(genome, &num_leaves, &walked_distance);
  switch (status) {
    case NOTF_CYCLES:
      fitness = num_nodes * 10.f;
      return fitness;
    case NOTF_DISTANCES_CONST:
    case NOTF_DOUBLE_EDGE:
    case NOTF_NOT_REACH_SCHOOL:
      fitness = num_leaves + num_nodes;
      return fitness;
    case FEASIBLE:
      fitness = num_leaves;
      return fitness;
  }
  // This line should never be reached
  assert(false);
  return fitness;
}

}

#endif  // __FOR_CH__GA_SOLVER__HPP
