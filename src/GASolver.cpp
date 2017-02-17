// Copyright 2017 <Biagio Festa>
#include <memory>
#include <cassert>
#include <vector>
#include <utility>
#include <set>
#include <unordered_map>
#include "GASolver.hpp"

namespace for_ch {

const GASolver* GASolver::mps_running_solver = nullptr;

GASolver::GASolver(const ProblemDatas& problem) noexcept :
    mp_problem(&problem) {
  set_default_parameter();
}

Solution GASolver::run(
    int argc, char** argv, const std::set<Solution>& initial_solutions) {
  assert(argv != nullptr);

  // Get some information about the problem
  const unsigned num_edges = mp_problem->m_numEdges;

  // Set global pointer at this solver
  mps_running_solver = this;

  // Initalize Random seeds
  GARandomSeed();

  // Initialize the first population
  GAPopulation initial_population;
  initial_population.order(GAPopulation::LOW_IS_BEST);
  for (const auto& solution : initial_solutions) {
    Genome genome(num_edges, &EvaluatorGenome);
    genome.mutator(&MutatorGenome);
    if (m_custom_crossover == true) {
      genome.crossover(&CrossoverSexGenome);
    }
    init_genome_with_solution(solution, &genome);
    assert(check_feasibility(genome) == FEASIBLE);
    initial_population.add(genome);
  }

  // Initialize genetic algorithm
  GASimpleGA ga(initial_population);
  // ga.set(gaNpopulationSize, m_sizePopulation);
  ga.set(gaNpCrossover, m_pCrossover);
  ga.set(gaNpMutation, m_pMutation);
  ga.minimize();
  // ga.initialize();
  ga.terminator(&TerminatorGa);
  ga.parameters(argc, argv);

  // Set the start time
  m_time_start = Clock::now();

  // Execute algorithm
  if (m_display_info == true) {
    unsigned index_step = 0;
    while (ga.done() != gaTrue) {
      ga.step();
      PrintCurrentState(ga, &std::cout);
      ++index_step;
    }
  } else {
    ga.evolve();
  }

  // Get the best genome and construct best solution
  const GAStatistics& stats = ga.statistics();
  const GAPopulation& best_pop = stats.bestPopulation();
  const Genome& best_genome = dynamic_cast<const Genome&>(best_pop.best());
  Solution solution;
  solution.m_active_edges.resize(best_genome.size());
  for (int b = 0; b < best_genome.size(); ++b) {
    solution.m_active_edges[b] = best_genome.gene(b);
  }
  assert(solution.compute_feasibility(*mp_problem) == true);
  unsigned best_num_leaves;
  check_feasibility(best_genome, &best_num_leaves,
                    &solution.m_danger);
  solution.m_num_leaves = best_num_leaves;
  return solution;
}

void GASolver::set_default_parameter() noexcept {
  const unsigned num_nodes = mp_problem->m_numNodes;
  if (num_nodes <= 11) {
    m_pCrossover = 0.9f;
    m_pMutation = 0.05f;
    m_sizePopulation = 10;
    m_custom_crossover = true;
  } else if (num_nodes <= 21) {
    m_pCrossover = 0.7f;
    m_pMutation = 0.002f;
    m_sizePopulation = 4;
    m_custom_crossover = true;
  } else if (num_nodes <= 31) {
    m_pCrossover = 0.99f;
    m_pMutation = 0.05f;
    m_sizePopulation = 12;
    m_custom_crossover = true;
  } else if (num_nodes <= 51) {
    m_pCrossover = 0.9f;
    m_pMutation = 0.05f;
    m_sizePopulation = 4;
    m_custom_crossover = true;
  } else if (num_nodes <= 81) {
    m_pCrossover = 0.7f;
    m_pMutation = 0.02f;
    m_sizePopulation = 6;
    m_custom_crossover = true;
  } else if (num_nodes <= 101) {
    m_pCrossover = 0.7f;
    m_pMutation = 0.02f;
    m_sizePopulation = 6;
    m_custom_crossover = true;
  } else if (num_nodes <= 151) {
    m_pCrossover = 0.7f;
    m_pMutation = 0.04f;
    m_sizePopulation = 6;
    m_custom_crossover = true;
  } else if (num_nodes <= 201) {
    m_pCrossover = 0.7f;
    m_pMutation = 0.04f;
    m_sizePopulation = 6;
    m_custom_crossover = true;
  } else if (num_nodes <= 251) {
    m_pCrossover = 0.7f;
    m_pMutation = 0.02f;
    m_sizePopulation = 12;
    m_custom_crossover = true;
  } else if (num_nodes <= 301) {
    m_pCrossover = 0.9f;
    m_pMutation = 0.01f;
    m_sizePopulation = 8;
    m_custom_crossover = true;
  } else {
    assert(false);  // you should not reach this point
  }

  // Set the fitness beta
  // TODO(biagio): check range
  if (num_nodes > 11 && num_nodes <= 101) {
    m_fitness_beta = 0.01;
  } else if (num_nodes > 101 && num_nodes <= 1001) {
    m_fitness_beta = 0.001;
  } else if (num_nodes > 1001) {
    m_fitness_beta = 0.0001;
  }
}

void GASolver::init_genome_with_solution(const Solution& solution,
                                         Genome* out_genome) const {
  const unsigned num_edges = mp_problem->m_numEdges;
  for (EdgeIndex e = 0; e < num_edges; ++e) {
    assert(e < solution.m_active_edges.size());
    assert(static_cast<int>(e) < out_genome->size());
    out_genome->gene(e, solution.m_active_edges[e]);
  }
}

GASolver::FeasibilityStatus GASolver::check_feasibility(
    const Genome& genome,
    unsigned* num_leaves,
    RealNumber* tot_danger,
    std::vector<RealNumber>* out_walked_distances) const {
  assert(mp_problem != nullptr);
  assert(SCHOOL_INDEX == 0);

  const unsigned num_nodes = mp_problem->m_numNodes;
  const unsigned num_nodes_minusone = num_nodes - 1;

  std::set<VertexIndex> inner_vertices;
  std::vector<std::set<VertexIndex>> nexts(num_nodes);
  std::vector<std::set<VertexIndex>> prevs(num_nodes);
  std::vector<RealNumber> walked_distance(num_nodes, 0);
  RealNumber danger = 0;
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
      assert(static_cast<int>(index_e) < genome.size());
      if (genome.gene(index_e) == 1) {
        if (has_outgoing_edge == true) {
          return FeasibilityStatus::NOTF_DOUBLE_EDGE;
        }
        has_outgoing_edge = true;

        // Since this is the first outgoing edge the walked distance is 0
        assert(walked_distance[n] == 0);

        // Get the target of that edge
        target = mp_problem->m_mapEdge_index2link.at(index_e).second;

        // Insert target as inner
        inner_vertices.insert(target);

        // Compute the delta distance
        delta_d = mp_problem->m_distance_matrix[n][target] +
            walked_distance[target];

        // Increase the total danger
        danger += mp_problem->m_dangerousness_matrix[n][target];

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

  // Set output if any
  if (num_leaves) {
    *num_leaves = num_nodes - inner_vertices.size();
  }
  if (tot_danger) {
    *tot_danger = danger;
  }

  // Check of distance constraint
  assert(walked_distance[SCHOOL_INDEX] == 0);
  for (VertexIndex n = 1; n < num_nodes; ++n) {
    if (walked_distance[n] > mp_problem->m_max_distances[n]) {
      if (out_walked_distances) {
        *out_walked_distances = std::move(walked_distance);
      }
      return NOTF_DISTANCES_CONST;
    }
  }

  // Check for cycles and school reachability
  for (VertexIndex n = 1; n < num_nodes; ++n) {
    if (nexts[n].find(n) != nexts[n].cend()) {
      if (out_walked_distances) {
          *out_walked_distances = std::move(walked_distance);
      }
      return NOTF_CYCLES;
    }

    if (nexts[n].find(SCHOOL_INDEX) == nexts[n].cend()) {
      if (out_walked_distances) {
          *out_walked_distances = std::move(walked_distance);
      }
      return NOTF_NOT_REACH_SCHOOL;
    }
  }

  if (out_walked_distances) {
    *out_walked_distances = std::move(walked_distance);
  }

  assert(prevs[SCHOOL_INDEX].size() == num_nodes - 1);
  return FEASIBLE;
}

float GASolver::EvaluatorGenome(GAGenome& genome) {
  const Genome& g = dynamic_cast<const Genome&>(genome);

  const unsigned num_nodes = mps_running_solver->mp_problem->m_numNodes;

  unsigned num_leaves;
  std::vector<RealNumber> walked_distance;
  RealNumber danger;

  FeasibilityStatus status =
      mps_running_solver->check_feasibility(g, &num_leaves,
                                            &danger, &walked_distance);

  switch (status) {
    case NOTF_CYCLES:
      return num_nodes * 10.f;
    case NOTF_DISTANCES_CONST:
    case NOTF_DOUBLE_EDGE:
    case NOTF_NOT_REACH_SCHOOL:
      return num_leaves + num_nodes;
    case FEASIBLE:
      float danger_scaled = danger * mps_running_solver->m_fitness_beta;
      return num_leaves + danger_scaled;
  }

  assert(false);  // you should not reach this point
  return num_nodes;
}

int GASolver::MutatorGenome(GAGenome& genome, float pmut) {
  Genome& g = dynamic_cast<Genome&>(genome);

  int num_mutation = 0;

  if (pmut > 0.f) {
    const unsigned num_nodes = mps_running_solver->mp_problem->m_numNodes;
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
          assert(static_cast<int>(index_e) < g.size());
          const auto gene = g.gene(index_e);

          if (gene == 0 && e == edge_to_activate) {
            // The gene is 0 but we want to activate it
            ++num_mutation;
            g.gene(index_e, 1);
          } else if (gene == 1 && e != edge_to_activate) {
            // The gene is 1 but we want to deactivate it
            ++num_mutation;
            g.gene(index_e, 0);
          }
        }  // for each edge of that node
      }  // if mutation to activate
    }  // for each node except school
  }  // if probability mutation > 0

  return num_mutation;
}

int GASolver::CrossoverSexGenome(const GAGenome& dad,
                                 const GAGenome& mom,
                                 GAGenome* bro,
                                 GAGenome* sis) {
  const Genome& g_dad = dynamic_cast<const Genome&>(dad);
  const Genome& g_mom = dynamic_cast<const Genome&>(mom);
  Genome* g_bro = dynamic_cast<Genome*>(bro);
  Genome* g_sis = dynamic_cast<Genome*>(sis);

  const unsigned num_nodes = mps_running_solver->mp_problem->m_numNodes;
  const unsigned num_nodes_minusone = num_nodes - 1;
  assert(num_nodes > 2);

  const int cross = GARandomInt(2, num_nodes_minusone);
  const EdgeIndex cross_edge = num_nodes_minusone * (cross - 1);
  const unsigned size_genome = g_dad.size();

  // Same len
  assert(g_dad.size() == g_mom.size());

  // *g_bro = g_dad;
  // *g_sis = g_mom;

  if (g_bro != nullptr) {
    g_bro->copy(g_dad, 0, 0, cross_edge);
    g_bro->copy(g_mom, cross_edge, cross_edge, size_genome - cross_edge);
  }

  if (g_sis != nullptr) {
    g_sis->copy(g_mom, 0, 0, cross_edge);
    g_sis->copy(g_dad, cross_edge, cross_edge, size_genome - cross_edge);
  }

  return 2;
}

GABoolean GASolver::TerminatorGa(GAGeneticAlgorithm& ga) {
  using Resolution = std::chrono::seconds;

  const auto time_now = Clock::now();
  const auto time_elapsed =
      time_now - mps_running_solver->m_time_start;

  const unsigned elapsed_count =
      std::chrono::duration_cast<Resolution>(time_elapsed).count();

  if (elapsed_count > mps_running_solver->m_timeMax_seconds) {
    return gaTrue;
  }

  return gaFalse;
}

void GASolver::PrintCurrentState(const GAGeneticAlgorithm& ga,
                                 std::ostream* os) {
  const GAStatistics& stats = ga.statistics();
  const GAPopulation& best_pop = stats.bestPopulation();
  const GAGenome& best_genome = best_pop.best();
  float best_score = best_genome.score();
  *os << best_score << std::endl;
}

}  // namespace for_ch
