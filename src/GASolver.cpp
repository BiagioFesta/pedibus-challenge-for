// Copyright 2017 <Biagio Festa>
#include <memory>
#include <cassert>
#include "GASolver.hpp"

namespace for_ch {

const GASolver* GASolver::mps_running_solver = nullptr;

GASolver::GASolver(const std::shared_ptr<ProblemDatas>& problem) noexcept :
    mp_problem(problem) {
}

void GASolver::run(int argc, char** argv) {
  assert(argv != nullptr);

  // Set global pointer as this solver
  mps_running_solver = this;

  // Get the numer of edges
  const unsigned num_edges = mp_problem->m_numEdges;

  // Initialize a genome
  Genome genome(num_edges, &GASolver::ga_genome_fitness);
  genome.initializer(&GASolver::ga_genome_init);
  genome.mutator(&GASolver::ga_genome_mutator);

  // Initialize genetic algorithm
  GASimpleGA ga(genome);
  ga.set(gaNpopulationSize, 10);
  ga.set(gaNpCrossover, 1);
  ga.set(gaNpMutation, 0.05);
  ga.set(gaNnGenerations, 100000);
  ga.minimize();
  // ga.terminator(&GASolver::ga_algorithm_terminator);
  // ga.parameters(argc, argv);

  // Run algorithm
  if (m_displayInfo == true) {
    while (ga.done() != gaTrue) {
      print_current_ga_state(ga, &std::cout);
      ga.step();
    }
  } else {
    ga.evolve();
  }

  const auto& result_stats = ga.statistics();
  const Genome& best_result = (const Genome&) result_stats.bestIndividual();
  std::cout << best_result.score() << std::endl;
}

void GASolver::print_current_ga_state(const GAGeneticAlgorithm& ga,
                                      std::ostream* os) const noexcept {
  assert(os != nullptr);
}

void GASolver::ga_genome_init(GAGenome& g) noexcept {
  Genome& genome = (Genome&) g;
  mps_running_solver->init_genome_w_trivial_solution<Genome>(&genome);
}
float GASolver::ga_genome_fitness(GAGenome& g) noexcept {
  Genome& genome = (Genome&) g;
  return mps_running_solver->evaluator_genome<Genome>(genome);
}
int GASolver::ga_genome_mutator(GAGenome& g, float mp) noexcept {
  Genome& genome = (Genome&) g;
  return mps_running_solver->mutator_genome<Genome>(&genome, mp);
}
GABoolean GASolver::ga_algorithm_terminator(GAGeneticAlgorithm & ga) noexcept {
  return gaFalse;
}

}  // namespace for_ch
