// Copyright 2017 <Biagio Festa>
#include <memory>
#include <cassert>
#include <vector>
#include "GASolver.hpp"

namespace for_ch {

const GASolver* GASolver::mps_running_solver = nullptr;
const std::vector<bool>* GASolver::mps_hint = nullptr;

GASolver::GASolver(const std::shared_ptr<ProblemDatas>& problem) noexcept :
    mp_problem(problem) {
}

void GASolver::run(int argc, char** argv, const std::vector<bool>* hint) {
  assert(argv != nullptr);

  // Set global pointer as this solver
  mps_running_solver = this;

  // Get the numer of edges
  const unsigned num_edges = mp_problem->m_numEdges;

  // Set hit if any
  mps_hint = (hint != nullptr ? hint : nullptr);

  // Initalize Random seeds
  GARandomSeed();

  // Initialize a genome
  Genome genome(num_edges, &GASolver::ga_genome_fitness);
  genome.initializer(&GASolver::ga_genome_init);
  genome.mutator(&GASolver::ga_genome_mutator);
  genome.crossover(&GASolver::ga_genome_crossover);

  // Initialize genetic algorithm
  GASimpleGA ga(genome);
  ga.set(gaNpopulationSize, m_sizePopulation);
  ga.set(gaNpCrossover, m_pCrossover);
  ga.set(gaNpMutation, m_pMutation);
  ga.set(gaNnGenerations, 100000);
  ga.minimize();
  ga.initialize();
  ga.terminator(&GASolver::ga_algorithm_terminator);

  // Save log on disk
  ga.scoreFilename("log_genetic.txt");
  ga.selectScores(GAStatistics::Minimum |
                  GAStatistics::Mean |
                  GAStatistics::Maximum |
                  GAStatistics::Deviation |
                  GAStatistics::Diversity);
  ga.scoreFrequency(1);
  ga.flushFrequency(10);

  // Override parameters from cmd
  ga.parameters(argc, argv);

  // Set the start time
  m_time_start = Clock::now();

  // Run algorithm
  print_ga_parameters(ga, &std::cout);

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
  if (mps_hint != nullptr) {
    mps_running_solver->init_genome_w_solution<Genome>(
        &genome, *mps_hint);
  } else {
    mps_running_solver->init_genome_w_trivial_solution<Genome>(&genome);
  }
}

float GASolver::ga_genome_fitness(GAGenome& g) noexcept {
  Genome& genome = (Genome&) g;
  return mps_running_solver->evaluator_genome<Genome>(genome);
}

int GASolver::ga_genome_mutator(GAGenome& g, float mp) noexcept {
  Genome& genome = (Genome&) g;
  return mps_running_solver->mutator_genome<Genome>(&genome, mp);
}

int GASolver::ga_genome_crossover(const GAGenome& dad,
                                  const GAGenome& mom,
                                  GAGenome* bro,
                                  GAGenome* sis) noexcept {
  assert(bro != nullptr);
  assert(sis != nullptr);
  Genome& genome_dad = (Genome&) dad;
  Genome& genome_mom = (Genome&) mom;
  Genome& genome_bro = dynamic_cast<Genome&>(*bro);
  Genome& genome_sis = dynamic_cast<Genome&>(*sis);

  return mps_running_solver->crossover_genome(genome_dad, genome_mom,
                                              &genome_bro, &genome_sis);
}

GABoolean GASolver::ga_algorithm_terminator(GAGeneticAlgorithm & ga) noexcept {
  using Resolution = std::chrono::seconds;

  const auto time_now = Clock::now();
  const auto time_elapsed =
      time_now - mps_running_solver->m_time_start;

  const unsigned elapsed_count =
      std::chrono::duration_cast<Resolution>(time_elapsed).count();

  if (elapsed_count > 60 * 10) {  // 10 minutes
    return gaTrue;
  }

  return gaFalse;
}

void GASolver::print_ga_parameters(const GAGeneticAlgorithm& ga,
                                  std::ostream* os) const noexcept {
  assert(os != nullptr);
  *os << "----------------------------------------\n"
      << "GENETIC ALGORITHM PARAMETERS:\n"
      << "Population Size: " << ga.populationSize() << "\n"
      << "Num. Generations: " << ga.nGenerations() << "\n"
      << "Num. Convergence: " << ga.nConvergence() << "\n"
      << "Prob. Convergence: " << ga.pConvergence() << "\n"
      << "Prob. Mutation: " << ga.pMutation() << "\n"
      << "Prob. Crossover: " << ga.pCrossover() << "\n"
      << "-----------------------------------------\n";
  os->flush();
}

}  // namespace for_ch
