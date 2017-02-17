// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH__GA_SOLVER__HPP
#define __FOR_CH__GA_SOLVER__HPP
#include <ga/ga.h>
#include <memory>
#include <algorithm>
#include <chrono>
#include <vector>
#include <set>
#include "ProblemDatas.hpp"

namespace for_ch {

class GASolver {
 public:
  explicit GASolver(const ProblemDatas& problem) noexcept;

  void set_flag_custom_crossover(bool flag) noexcept {
    m_custom_crossover = flag;
  }

  void set_max_time_seconds(int seconds) noexcept {
    m_timeMax_seconds = seconds;
  }

  void set_verbose(bool flag) noexcept {
    m_verbose = flag;
  }

  Solution run(int argc,
               char** argv,
               const std::set<Solution>& initial_solutions);
 private:
  using Genome = GA1DBinaryStringGenome;
  using Clock = std::chrono::system_clock;
  using TimePoint = Clock::time_point;

  enum FeasibilityStatus {
    FEASIBLE,
    NOTF_CYCLES,
    NOTF_DOUBLE_EDGE,
    NOTF_DISTANCES_CONST,
    NOTF_NOT_REACH_SCHOOL
  };

  struct ComparatorSolutionNumLeaves {
    ComparatorSolutionNumLeaves() = default;
    bool operator()(const Solution& s1,
                    const Solution& s2) const noexcept {
      return s1.m_num_leaves > s2.m_num_leaves;
    }
  };

  const ProblemDatas* mp_problem;
  static const GASolver* mps_running_solver;
  TimePoint m_time_start;
  float m_fitness_beta = 0.1f;
  bool m_display_info = true;

  /// parameters genetic algorithm
  float m_pCrossover;
  float m_pMutation;
  unsigned m_sizePopulation;
  unsigned m_timeMax_seconds;
  bool m_verbose = false;
  bool m_custom_crossover = true;

  // Set the paramters in accordance with the size problem
  void set_default_parameter() noexcept;

  FeasibilityStatus check_feasibility(
    const Genome& genome,
    unsigned* num_leaves = nullptr,
    RealNumber* tot_danger = nullptr,
    std::vector<RealNumber>* out_walked_distances = nullptr) const;

  void init_genome_with_solution(const Solution& solution,
                                 Genome* out_genome) const;

  static float EvaluatorGenome(GAGenome& genome);
  static int MutatorGenome(GAGenome& genome, float pmut);
  static int CrossoverSexGenome(const GAGenome& dad,
                                const GAGenome& mom,
                                GAGenome* bro,
                                GAGenome* sis);
  static GABoolean TerminatorGa(GAGeneticAlgorithm& ga);
  static void PrintCurrentState(const GAGeneticAlgorithm& ga,
                                std::ostream* os);

  template<typename T>
  static void AddSet(const std::set<T>& source,
                     std::set<T>* destination) noexcept;

};  // class GASolver

template<typename T>
void GASolver::AddSet(const std::set<T>& source,
                      std::set<T>* destination) noexcept {
  std::for_each(source.cbegin(),
                source.cend(),
                [&destination] (const T& e) {
                  destination->insert(e);
                });
}

}  // namespace for_ch

#endif  // __FOR_CH__GA_SOLVER__HPP
