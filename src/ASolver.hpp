// Copyright 2017 <Biagio Festa> <Alessandro Erba>
#ifndef __FOR_CH_A_SOLVER_HPP
#define __FOR_CH_A_SOLVER_HPP
#include <random>
#include <set>
#include <vector>
#include <algorithm>
#include <queue>
#include "ProblemDatas.hpp"

namespace for_ch {

class ASolver {
 public:
  ASolver(const ProblemDatas* problem);

  /// @brief run the solver
  /// @param [out] active_edges  The active active in the solution found
  /// @return the number of leaves found in that solution
  /// @note The algorithm is randomic, that means each run may find different
  /// solution
  void run(Solution* out_soution);

 private:
  struct Path {
    Path(const ProblemDatas* problem);
    bool add_node(VertexIndex node);
    
    std::vector<VertexIndex> m_nodes;
    std::vector<RealNumber> m_walked_distances;
    const ProblemDatas* mp_problem;
  };

  struct CompareDistanceWithSchool {
    CompareDistanceWithSchool(const ProblemDatas* problem);
    bool operator()(const VertexIndex& v1,
                    const VertexIndex& v2) const {
      return mp_problem->m_distance_matrix.at(v1).at(SCHOOL_INDEX) >
        mp_problem->m_distance_matrix.at(v2).at(SCHOOL_INDEX);
    }
    const ProblemDatas* mp_problem;
  };

  // Problem Datas
  const ProblemDatas* mp_problem;

  // A set of node which have to be linked
  // ordered in according to the distance with the school
  std::priority_queue<VertexIndex,
                      std::vector<VertexIndex>,
                      CompareDistanceWithSchool> m_freeNodes;

  // Paths found
  std::vector<Path> m_found_paths;

  // Cache vector for selecting the node to link
  std::vector<VertexIndex> m_nodes_discarded;
  
  // Random Engine Generator
  std::mt19937_64 m_rnd_engine;
  
  // Uniform Random variable (0, 1)
  std::uniform_real_distribution<RealNumber> m_rnd_variable;

  // Probability to swap in perturbation
  RealNumber m_inital_prob_swap = 0.2;
  
  // Scale factor of probability perturbation
  RealNumber m_prob_scale = 0.8;

  void buildPath();
};

}  // namespace for_ch

#endif  // __FOR_CH_A_SOLVER_HPP
