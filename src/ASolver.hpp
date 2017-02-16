// Copyright 2017 <Biagio Festa> <Alessandro Erba>
#ifndef __FOR_CH_A_SOLVER_HPP
#define __FOR_CH_A_SOLVER_HPP
#include "ProblemDatas.hpp"

namespace for_ch {

class ASolver {
 public:
  ASolver(const ProblemDatas* problem);
  void run(std::vector<bool>* active_edges);

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
      return mp_problem->m_distance_matrix.at(v1).at(SCHOOL_INDEX) <
        mp_problem->m_distance_matrix.at(v2).at(SCHOOL_INDEX);
    }
    const ProblemDatas* mp_problem;
  };

  // Problem Datas
  const ProblemDatas* mp_problem;

  // A set of node which have to be linked
  // ordered in according to the distance with the school
  std::set<VertexIndex, CompareDistanceWithSchool> m_freeNodes;

  // Paths found
  std::vector<Path> m_found_paths;

  void buildPath();
};

}  // namespace for_ch

#endif  // __FOR_CH_A_SOLVER_HPP
