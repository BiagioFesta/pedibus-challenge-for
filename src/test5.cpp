// Copyright 2017 <Biagio Festa>
#include <vector>
#include "ProblemDatas.hpp"
#include "ProblemState.hpp"

int main(int argc, char *argv[]) {
  ProblemDatas problem;
  problem.parse_problem_dat("pedibus_300.dat");
  ProblemState state(problem);
  std::vector<EdgeIndex> edges;
  state.compute_admissible_edges_with_source(1, &edges);

  state.add_admissible_edge(65);  // 7 --> 5
  state.compute_admissible_edges_with_source(5, &edges);

  return 0;
}
