// Copyright 2017 <Biagio Festa>
#include <iostream>
#include <memory>
#include "ProblemDatas.hpp"
#include "GASolver.hpp"

int main(int argc, char *argv[]) {
  auto ptr_problem = std::make_shared<for_ch::ProblemDatas>();
  ptr_problem->parse_problem_dat("pedibus_100.dat");

  for_ch::GASolver solver(ptr_problem);
  solver.run(argc, argv);

  return 0;
}
