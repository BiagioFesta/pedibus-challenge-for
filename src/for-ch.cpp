// Copyright 2017 <Biagio Festa>
#include <tclap/CmdLine.h>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <utility>
#include <string>
#include "ProblemDatas.hpp"
#include "GASolver.hpp"
#include "HESolver.hpp"

#define _CMD_HEADER                                             \
  "Foundation Operational Research Challenge\n"                 \
  "   Developed by Biagio Festa and Alessandro Erba <2017>"

#define _PROGRAM_VERSION "0.0.1"

class FORCH_Program {
 public:
  FORCH_Program();
  void run(int argc, char** argv);

 private:
  using PTR_CmdArg = std::unique_ptr<TCLAP::Arg>;
  using StringArg = TCLAP::ValueArg<std::string>;
  using FloatArg = TCLAP::ValueArg<float>;
  using IntArg = TCLAP::ValueArg<int>;
  using FlagArg = TCLAP::SwitchArg;
  static const char CMD_HEADER[];
  static const char VERSION[];

  TCLAP::CmdLine m_cmd_parser;
  std::map<std::string, PTR_CmdArg> m_cmd_args;
  std::shared_ptr<for_ch::ProblemDatas> mp_problem;

  void init_command_line_parse();
  void print_solution(const std::string& filename,
                      const std::vector<bool>& solution) const;
};

const char FORCH_Program::CMD_HEADER[] = _CMD_HEADER;
const char FORCH_Program::VERSION[] = _PROGRAM_VERSION;

int main(int argc, char *argv[]) {
  FORCH_Program program;
  try {
    program.run(argc, argv);
  } catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    return -1;
  }
  return 0;
}

FORCH_Program::FORCH_Program() :
    m_cmd_parser(CMD_HEADER, ' ', VERSION) {
}

void FORCH_Program::init_command_line_parse() {
  PTR_CmdArg modelFile(
      new StringArg("d",
                    "data-file",
                    "The dat file which describes the problem data",
                    true,
                    "none",
                    "string"));
  m_cmd_args.insert(std::make_pair("data", std::move(modelFile)));

  PTR_CmdArg custom_crossover(
      new FlagArg("z",
                  "disable-custom-crossover",
                  "Flag to DISABLE custom crossover function",
                  false));
  m_cmd_args.insert(std::make_pair("disable-personal-crossover",
                                   std::move(custom_crossover)));

  PTR_CmdArg seconds_max(
      new IntArg("t",
                 "time-genetic",
                 "The time in SECONDS to wait for genetic algorithm",
                 false,
                 10 * 60,
                 "integer"));
  m_cmd_args.insert(std::make_pair("time-genetic",
                                   std::move(seconds_max)));

  for (const auto& arg : m_cmd_args) {
    m_cmd_parser.add(*(arg.second));
  }
}

void FORCH_Program::run(int argc, char** argv) {
  init_command_line_parse();
  m_cmd_parser.parse(argc, argv);

  // Get parameter
  const std::string& dat_filename = dynamic_cast<StringArg*>(
        m_cmd_args.at("data").get())->getValue();
  mp_problem = std::make_shared<for_ch::ProblemDatas>();
  mp_problem->parse_problem_dat(dat_filename);

  // First phase launch HESolver
  std::vector<bool> he_solution;
  for_ch::HESolver hesolver(mp_problem);
  bool he_found = hesolver.run(&he_solution);

  // Personal crossove flag
  for_ch::GASolver gasolver(mp_problem);
  const bool& disable_custom_crossover = dynamic_cast<FlagArg*>(
      m_cmd_args.at("disable-personal-crossover").get())->getValue();
  const int& time = dynamic_cast<IntArg*>(
      m_cmd_args.at("time-genetic").get())->getValue();
  gasolver.set_max_time_seconds(time);
  gasolver.set_flag_custom_crossover(~disable_custom_crossover);

  // Second pgase launch GASolver
  std::vector<bool> ga_solution = gasolver.run(
      argc, argv,
      (he_found == true ? &he_solution : nullptr));

  // Print solution
  std::string solution_filename =
      dat_filename.substr(0, dat_filename.find_last_of('.')) + ".sol";
  print_solution(solution_filename, ga_solution);
  std::cout << "Solution written in the file: '" << solution_filename << "'\n";
}

void FORCH_Program::print_solution(const std::string& filename,
                                   const std::vector<bool>& solution) const {
  const unsigned num_edges = mp_problem->m_numEdges;
  assert(solution.size() == num_edges);

  std::ofstream file;
  file.open(filename);
  if (file.fail()) {
    throw std::runtime_error(std::string(
        "The output file '" + filename + "' cannot be open"));
  }

  for (EdgeIndex i = 0; i < num_edges; ++i) {
    assert(i < solution.size());
    bool edge_active = solution[i];
    if (edge_active) {
      const EdgeLink& link = mp_problem->m_mapEdge_index2link.at(i);
      file << link.first << ' ' << link.second << '\n';
    }
  }

  file.close();
}
