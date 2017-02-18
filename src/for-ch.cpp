// Copyright 2017 <Biagio Festa>
#include <tclap/CmdLine.h>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <limits>
#include <set>
#include "ProblemDatas.hpp"
#include "GASolver.hpp"
#include "HESolver.hpp"
#include "ASolver.hpp"

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
  using StringArgUnlabel = TCLAP::UnlabeledValueArg<std::string>;
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
                      const Solution& solution) const;
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
      new StringArgUnlabel("data-file",
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

  PTR_CmdArg verbose(
      new FlagArg("v",
                  "verbose-output",
                  "Flag to display verbose output",
                  false));
  m_cmd_args.insert(std::make_pair("verbose",
                                   std::move(verbose)));

  PTR_CmdArg heuristic_a(
      new FloatArg("a",
                   "heuristic-a",
                   "Set the 'a' parameter for heuristic",
                   false,
                   0.f,
                   "float"));
  m_cmd_args.insert(std::make_pair("heuristic-a",
                                   std::move(heuristic_a)));
  PTR_CmdArg heuristic_b(
      new FloatArg("b",
                   "heuristic-b",
                   "Set the 'b' parameter for heuristic",
                   false,
                   0.f,
                   "float"));
  m_cmd_args.insert(std::make_pair("heuristic-b",
                                   std::move(heuristic_b)));
  PTR_CmdArg heuristic_c(
      new FloatArg("c",
                   "heuristic-c",
                   "Set the 'c' parameter for heuristic",
                   false,
                   0.f,
                   "float"));
  m_cmd_args.insert(std::make_pair("heuristic-c",
                                   std::move(heuristic_c)));
  PTR_CmdArg heuristic_d(
      new FloatArg("d",
                   "heuristic-d",
                   "Set the 'd' parameter for heuristic",
                   false,
                   0.f,
                   "float"));
  m_cmd_args.insert(std::make_pair("heuristic-d",
                                   std::move(heuristic_d)));
  PTR_CmdArg heuristic_e(
      new FloatArg("e",
                   "heuristic-e",
                   "Set the 'e' parameter for heuristic",
                   false,
                   0.f,
                   "float"));
  m_cmd_args.insert(std::make_pair("heuristic-e",
                                   std::move(heuristic_e)));

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
  const bool verbose_output = dynamic_cast<FlagArg*>(
      m_cmd_args.at("verbose").get())->getValue();

  // Construct problem
  mp_problem = std::make_shared<for_ch::ProblemDatas>();
  mp_problem->parse_problem_dat(dat_filename);

  // Set of solution
  std::set<Solution> solutions;

  // #############################################
  // Launch Heuristic Algorithm
  Solution he_solution;
  for_ch::HESolver hesolver(*mp_problem);

  // Parameters for heuristic solver
  std::vector<char> params_he = {'a', 'b', 'c', 'd', 'e'};
  for (const auto& cp : params_he) {
    const std::string key_flag = std::string("heuristic-") + cp;
    FloatArg* parg =
        dynamic_cast<FloatArg*>(m_cmd_args.at(key_flag).get());
    if (parg->isSet()) {
      float value = parg->getValue();
      switch (cp) {
        case 'a':
          hesolver.set_param_a(value);
          break;
        case 'b':
          hesolver.set_param_b(value);
          break;
        case 'c':
          hesolver.set_param_c(value);
          break;
        case 'd':
          hesolver.set_param_d(value);
          break;
        case 'e':
          hesolver.set_param_e(value);
          break;
        default:
          // No default should be reached
          assert(false);
      }
    }
  }

  // Launch heuristic algorithm
  bool he_found = hesolver.run(&he_solution);
  assert(he_found == true);
  if (he_found == true) {
    solutions.insert(he_solution);
  }
  // #############################################

  // #############################################
  // Launch the Algorithmic solver
#ifdef NDEBUG
  constexpr unsigned NUM_SOLU_REQ = 10000;
#else
  constexpr unsigned NUM_SOLU_REQ = 100;
#endif
  constexpr unsigned NUM_TRIALS = NUM_SOLU_REQ * 10;

  std::cout << "ASolver running...\n";
  for_ch::ASolver asolver(mp_problem.get());
  Solution temp_solution;
  for (unsigned i = 0; i < NUM_TRIALS && solutions.size() < NUM_SOLU_REQ; ++i) {
    asolver.run(&temp_solution);
    assert(temp_solution.compute_feasibility(*mp_problem) == true);
    solutions.insert(std::move(temp_solution));
  }
  std::cout << "ASolver completed\n"
      "Num Feas. Solution found: " << solutions.size() << "\n";
  // #############################################

  // #############################################
  // Launch Genetic Algorithm
  for_ch::GASolver gasolver(*(mp_problem.get()));
  const bool& disable_custom_crossover = dynamic_cast<FlagArg*>(
      m_cmd_args.at("disable-personal-crossover").get())->getValue();
  const int& time = dynamic_cast<IntArg*>(
      m_cmd_args.at("time-genetic").get())->getValue();
  gasolver.set_max_time_seconds(time);
  gasolver.set_flag_custom_crossover(~disable_custom_crossover);
  gasolver.set_verbose(verbose_output);

  Solution ga_solution;
  if (time > 0) {
    std::cout << "GASolver running...\n";
    ga_solution = gasolver.run(argc, argv, solutions);
    std::cout << "\nGASolver completed\n";
  }
  // #############################################

  // Construct final solution
  Solution final_solution;
  if (time > 0) {
    final_solution = ga_solution;
  } else {
    final_solution = *std::min_element(
        solutions.cbegin(), solutions.cend(),
        [] (const Solution& s1, const Solution& s2) {
          return s1.m_num_leaves < s2.m_num_leaves;
        });
  }

  // Print solution
  const std::string path_complete_binary = std::string(argv[0]);
  const std::string path_binary = path_complete_binary.substr(
      0, path_complete_binary.find_last_of('/'));
  const std::string solution_filename =
      path_binary + "/"
      + dat_filename.substr(dat_filename.find_last_of('/') + 1,
                            dat_filename.find_last_of('.') -
                            dat_filename.find_last_of('/') - 1)
      + ".sol";
  print_solution(solution_filename, final_solution);
  std::cout << "\nSolution written in the file: '" <<
      solution_filename << "'\n";
}

void FORCH_Program::print_solution(const std::string& filename,
                                   const Solution& solution) const {
  const unsigned num_edges = mp_problem->m_numEdges;
  assert(solution.m_active_edges.size() == num_edges);

  std::ofstream file;
  file.open(filename);
  if (file.fail()) {
    throw std::runtime_error(std::string(
        "The output file '" + filename + "' cannot be open"));
  }

  for (EdgeIndex i = 0; i < num_edges; ++i) {
    assert(i < solution.m_active_edges.size());
    bool edge_active = solution.m_active_edges[i];
    if (edge_active) {
      const EdgeLink& link = mp_problem->m_mapEdge_index2link.at(i);
      file << link.first << ' ' << link.second << '\n';
    }
  }

  file.close();
}
