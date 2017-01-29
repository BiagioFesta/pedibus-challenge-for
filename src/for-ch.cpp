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
  static const char CMD_HEADER[];
  static const char VERSION[];

  TCLAP::CmdLine m_cmd_parser;
  std::map<std::string, PTR_CmdArg> m_cmd_args;
  std::shared_ptr<for_ch::ProblemDatas> mp_problem;

  void init_command_line_parse();
};

const char FORCH_Program::CMD_HEADER[] = _CMD_HEADER;
const char FORCH_Program::VERSION[] = _PROGRAM_VERSION;

int main(int argc, char *argv[]) {
  FORCH_Program program;
  try {
    program.run(argc, argv);
  } catch (const std::exception& err) {
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

  for (const auto& arg : m_cmd_args) {
    m_cmd_parser.add(*(arg.second));
  }
}

void FORCH_Program::run(int argc, char** argv) {
  init_command_line_parse();
  m_cmd_parser.parse(argc, argv);

  const std::string& dat_filename = dynamic_cast<StringArg*>(
        m_cmd_args.at("data").get())->getValue();
  mp_problem = std::make_shared<for_ch::ProblemDatas>();
  mp_problem->parse_problem_dat(dat_filename);

  // First phase launch HESolver
  std::vector<bool> he_solution;
  for_ch::HESolver hesolver(mp_problem);
  bool he_found = hesolver.run(&he_solution);

  for_ch::GASolver gasolver(mp_problem);
  gasolver.run(argc, argv,
               (he_found == true ? &he_solution : nullptr));
}
