#ifndef BARRIERS_PARAMETERS_H_INCLUDED
#define BARRIERS_PARAMETERS_H_INCLUDED

#include <string>
#include "cmd_line.h"
#include <fstream>

struct param_t
{
  int pop_size = 10000;
  int Gmax = 10000;
  double mutation_rate = 0.001;
  double mutation_shape = 0.01;
  double tolerance = 0.1;
  int genome_length = 4;
  bool chain = false;
  std::string outdir;

};

cmd::cmd_line_parser config_file_parser(const std::string& config)
{
  std::ifstream is(config);
  if (!is) throw cmd::parse_error("can't open config file");
  std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
  return cmd::cmd_line_parser(str);
}



#endif