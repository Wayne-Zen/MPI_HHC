#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>

namespace hhc {

struct Parameter {
  int kmer;
  int cluster_size_threshold;
  std::string input_fp;
  std::string output_fp;
};

namespace parameter {
void parse_cmd_arguments(int argc, char** argv, Parameter& param);
}

}

#endif