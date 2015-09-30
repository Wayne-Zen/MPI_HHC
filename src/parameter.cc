#include "parameter.h"
#include <string>
#include "tclap/Cmdline.h"

namespace hhc {

namespace parameter {
void parse_cmd_arguments(int argc, char** argv, Parameter& param) {
  try {
    TCLAP::CmdLine cmd("Command description message", ' ', "0.1.0");
    TCLAP::ValueArg<int> kmer_arg("k", "kmer", "Kmer length", false, 6, "int");
    TCLAP::ValueArg<int> cluster_size_threshold_arg(
        "s", "stop_threshold",
        "Stop splitting when cluster size is smaller than stop_threshold",  
        false, 2000, "int");
    TCLAP::ValueArg<std::string> input_fp_arg(
        "i", "input", "input file path", 
        true, "", "string");
    TCLAP::ValueArg<std::string> output_fp_arg(
        "o", "output", "output file path prefix",
        false, "", "string");
    cmd.add(kmer_arg);
    cmd.add(cluster_size_threshold_arg);
    cmd.add(input_fp_arg);
    cmd.add(output_fp_arg);
    cmd.parse(argc, argv);
    param.kmer = kmer_arg.getValue();
    param.cluster_size_threshold = cluster_size_threshold_arg.getValue();
    param.input_fp = input_fp_arg.getValue();
    param.output_fp = output_fp_arg.getValue();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error()
              << "for arg " << e.argId() << std::endl;
  } // try catch
} // parse_cmd_arguments

} // parameter

} // hhc