#include <mpi.h>

#include "type.h"
#include "fasta_reader.h"
#include "lahdc.h"
#include "parameter.h"
#include "logger.h"

int main(int argc, char** argv) {
  // initi
  MPI_Init(NULL, NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  hhc::Logger logger(world_size, world_rank);
  logger.open();
  logger.greet();
  struct hhc::Parameter param  = {};
  hhc::parameter::parse_cmd_arguments(argc, argv, param);
  
  // load file
  hhc::ReadIndex total_seq_num = 0;
  hhc::ReadLength max_len = 0;
  std::vector<std::string> headers;
  std::vector<std::string> reads;
  const char* input_fp = param.input_fp.c_str();
  hhc::FastaReader fa_reader(input_fp, world_size, world_rank);
  fa_reader.load_data(total_seq_num, max_len, headers, reads, logger);


  // run LAHDC
  std::vector<hhc::ClusterIndex> assignments(reads.size(), 0);
  hhc::LAHDC lahdc(reads, assignments, 
                   param.kmer, param.cluster_size_threshold,
                   world_size, world_rank);
  lahdc.run(logger);


  // write file


  // clean up
  logger.close();
  MPI_Finalize();
  return 0;
}