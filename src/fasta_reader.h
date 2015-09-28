#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <vector>
#include <string>

#include "logger.h"
#include "type.h"

namespace hhc {

class FastaReader {
 public:
  FastaReader(char* input_fp, int world_size, int world_rank);
  void load_data(ReadIndex& total_seq_num, ReadLength& max_len, 
                 std::vector<std::string>& headers, 
                 std::vector<std::string>& reads,
                 Logger& logger);
 
 private:
  void get_global_info(ReadIndex& total_seq_num, ReadLength& max_len);
  void load_chunk(std::vector<std::string>& headers, 
                  std::vector<std::string>& reads,
                  std::vector<ReadIndex>& starting_index);
  char* input_fp_;
  int world_size_;
  int world_rank_;
};

}

#endif