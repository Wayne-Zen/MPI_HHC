#include "fasta_reader.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>

#include "type.h"
#include "logger.h"



namespace hhc {

FastaReader::FastaReader(char* input_fp, int world_size, int world_rank)
    : input_fp_(input_fp),
      world_size_(world_size),
      world_rank_(world_rank) {}

void FastaReader::load_data(ReadIndex& total_seq_num, ReadLength& max_len,
                             std::vector<std::string>& headers, 
                             std::vector<std::string>& reads,
                             Logger& logger) {
  get_global_info(total_seq_num, max_len);
  logger.log << "Total sequence number(global): "
                      << total_seq_num << std::endl;
  logger.log << "Max read length(global): "
                      << max_len << std::endl;
  ReadIndex quota = total_seq_num / world_size_;
  ReadIndex extra = total_seq_num % world_size_;
  std::vector<ReadIndex> starting_index(1, 0);
  for (int i = 0; i != world_size_; ++i) {
    ReadIndex num = quota + (i < extra ? 1 : 0);
    starting_index.push_back(
        num + starting_index[starting_index.size() - 1]);
  }
  load_chunk(headers, reads, starting_index);
  logger.log << "local header,reads number: "
                      << headers.size() << "," << reads.size() << std::endl;
}

void FastaReader::get_global_info(ReadIndex& total_seq_num, 
                                  ReadLength& max_len) {
  std::ifstream in;
  in.open(input_fp_);
  if (!in) {
    std::cerr << "error opening file" << std::string(input_fp_) << std::endl;
  }
  std::string line;
  ReadIndex line_cnt = 0;
  total_seq_num = 0;
  max_len = 0;
  while (std::getline(in, line)) {
    if (line_cnt % 2 == 1) {
      ++total_seq_num;
      max_len = line.size() > max_len ? line.size() : max_len;
    }
    ++line_cnt;
  }
  in.close();
  in.clear();
}

void FastaReader::load_chunk(std::vector<std::string>& headers, 
                             std::vector<std::string>& reads,
                             std::vector<ReadIndex>& starting_index) {
  std::ifstream in;
  in.open(input_fp_);
  if (!in) {
    std::cerr << "error opening file" << std::string(input_fp_) << std::endl;
  }
  std::string line;
  ReadIndex line_cnt = 0;
  ReadIndex seq_cnt = 0;
  while (std::getline(in, line)) {
    if (starting_index[world_rank_] <= seq_cnt &&
        seq_cnt < starting_index[world_rank_ + 1]) {
      if (line_cnt % 2 == 0) {
        headers.push_back(line);
      } else {
        reads.push_back(line);
      }
    }
    ++line_cnt;
    seq_cnt = line_cnt / 2;
  }
  in.close();
  in.clear();
}

}