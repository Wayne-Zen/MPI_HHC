#include "logger.h"

#include <cassert>

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "type.h"

namespace hhc {

Logger::Logger(int world_size, int world_rank)
    : world_size_(world_size),
      world_rank_(world_rank) {
  std::ostringstream log_fp_ss;
  log_fp_ss << "log/node_" << world_rank
            << "_of_" << world_size
            << ".txt" << std::flush;
  log_fp_ = log_fp_ss.str();
  
}

void Logger::open() {
  log.open(log_fp_.c_str());
  if (!log) {
      std::cerr << "error opening file:" 
                << log_fp_
                << std::endl;
  }
}

void Logger::close() {
  log << "##### END #####" << std::endl;
  log << std::endl;
  log << std::endl;
  log.close();
  log.clear();
}

void Logger::greet() {
  log << "### Logger starts on node "
      << world_rank_ << " of " 
      << world_size_ << " ###" << std::endl; 
}

template<>
void Logger::show_vector_oneline<ReadLocation>(
    std::string header,
    std::vector<ReadLocation>& vect) {
  log << std::left << std::setw(kHeaderWidth) 
      << header  << " :";
  for (std::vector<ReadLocation>::size_type i = 0;
       i != vect.size(); ++i) {
    std::ostringstream ss;
    ss << "("
       << "m" << vect[i].machine_index << ","
       << "r" << vect[i].read_index_in_cluster << ")" << std::flush;
    log << std::setw(kCellWidth) << ss.str();
        
  }  
  log << std::endl;
}

template<>
void Logger::show_vector_oneline<ReadInfo>(
    std::string header,
    std::vector<ReadInfo>& vect) {
  log << std::left << std::setw(kHeaderWidth) 
      << header  << " :";
  for (std::vector<ReadInfo>::size_type i = 0;
       i != vect.size(); ++i) {
    std::ostringstream ss;
    ss << "(" 
       << "m" << vect[i].machine_index << ","
       << "c" << vect[i].cluster_label << "," 
       << "s" << vect[i].read_len
       << ")" << std::flush;
    log << std::setw(kCellWidth) << ss.str();
        
  }  
  log << std::endl;
}

template<>
void Logger::show_vector_oneline< std::vector<std::string> >(
    std::string header,
    std::vector< std::vector<std::string> >& vect) {
  log << std::left << std::setw(kHeaderWidth) 
      << header  << " :";
  for (std::vector< std::vector<std::string> >::size_type i = 0;
       i != vect.size(); ++i) {
    std::vector<std::string>& landmark_collection = vect[i];
    assert(landmark_collection.size() >= 1);
    std::string& landmark = landmark_collection[landmark_collection.size() - 1];
    std::ostringstream ss;
    ss << landmark.substr(100, 15) << "..." << std::flush;
    log << std::setw(kCellWidth) << ss.str();
  }  
  log << std::endl;
}

}