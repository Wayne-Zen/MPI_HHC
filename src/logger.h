#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>

#include "type.h"

namespace hhc {

class Logger {
 public:
  Logger(int world_size, int world_rank);
  void open();
  void close();
  void greet();

  template <typename T>
  void show_vector_oneline(std::string header,
                           std::vector<T>& vect) {
    log << std::left << std::setw(kHeaderWidth) 
        << header  << " :";
    for (typename std::vector<T>::size_type i = 0;
         i != vect.size(); ++i) {
      log << std::setw(kCellWidth) << vect[i];
    }  
    log << std::endl;
  }

  void show_vector_oneline(
      std::string header,
      std::vector<ReadLocation>& vect);

  void show_vector_oneline(
    std::string header,
    std::vector<ReadInfo>& vect);

  void show_vector_oneline(
    std::string header,
    std::vector< std::vector<std::string> >& vect);

  template <typename T_key, typename T>
  void show_vmap_multiline(std::string header,
                           std::map< T_key, std::vector<T> >& vmap) {
    std::string item_sep = "\n";
    std::string value_sep = ", ";
    log << header  << " {" << std::endl;
    //typename std::map< T_key, std::vector<T> >::iterator iter_last = --vmap.end();
    for (typename std::map< T_key, std::vector<T> >::iterator iter = vmap.begin();
         iter != vmap.end(); ++iter) {
      log << "  " << iter->first << ":";
      typename std::vector<T>& value = iter->second;
      typename std::vector<T>::iterator it_last = --value.end();
      log << "[";
      for (typename std::vector<T>::iterator it = value.begin();
           it != value.end(); ++it) {
        log << *(it);
        if (it != it_last) {
          log << value_sep;
        }
      }
      log << "]";
      log << item_sep;
    }
    log << "}" << std::endl;
  }

  std::ofstream log;

 private:
  std::string log_fp_;
  int world_size_;
  int world_rank_;
  const static int kCellWidth = 20; 
  const static int kHeaderWidth = 30;
};

}

#endif