#ifndef KMER_DISTANCE_H
#define KMER_DISTANCE_H

#include <string>

namespace hhc {

class KmerDistance {
 public:
  explicit KmerDistance(int K);
  float operator() (const std::string s0, const std::string s1);
 private:
  int K_;
};

}

#endif