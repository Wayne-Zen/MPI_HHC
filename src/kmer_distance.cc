#include "kmer_distance.h"
#include <string>
#include <unordered_map>

namespace hhc {

KmerDistance::KmerDistance(int K) : K_(K) {}

float
KmerDistance::operator() (const std::string s0, const std::string s1) {
  std::unordered_map<std::string, int> kmer_cnt;
  std::string kmer;
  int l0 = s0.size();
  int l1 = s1.size();
  for (int i = 0; i < l0 - K_ + 1; ++i) {
    kmer = s0.substr(i, K_);
    ++kmer_cnt[kmer];
  } 
  int cnt = 0;
  for (int i = 0; i < l1 - K_ + 1; ++i) {
    kmer = s1.substr(i, K_);
    if (kmer_cnt[kmer] > 0) {
      ++cnt;
      --kmer_cnt[kmer];
    }
  }
  return static_cast<float>(cnt)/(std::min(l0, l1) - K_ + 1);
}

}