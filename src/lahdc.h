#ifndef LAHDC_H
#define LAHDC_H

#include <string>
#include <vector>
#include <map>
#include <cassert>

#include "kmer_distance.h"
#include "logger.h"
#include "type.h"

namespace hhc {

class LAHDC {
 public:  
  LAHDC(std::vector<std::string>& reads,
        std::vector<ClusterIndex>& assignments, 
        int kmer,
        int cluster_size_threshold,
        int world_size, 
        int world_rank);
  void run(Logger &logger);
 
 private:
  void get_splitting_cluster_order(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::map<ClusterIndex, int>& splitting_cluster_order);

  void get_splitting_cluster(
      std::vector<ClusterIndex>& ready_cluster_vect,
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
      std::vector<ReadIndex>& read_cnt_hist,
      std::vector<ReadIndex>& read_cnt_hist_global);

  void get_seed_candidate(
      std::vector<ClusterIndex>& splitting_cluster_vect, 
      std::vector<ReadIndex>& seed_candidate_count_vect,
      std::vector< std::vector<ReadIndex> >& seed_candidate_index_list_vect);

  void select_from_candidates(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector<ReadIndex>& candidate_count_vect,
      std::vector< std::vector<ReadIndex> >& candidate_index_list_vect,
      std::vector< std::vector<std::string> >& cluster_landmarks_vect,
      Logger& logger);

  void select_read_location(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector<ReadIndex>& candidate_count_vect_global,
      std::vector<ReadLocation>& selected_read_location_vect);

  void get_local_read_info(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector<ReadLocation>& selected_read_location_vect,
      std::vector< std::vector<ReadIndex> >& candidate_index_list_vect,
      std::vector<ReadInfo>& read_info_vect_local,
      std::vector<int>& read_info_recvcount,
      std::vector<int>& read_info_displs,
      std::string& read_collection_content);

  void get_ready_for_transmission(
      std::vector<ReadInfo>& read_info_vect_global,
      std::string& selected_read_local,
      std::vector<int>& read_recvcount,
      std::vector<int>& read_displs,
      std::string& selected_read_global);

  void append_to_landmarks(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector<ReadInfo>& read_info_vect_global,
      std::string& selected_read_global,
      std::vector< std::vector<std::string> >& cluster_landmarks_vect);

  ReadIndex get_landmarks_num_target(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector<ClusterIndex>& ready_cluster_vect,
      std::vector<ReadIndex>& read_cnt_hist_global);

  void calculate_dist(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector< std::vector<std::string> >& cluster_landmarks_vect,
      std::vector< std::vector<float> >& dist_matrix, 
      std::vector<float>& sum_dist);

  void get_median_of_medians(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector<float>& sum_dist,
      std::vector<float>& mom);

  void get_landmark_candidate(
      std::vector<ClusterIndex>& splitting_cluster_vect, 
      std::vector<float>& mom,
      std::vector<float>& sum_dist,
      std::vector<ReadIndex>& landmark_candidate_count_vect,
      std::vector< std::vector<ReadIndex> >& landmark_candidate_index_list_vect);

  void split(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
      std::vector< std::vector<std::string> >& cluster_landmarks_vect,
      ReadIndex landmarks_num_target,
      std::vector< std::vector<float> >& dist_matrix,
      Logger& logger);

  void split_landmarks(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
      std::vector< std::vector<std::string> >& cluster_landmarks_vect,
      ReadIndex landmarks_num_target,
      std::vector<ClusterIndex>& landmark_assignment_local,
      std::vector<int>& landmark_assignment_recvcount,
      std::vector<int>& landmark_assignment_displs);

  void spectral_clustering(
      std::vector<std::string>& landmarks,
      std::vector<ClusterIndex>& label_choice,
      std::vector<ClusterIndex>& landmark_assignments);

  void split_reads(
      std::vector<ClusterIndex>& splitting_cluster_vect,
      std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
      ReadIndex landmarks_num_target,
      std::vector<ClusterIndex>& landmark_assignment_global,
      std::vector< std::vector<float> >& dist_matrix);

  float get_ave_dist(
      std::vector< std::vector<float> >& dist_matrix,
      ReadIndex read_index,
      std::vector<ReadIndex>& landmark_index_list);

  void update_ready_clusters(
      std::vector<ClusterIndex>& ready_cluster_vect,
      std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect);

  template <typename T_key, typename T_value>
  void build_map_from_vectors(std::vector<T_key>& key_vect, 
                              std::vector<T_value>& value_vect,
                              std::map<T_key, T_value>& res_map) {
    assert(key_vect.size() == value_vect.size());
    for (typename std::vector<T_key>::size_type i = 0;
         i != key_vect.size(); ++i) {
      res_map[key_vect[i]] = value_vect[i];
    }
  }

  std::vector<std::string> reads_;
  std::vector<ClusterIndex> assignments_;  
  int world_size_;
  int world_rank_;
  int max_cluster_label_;
  int cluster_size_threshold_;
  KmerDistance kmer_distance;
  };

}

#endif