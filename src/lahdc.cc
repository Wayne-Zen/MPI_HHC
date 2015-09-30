#include "lahdc.h"

#include <cstdlib> // rand
#include <cassert>
#include <cmath>

#include <mpi.h>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "armadillo/armadillo"
#include "kmer_distance.h"
#include "type.h"
#include "logger.h"


namespace hhc {

LAHDC::LAHDC(std::vector<std::string>& reads,
             std::vector<ClusterIndex>& assignments, 
             int kmer,
             int cluster_size_threshold,
             int world_size, 
             int world_rank)
    : reads_(reads), assignments_(assignments), 
      world_size_(world_size), world_rank_(world_rank),
      cluster_size_threshold_(cluster_size_threshold),
      kmer_distance(KmerDistance(kmer)) {
  max_cluster_label_ = 0;
  srand(0);
}

void LAHDC::get_splitting_cluster_order(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::map<ClusterIndex, int>& splitting_cluster_order) {
  for (std::vector<ClusterIndex>::size_type i = 0; 
       i != splitting_cluster_vect.size(); ++i) {
    splitting_cluster_order[splitting_cluster_vect[i]] = i;
  }
}


void LAHDC::get_splitting_cluster(
    std::vector<ClusterIndex>& ready_cluster_vect,
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
    std::vector<ReadIndex>& read_cnt_hist,
    std::vector<ReadIndex>& read_cnt_hist_global) {
  std::map<ClusterIndex, ReadIndex> read_cnt_map;
  for (std::vector<std::string>::size_type i = 0; i != reads_.size(); ++i) {
    ++read_cnt_map[assignments_[i]];
  }
  for (std::vector<ClusterIndex>::size_type i = 0; 
       i != ready_cluster_vect.size(); ++i) {
    read_cnt_hist.push_back(read_cnt_map[ready_cluster_vect[i]]);
  }
  MPI_Allreduce(&read_cnt_hist[0],
                &read_cnt_hist_global[0],
                ready_cluster_vect.size(),
                MPI_LONG,
                MPI_SUM,
                MPI_COMM_WORLD);
  for (std::vector<ClusterIndex>::size_type i = 0;
       i != ready_cluster_vect.size(); ++i) {
    if (read_cnt_hist_global[i] > cluster_size_threshold_) {
      splitting_cluster_vect.push_back(ready_cluster_vect[i]);
      std::vector<ClusterIndex> splitted_clusters;
      splitted_clusters.push_back(++max_cluster_label_);
      splitted_clusters.push_back(++max_cluster_label_);
      splitted_clusters_vect.push_back(splitted_clusters);
    }
  }
}

void LAHDC::get_seed_candidate(
    std::vector<ClusterIndex>& splitting_cluster_vect, 
    std::vector<ReadIndex>& seed_candidate_count_vect,
    std::vector< std::vector<ReadIndex> >& seed_candidate_index_list_vect) {
  std::map<ClusterIndex, int> splitting_cluster_order;
  get_splitting_cluster_order(splitting_cluster_vect, splitting_cluster_order);

  for (std::vector<std::string>::size_type i = 0; i != reads_.size(); ++i) {
    ClusterIndex cluster_label = assignments_[i];
    std::vector<ClusterIndex>::size_type order = 
        splitting_cluster_order[cluster_label];
    ++seed_candidate_count_vect[order];
    seed_candidate_index_list_vect[order].push_back(i);
  }
}

void LAHDC::select_from_candidates(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector<ReadIndex>& candidate_count_vect,
    std::vector< std::vector<ReadIndex> >& candidate_index_list_vect,
    std::vector< std::vector<std::string> >& cluster_landmarks_vect,
    Logger& logger) {

  // gather to root to make consistent decision on selection
  ClusterIndex splitting_cluster_num = splitting_cluster_vect.size();
  std::vector<ReadIndex> 
      candidate_count_vect_global(world_size_ * splitting_cluster_num);
  MPI_Gather(&candidate_count_vect[0],
              splitting_cluster_num,
              MPI_LONG,
              &candidate_count_vect_global[0],
              splitting_cluster_num,
              MPI_LONG,
              0,
              MPI_COMM_WORLD);
  
  std::vector<ReadLocation> selected_read_location_vect(splitting_cluster_num);
  if (world_rank_ == 0) {
    select_read_location(splitting_cluster_vect,
                         candidate_count_vect_global,
                         selected_read_location_vect);
  }

  // broadcast decision
  int location_blocklengths[2] = {1, 1};
  MPI_Datatype location_types[2] = {MPI_INT, MPI_LONG};
  MPI_Aint location_displs[2];
  MPI_Aint location_int_extent;
  MPI_Type_extent(MPI_INT, &location_int_extent);
  MPI_Datatype MPI_READ_LOCATION;
  location_displs[0] = (MPI_Aint)0;
  location_displs[1] = location_int_extent;
  MPI_Type_struct(2, location_blocklengths, location_displs, location_types, 
                  &MPI_READ_LOCATION);
  MPI_Type_commit(&MPI_READ_LOCATION);
  MPI_Bcast(&selected_read_location_vect[0],
            splitting_cluster_num,
            MPI_READ_LOCATION,
            0,
            MPI_COMM_WORLD);
  MPI_Type_free(&MPI_READ_LOCATION);
  logger.show_vector_oneline("Selected read location", selected_read_location_vect);
  
  // exchange read info; 
  std::vector<ReadInfo> read_info_vect_local;
  std::vector<int> read_info_recvcount(world_size_, 0);
  std::vector<int> read_info_displs(world_size_, 0);
  std::string read_collection_content;
  std::vector<ReadInfo> read_info_vect_global(splitting_cluster_num);
  get_local_read_info(splitting_cluster_vect,
                      selected_read_location_vect,
                      candidate_index_list_vect,
                      read_info_vect_local,
                      read_info_recvcount,
                      read_info_displs,
                      read_collection_content);
  MPI_Datatype MPI_READ_INFO;
  MPI_Type_contiguous(3, MPI_INT, &MPI_READ_INFO);
  MPI_Type_commit(&MPI_READ_INFO);
  MPI_Allgatherv(&read_info_vect_local[0],
                read_info_recvcount[world_rank_],
                MPI_READ_INFO,
                &read_info_vect_global[0],
                &read_info_recvcount[0],
                &read_info_displs[0],
                MPI_READ_INFO,
                MPI_COMM_WORLD);
  MPI_Type_free(&MPI_READ_INFO);
  // read info order maybe different from splitting_cluster_vect
  logger.show_vector_oneline("Read info", read_info_vect_global);

  // transmit read
  std::string& selected_read_local = read_collection_content;
  std::vector<int> read_recvcount(world_size_, 0);
  std::vector<int> read_displs(world_size_, 0);
  std::string selected_read_global;
  get_ready_for_transmission(read_info_vect_global,
                             selected_read_local,
                             read_recvcount,
                             read_displs,
                             selected_read_global);

  MPI_Allgatherv(&selected_read_local[0],
                 read_recvcount[world_rank_],
                 MPI_CHAR,
                 &selected_read_global[0],
                 &read_recvcount[0],
                 &read_displs[0],
                 MPI_CHAR,
                 MPI_COMM_WORLD);

  // append read to landmark collection
  append_to_landmarks(splitting_cluster_vect,
                      read_info_vect_global,
                      selected_read_global,
                      cluster_landmarks_vect);
  logger.show_vector_oneline("Read content", cluster_landmarks_vect);
}

void LAHDC::append_to_landmarks(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector<ReadInfo>& read_info_vect_global,
    std::string& selected_read_global,
    std::vector< std::vector<std::string> >& cluster_landmarks_vect) {

  std::map<ClusterIndex, int> splitting_cluster_order;
  get_splitting_cluster_order(splitting_cluster_vect, splitting_cluster_order);
  std::vector<std::string> landmarks;
  std::string tmp;
  for (std::string::size_type i = 0; i != selected_read_global.size(); ++i) {
    if (selected_read_global[i] == '^') {
      tmp.clear();
    } else if (selected_read_global[i] == '$') {
      landmarks.push_back(tmp);
    } else {
      tmp.push_back(selected_read_global[i]);
    }
  }
  assert(landmarks.size() == read_info_vect_global.size());
  for (std::vector<ReadInfo>::size_type i = 0; 
       i != read_info_vect_global.size(); ++i) {
    ClusterIndex cluster_label = read_info_vect_global[i].cluster_label;
    cluster_landmarks_vect[splitting_cluster_order[cluster_label]].push_back(
        landmarks[i]);
  }
}

void LAHDC::get_ready_for_transmission(
    std::vector<ReadInfo>& read_info_vect_global,
    std::string& selected_read_local,
    std::vector<int>& read_recvcount,
    std::vector<int>& read_displs,
    std::string& selected_read_global) {

  int total_char_num = 0;
  for (std::vector<ReadInfo>::size_type i = 0;
       i != read_info_vect_global.size(); ++i) {
    int machine_index = read_info_vect_global[i].machine_index;
    read_recvcount[machine_index] += read_info_vect_global[i].read_len;
    total_char_num += read_info_vect_global[i].read_len;
  }
  selected_read_global.resize(total_char_num);
  for (int i = 1; i != world_size_; ++i) {
    read_displs[i] = read_displs[i - 1] + read_recvcount[i - 1];
  }
}

void LAHDC::get_local_read_info(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector<ReadLocation>& selected_read_location_vect,
    std::vector< std::vector<ReadIndex> >& candidate_index_list_vect,
    std::vector<ReadInfo>& read_info_vect_local,
    std::vector<int>& read_info_recvcount,
    std::vector<int>& read_info_displs,
    std::string& read_collection_content) {

  std::ostringstream ss;
  for (std::vector<ReadLocation>::size_type i = 0; 
       i != selected_read_location_vect.size(); ++i) {
    int machine_index = selected_read_location_vect[i].machine_index;
    ReadIndex read_index_in_cluster = 
                        selected_read_location_vect[i].read_index_in_cluster;
    if (machine_index == world_rank_) {
      ReadIndex selected_read_index = 
          candidate_index_list_vect[i][read_index_in_cluster];
      ss << "^" << reads_[selected_read_index] << "$" << std::flush;
      ReadInfo read_info = {};
      read_info.machine_index = world_rank_;
      read_info.cluster_label = splitting_cluster_vect[i];
      read_info.read_len = reads_[selected_read_index].size() + 2; // '^' '$'
      read_info_vect_local.push_back(read_info);
    }
    ++read_info_recvcount[machine_index];
  }
  for (int i = 1; i != world_size_; ++i) {
    read_info_displs[i] = read_info_displs[i - 1] + read_info_recvcount[i - 1];
  }
  read_collection_content = ss.str();
}

void LAHDC::select_read_location(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector<ReadIndex>& candidate_count_vect_global,
    std::vector<ReadLocation>& selected_read_location_vect) {

  ClusterIndex splitting_cluster_num = splitting_cluster_vect.size();
  for (std::vector<ClusterIndex>::size_type i = 0;
       i != splitting_cluster_vect.size(); ++i) {
    ReadIndex sum_count = 0;
    for (int j = 0; j != world_size_; ++j) {
      sum_count += candidate_count_vect_global[j * splitting_cluster_num + i];
    }
    ReadIndex rand_index = rand() % sum_count;
    ReadLocation loc = {};
    for (int j = 0; j != world_size_; ++j) {
      ReadIndex cur_machine_num = 
          candidate_count_vect_global[j * splitting_cluster_num + i];
      if (rand_index < cur_machine_num) {
        loc.machine_index = j;
        loc.read_index_in_cluster = rand_index;
        std::cout << "Decision made at Root: Selected candidate for cluster["
                  << splitting_cluster_vect[i] << "], location: ["
                  << loc.machine_index << ", "
                  << loc.read_index_in_cluster << "]" << std::endl;
        break;
      } else {
        rand_index -= cur_machine_num;
      }
    }
    selected_read_location_vect[i] = loc;
  }
}

ReadIndex LAHDC::get_landmarks_num_target(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector<ClusterIndex>& ready_cluster_vect,
    std::vector<ReadIndex>& read_cnt_hist_global) {

  std::map<ClusterIndex, ReadIndex> hist_map;
  build_map_from_vectors(ready_cluster_vect, 
                         read_cnt_hist_global,
                         hist_map);
  ReadIndex max_num = 0;
  for (std::vector<ClusterIndex>::size_type i = 0;
       i != splitting_cluster_vect.size(); ++i) {
    if (hist_map[splitting_cluster_vect[i]] > max_num) {
      max_num = hist_map[splitting_cluster_vect[i]];
    }
  }
  return static_cast<int>(ceil(log2(max_num)));
}

void LAHDC::calculate_dist(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector< std::vector<std::string> >& cluster_landmarks_vect,
    std::vector< std::vector<float> >& dist_matrix, 
    std::vector<float>& sum_dist) {
  ClusterIndex splitting_cluster_num = splitting_cluster_vect.size();
  std::vector<std::string> latest_landmark_vect(splitting_cluster_num);
  assert(cluster_landmarks_vect.size() == splitting_cluster_num);
  for (std::vector<ClusterIndex>::size_type i = 0;
       i != splitting_cluster_vect.size(); ++i) {
    std::vector<std::string>& landmarks = cluster_landmarks_vect[i];
    latest_landmark_vect[i] = landmarks[landmarks.size() - 1];
  }

  std::map<ClusterIndex, std::string> latest_landmark_map;
  build_map_from_vectors(splitting_cluster_vect, 
                         latest_landmark_vect,
                         latest_landmark_map);

  std::vector<float> new_dist(reads_.size(), 0);
  for (std::vector<std::string>::size_type i = 0; i != reads_.size(); ++i) {
    ClusterIndex cluster_label = assignments_[i];
    new_dist[i] = kmer_distance(reads_[i], latest_landmark_map[cluster_label]);
    sum_dist[i] += new_dist[i];
  }
  dist_matrix.push_back(new_dist);
}

void LAHDC::get_median_of_medians(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector<float>& sum_dist,
    std::vector<float>& mom) {

  ClusterIndex splitting_cluster_num = splitting_cluster_vect.size();
  std::map<ClusterIndex, int> splitting_cluster_order;
  get_splitting_cluster_order(splitting_cluster_vect, splitting_cluster_order);

  assert(sum_dist.size() == reads_.size());
  std::vector< std::vector<float> > cluster_dists_vect(splitting_cluster_num);
  for (std::vector<float>::size_type i = 0; i != sum_dist.size(); ++i) {
    ClusterIndex cluster_label = assignments_[i];
    cluster_dists_vect[splitting_cluster_order[cluster_label]].push_back(sum_dist[i]);
  }
  std::vector<float> cluster_median_vect_local(splitting_cluster_num);
  std::vector<float> cluster_median_vect_global(splitting_cluster_num * world_size_);
  for (std::vector<ClusterIndex>::size_type i = 0; 
       i != splitting_cluster_vect.size(); ++i) {
    std::vector<float>& dist_vect = cluster_dists_vect[i];
    std::nth_element(dist_vect.begin(), 
                     dist_vect.begin() + dist_vect.size() / 2, 
                     dist_vect.end());
    cluster_median_vect_local[i] = dist_vect[dist_vect.size() / 2 - 1];
  }
  MPI_Allgather(&cluster_median_vect_local[0],
                splitting_cluster_num,
                MPI_FLOAT,
                &cluster_median_vect_global[0],
                splitting_cluster_num,
                MPI_FLOAT,
                MPI_COMM_WORLD);
  for (std::vector<ClusterIndex>::size_type i = 0; 
       i != splitting_cluster_vect.size(); ++i) {
    std::vector<float> tmp;
    for (int j = 0; j != world_size_; ++j) {
      tmp.push_back(cluster_median_vect_global[j * splitting_cluster_num + i]);
    }
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
    mom[i] = tmp[tmp.size() / 2 - 1];
  }
}

void LAHDC::get_landmark_candidate(
    std::vector<ClusterIndex>& splitting_cluster_vect, 
    std::vector<float>& mom,
    std::vector<float>& sum_dist,
    std::vector<ReadIndex>& landmark_candidate_count_vect,
    std::vector< std::vector<ReadIndex> >& landmark_candidate_index_list_vect) {
  std::map<ClusterIndex, int> splitting_cluster_order;
  get_splitting_cluster_order(splitting_cluster_vect, splitting_cluster_order);

  assert(reads_.size() == sum_dist.size());
  for (std::vector<std::string>::size_type i = 0; i != reads_.size(); ++i) {
    ClusterIndex cluster_label = assignments_[i];
    std::vector<ClusterIndex>::size_type order = 
        splitting_cluster_order[cluster_label];
    if (sum_dist[i] > mom[order]) {
      ++landmark_candidate_count_vect[order];
      landmark_candidate_index_list_vect[order].push_back(i);
    }
  }
}

void LAHDC::split(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
    std::vector< std::vector<std::string> >& cluster_landmarks_vect,
    ReadIndex landmarks_num_target,
    std::vector< std::vector<float> >& dist_matrix,
    Logger& logger) {
  ClusterIndex splitting_cluster_num = splitting_cluster_vect.size();
  std::vector<ClusterIndex> landmark_assignment_local;
  std::vector<int> landmark_assignment_recvcount(world_size_, 0);
  std::vector<int> landmark_assignment_displs(world_size_, 0);
  std::vector<ClusterIndex> 
      landmark_assignment_global(landmarks_num_target * splitting_cluster_num);
  split_landmarks(splitting_cluster_vect,
                  splitted_clusters_vect,
                  cluster_landmarks_vect,
                  landmarks_num_target,
                  landmark_assignment_local,
                  landmark_assignment_recvcount,
                  landmark_assignment_displs);
  // logger.show_vector_oneline("landmark_assignment_local", landmark_assignment_local);
  // logger.show_vector_oneline("landmark_assignment_recvcount", landmark_assignment_recvcount);
  // logger.show_vector_oneline("landmark_assignment_displs", landmark_assignment_displs);
  MPI_Allgatherv(&landmark_assignment_local[0],
                 landmark_assignment_recvcount[world_rank_],
                 MPI_INT,
                 &landmark_assignment_global[0],
                 &landmark_assignment_recvcount[0],
                 &landmark_assignment_displs[0],
                 MPI_INT,
                 MPI_COMM_WORLD);
  for (std::vector<ClusterIndex>::size_type i = 0;
       i != splitting_cluster_vect.size(); ++i) {

  }
  split_reads(splitting_cluster_vect,
              splitted_clusters_vect,
              landmarks_num_target,
              landmark_assignment_global,
              dist_matrix);
}

void LAHDC::split_reads(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
    ReadIndex landmarks_num_target,
    std::vector<ClusterIndex>& landmark_assignment_global,
    std::vector< std::vector<float> >& dist_matrix) {

  std::map< ClusterIndex, 
            std::map<ClusterIndex, std::vector<ReadIndex> > > ref_map;
  for (std::vector<ClusterIndex>::size_type i = 0;
       i != landmark_assignment_global.size(); ++i) {
    ClusterIndex splitting_cluster_label = 
        splitting_cluster_vect[i / landmarks_num_target];
    ClusterIndex splitted_cluster_label = landmark_assignment_global[i];
    ReadIndex landmark_index = i % landmarks_num_target;
    ref_map[splitting_cluster_label][splitted_cluster_label].
        push_back(landmark_index);
  }  
  for (std::vector<std::string>::size_type i = 0; i != reads_.size(); ++i) {
    ClusterIndex splitting_cluster_label = assignments_[i];
    
    if (ref_map.count(splitting_cluster_label)) {
      assert(ref_map[splitting_cluster_label].size() == 2);
      std::map<ClusterIndex, std::vector<ReadIndex> >::iterator 
          iter = ref_map[splitting_cluster_label].begin();
      ClusterIndex splitted_cluster_label0 = iter->first;
      std::vector<ReadIndex>& landmark_index_list0 = iter->second;
      ++iter;
      ClusterIndex splitted_cluster_label1 = iter->first;
      std::vector<ReadIndex>& landmark_index_list1 = iter->second;

      float ave0 = get_ave_dist(dist_matrix, i, landmark_index_list0);
      float ave1 = get_ave_dist(dist_matrix, i, landmark_index_list1);
      assignments_[i] = ave0 < ave1 ?
                        splitted_cluster_label0 : splitted_cluster_label1;
    }
    
  }
}

float LAHDC::get_ave_dist(
    std::vector< std::vector<float> >& dist_matrix,
    ReadIndex read_index,
    std::vector<ReadIndex>& landmark_index_list) {
  float sum = 0;
  for (std::vector<ReadIndex>::size_type i = 0;
       i != landmark_index_list.size(); ++i) {
    sum += dist_matrix[i][read_index];
  }
  return sum / landmark_index_list.size();
}

void LAHDC::split_landmarks(
    std::vector<ClusterIndex>& splitting_cluster_vect,
    std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect,
    std::vector< std::vector<std::string> >& cluster_landmarks_vect,
    ReadIndex landmarks_num_target,
    std::vector<ClusterIndex>& landmark_assignment_local,
    std::vector<int>& landmark_assignment_recvcount,
    std::vector<int>& landmark_assignment_displs) {

  ClusterIndex splitting_cluster_num = splitting_cluster_vect.size();
  ClusterIndex quota = splitting_cluster_num / world_size_;
  ClusterIndex extra = splitting_cluster_num % world_size_;
  // assign job
  std::vector<int> starting_index(1, 0);
  for (int i = 0; i != world_size_; ++i) {
    int num = quota + (i < extra ? 1 : 0);
    landmark_assignment_recvcount[i] = num * landmarks_num_target;
    if (i > 0) {
      landmark_assignment_displs[i] = 
          landmark_assignment_displs[i - 1] + 
          landmark_assignment_recvcount[i - 1];
    }
    starting_index.push_back(
        num + starting_index[starting_index.size() - 1]);
  }

  for (std::vector<ClusterIndex>::size_type i = 0;
       i != splitting_cluster_vect.size(); ++i) {
    if (starting_index[world_rank_] <= i &&
        i < starting_index[world_rank_ + 1]) {
      spectral_clustering(cluster_landmarks_vect[i],
                          splitted_clusters_vect[i],
                          landmark_assignment_local);
    }
  }
}

void LAHDC::spectral_clustering(
    std::vector<std::string>& landmarks,
    std::vector<ClusterIndex>& label_choice,
    std::vector<ClusterIndex>& landmark_assignments) {
  assert(label_choice.size() == 2);
  int size = landmarks.size();
  arma::fmat W(size, size, arma::fill::zeros);
  for (int i = 0; i < size; ++i) {
      for (int j = i + 1; j < size; ++j) {
          W(i, j) = W(j, i) = kmer_distance(landmarks[i], landmarks[j]);
      }
  }

  arma::fmat D = arma::diagmat(arma::sum(W, 1));
  arma::fmat L = D.t() * (D - W);
  arma::fvec eigval;
  arma::fmat eigvec;
  arma::eig_sym(eigval, eigvec, L);
  arma::uvec index = arma::sort_index(eigval,"ascend");
  arma::fvec eig = eigvec.col(index[1]);

  for (int i = 0; i < eig.size(); ++i) {
      if (eig[i] <= 0) {
          landmark_assignments.push_back(label_choice[0]);
      } else {
          landmark_assignments.push_back(label_choice[1]);
      } 
  }

}

void LAHDC::update_ready_clusters(
    std::vector<ClusterIndex>& ready_cluster_vect,
    std::vector< std::vector<ClusterIndex> >& splitted_clusters_vect) {
  ready_cluster_vect.clear();
  for (std::vector< std::vector<ClusterIndex> >::iterator
       iter = splitted_clusters_vect.begin();
       iter != splitted_clusters_vect.end(); ++iter) {
    ready_cluster_vect.insert(ready_cluster_vect.begin() + 
                                  ready_cluster_vect.size(),
                              iter->begin(), iter->end());
  }
}

void LAHDC::run(Logger& logger) {
  std::vector<ClusterIndex> ready_cluster_vect(1, 0);
  int level = 0;
  while (true) {
    logger.log << "*************" << std::endl;
    logger.log << "level " << level << std::endl;
    logger.show_vector_oneline("Ready clusters", ready_cluster_vect);

    // check if need to futher go deep
    std::vector<ClusterIndex> splitting_cluster_vect;
    std::vector< std::vector<ClusterIndex> > splitted_clusters_vect;
    std::vector<ReadIndex> read_cnt_hist; // size: ready_cluster_vect
    std::vector<ReadIndex> read_cnt_hist_global(ready_cluster_vect.size());
    get_splitting_cluster(ready_cluster_vect,
                          splitting_cluster_vect,
                          splitted_clusters_vect,
                          read_cnt_hist,
                          read_cnt_hist_global);
    ClusterIndex splitting_cluster_num = splitting_cluster_vect.size();
    if (splitting_cluster_num == 0) {
      break;
    }

    logger.show_vector_oneline("Local reads count", read_cnt_hist);
    logger.show_vector_oneline("Global reads count", read_cnt_hist_global);
    
    std::map< ClusterIndex, std::vector<ClusterIndex> > split_to_map;
    build_map_from_vectors(splitting_cluster_vect, 
                           splitted_clusters_vect, 
                           split_to_map);
    logger.show_vmap_multiline("Split to map", split_to_map);
    logger.show_vector_oneline("Splitting clusters", splitting_cluster_vect);

    // select seed
    logger.log << "*SEED(landmark[0])*" << std::endl;
    std::vector<ReadIndex> 
        seed_candidate_count_vect(splitting_cluster_num, 0); // size:: splitting_cluster_vect
    std::vector< std::vector<ReadIndex> > 
        seed_candidate_index_list_vect(splitting_cluster_num);
    get_seed_candidate(splitting_cluster_vect, 
                       seed_candidate_count_vect, 
                       seed_candidate_index_list_vect);
    logger.show_vector_oneline("Candicate count local", seed_candidate_count_vect);

    std::vector< std::vector<std::string> > 
        cluster_landmarks_vect(splitting_cluster_num);
    select_from_candidates(splitting_cluster_vect,
                           seed_candidate_count_vect,
                           seed_candidate_index_list_vect,
                           cluster_landmarks_vect, logger);

    // select landmarks 
    ReadIndex landmarks_num_target = get_landmarks_num_target(
        splitting_cluster_vect, ready_cluster_vect, read_cnt_hist_global);
    logger.log << "Target Landmarks Number: " 
               << landmarks_num_target << std::endl;
    std::vector< std::vector<float> > dist_matrix;
    std::vector<float> sum_dist(reads_.size(), 0);
    for (ReadIndex landmarks_num = 1; 
         landmarks_num != landmarks_num_target; ++landmarks_num) {
      logger.log << "*landmarks[" << landmarks_num << "]*" << std::endl;
      calculate_dist(splitting_cluster_vect, 
                     cluster_landmarks_vect, 
                     dist_matrix, 
                     sum_dist);
      std::vector<float> mom(splitting_cluster_num);
      get_median_of_medians(splitting_cluster_vect, sum_dist, mom);
      logger.show_vector_oneline("Median of medians", mom);

      std::vector<ReadIndex> 
          landmark_candidate_count_vect(splitting_cluster_num, 0); // size:: splitting_cluster_vect
      std::vector< std::vector<ReadIndex> > 
          landmark_candidate_index_list_vect(splitting_cluster_num);
      get_landmark_candidate(splitting_cluster_vect, mom, sum_dist,
                             landmark_candidate_count_vect,
                             landmark_candidate_index_list_vect);
      logger.show_vector_oneline("Candicate count local", landmark_candidate_count_vect);
      select_from_candidates(splitting_cluster_vect,
                             landmark_candidate_count_vect,
                             landmark_candidate_index_list_vect,
                             cluster_landmarks_vect, logger);
    }

    // split
    split(splitting_cluster_vect,
          splitted_clusters_vect,
          cluster_landmarks_vect,
          landmarks_num_target,
          dist_matrix,
          logger);
    update_ready_clusters(ready_cluster_vect, splitted_clusters_vect);
    //break;
    ++level;
  }
}

}