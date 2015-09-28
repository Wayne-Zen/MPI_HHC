#ifndef TYPE_H
#define TYPE_H

namespace hhc {
  typedef int ClusterIndex;
  typedef long ReadIndex;
  typedef int ReadLength;

  typedef struct {
    int machine_index;
    ReadIndex read_index_in_cluster;
  } ReadLocation;

  typedef struct {
    int machine_index;
    ClusterIndex cluster_label;
    ReadLength read_len;
  } ReadInfo;
}

#endif