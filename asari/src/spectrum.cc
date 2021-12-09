#include "specToChrom.h"
#include <cstring>

namespace asaristc {

spectrum::spectrum(int &n_pts, double &RT, int &id, int &offset) {
  n_pts_ = n_pts;
  rt_ = RT;
  id_ = id;
  offset_ = offset;
}

void spectrum::pointers_to_offset(double* study_start) {
  mzs_ = study_start + offset_;
  intns_ = study_start + offset_;
}

void spectrum::copy_values(int &count, double* mzs, double * intns) {
  std::memcpy( (char*) mzs_, (char*) mzs, count);
  std::memcpy( (char*) intns_, (char*) intns, count);
}

}// namespace asaristc
