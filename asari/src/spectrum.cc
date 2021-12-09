#include "specToChrom.h"

spectrum::spectrum(int &n_specs, double &RT, int &id, double * rts, double * intns) {
  n_spec_ = n_specs;
  rt_ = RT;
  id_ = id;
  rts_ = rts;
  intns_ = intns;
}
