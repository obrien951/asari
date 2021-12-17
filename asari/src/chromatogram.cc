#include "specToChrom.h"

namespace asaristc {

chromatogram::chromatogram(int &count, point &peak, double * mzs, double * rts, double * intns) {
  mz_ = peak.mz_;
  data.resize(count);
  for (int i =0; i < count; i++) {
    data[i].mz_ = mzs[i];
    data[i].rt_ = rts[i];
    data[i].intn_ = intns[i];
  }
}

}
