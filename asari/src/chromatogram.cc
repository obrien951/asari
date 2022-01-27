#include "specToChrom.h"
#include <iostream>

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

void specToChrom::get_chrom(int &ind, py::list &mz_s, py::list &intnss, 
                            py::list &rts){
  double hold;
  if (chroms_[ind].l_start_ == nullptr && chroms_[ind].h_start_ == nullptr) {
    for (int i = 0; i < chroms_[ind].chrom_len_; i++) {
      mz_s[i] = chroms_[ind].s_start_[i].mz_;
      intnss[i] = chroms_[ind].s_start_[i].intensity_;
      rts[i] = spectra_[chroms_[ind].s_start_[i].spec_id_].get_rt();
    }
  } else {
    std::fill(mz_quai_.begin(), mz_quai_.end() , 0.0);
    std::fill(in_quai_.begin(), in_quai_.end() , 0.0);
    std::fill(rt_quai_.begin(), rt_quai_.end() , 0.0);
    int it = 0;
    int count = 0;
    while ( count < chroms_[ind].l_count_ ) {
      if ( chroms_[ind].l_start_[it].chrom_id_ == ind ) {
        mz_quai_[ chroms_[ind].l_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        += chroms_[ind].l_start_[it].mz_;
        in_quai_[ chroms_[ind].l_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        += chroms_[ind].l_start_[it].intensity_;
        rt_quai_[ chroms_[ind].l_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        = spectra_[chroms_[ind].l_start_[it].spec_id_].get_rt();
        count++;
      }
      it++;
    }
    it = 0;
    count = 0;
    while (count < chroms_[ind].s_count_) {
      if ( chroms_[ind].s_start_[it].chrom_id_ == ind ) {
        mz_quai_[ chroms_[ind].s_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        += chroms_[ind].s_start_[it].mz_;
        in_quai_[ chroms_[ind].s_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        += chroms_[ind].s_start_[it].intensity_;
        rt_quai_[ chroms_[ind].s_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        = spectra_[chroms_[ind].s_start_[it].spec_id_].get_rt();
        count++;
      }
      it++;
    }
    it = 0;
    count = 0;
    while (count < chroms_[ind].h_count_) {
      if ( chroms_[ind].h_start_[it].chrom_id_ == ind ) {
        mz_quai_[ chroms_[ind].h_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        += chroms_[ind].h_start_[it].mz_;
        in_quai_[ chroms_[ind].h_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        += chroms_[ind].h_start_[it].intensity_;
        rt_quai_[ chroms_[ind].h_start_[it].spec_id_ - chroms_[ind].start_id_ ]
        = spectra_[chroms_[ind].h_start_[it].spec_id_].get_rt();
        count++;
      }
      it++;
    }
    for (int i = 0; i < chroms_[ind].chrom_len_; i++) {
      mz_s[i] = mz_quai_[i];
      intnss[i] = in_quai_[i];
      rts[i] = rt_quai_[i];
    }
  }
} // specToChrom::get_chrom

} // namespace asaristc
