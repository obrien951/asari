#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <zlib.h>

#include "base64.h"
#include "specToChrom.h"

namespace asaristc {

specToChrom::specToChrom() {
  current_spec_ = nullptr;
  text_decoder_ = std::make_shared<b64_decoder>();

  tmp_mz_.resize(1048576);
  tmp_intns_.resize(1048576);
}

void specToChrom::parse_xml() {
 std::ifstream mzmlFile(readFilename_.c_str());

  /* Get the text from the file. We could use vector.end() to get the same
   * result afaik */
  std::cout << readFilename_.c_str() << std::endl;
  parsedMzML_.insert(parsedMzML_.begin(),
                     std::istreambuf_iterator<char>(mzmlFile),
                     std::istreambuf_iterator<char>());

  parsedMzML_.push_back('\0');

  lcms_DOC_.parse<0>(&parsedMzML_[0]);
  spec_count_ = std::stoi(lcms_DOC_.first_node()
                              ->first_node()
                              ->first_node("run")
                              ->first_node("spectrumList")
                              ->first_attribute("count")
                              ->value());

  std::cout << spec_count_ << std::endl;
  mzmlFile.close(); 
}

void specToChrom::set_filename(std::string filename) {
  readFilename_ = filename;
}

void specToChrom::readSpectra() {
  parse_xml();
  initiate_spectra();
  convert_spectra();
}

void specToChrom::initiate_spectra(){
  spectra_.reserve(spec_count_);
  int npts;
  double rt;
  int id;
  int ofst = 0;
  for (int i = 0; i < spec_count_; i++) {
    advance_run();
    xmlRtntnTime(rt);
    xmlSpecPts(npts);
    if (i==0) {
      std::cout << "rt is " << rt << std::endl;
      std::cout << "npts is " << npts << std::endl;
    }
    spectra_.push_back(spectrum(npts, rt, i, ofst));
    ofst += npts;
  }
  mzs_.resize(ofst);
  intns_.resize(ofst);
  for (int i = 0; i < spec_count_; i++) {
    spectra_[i].pointers_to_offset(&mzs_[0], &intns_[0]);
  }
  specRunToFront();
}

void specToChrom::convert_spectra(){
  double* mz;
  double* intns;
  for (int i = 0; i < spec_count_; i++) {
    advance_run();
    mz = spectra_[i].get_mzs();
    intns = spectra_[i].get_intns();
    spec_data_from_b64(mz, intns);
  }
  std::cout << "smallest_mz_ is " << smallest_mz_ << std::endl;
  std::cout << "biggest_mz_ is " << biggest_mz_ << std::endl;
  //check_spec(0);
}// specToChrom::convert_spectra

void specToChrom::check_spec(int spec_id){
  double * check_vals = spectra_[spec_id].get_mzs();
  for (int i = 0; i < spectra_[spec_id].get_n_pts(); i++){
    std::cout << "i is " << i << "check_vals is " << check_vals[i] << " ";
  }
  std::cout << std::endl;
}

void specToChrom::spec_data_from_b64(double* mz, double* intns) {
  set_mz_intns();
  decoding_step();
  min_max_mz();
  decompression_step(mz, intns);
}// void specToChrom::spec_data_from_b64

void specToChrom::set_mz_intns() {
  current_mz_ = nullptr;
  current_intns_ = nullptr;
  for (rapidxml::xml_node<> * ba_it = current_spec_
         ->first_node("binaryDataArrayList")
         ->first_node("binaryDataArray"); 
         ba_it; ba_it = ba_it->next_sibling("binaryDataArray") ) {
    hold_fpe_ = -1;
    hold_zlc_ = -1;
    /* loop over the controlled vocabulary for the current  */
    for (rapidxml::xml_node<> * cv_it = ba_it->first_node("cvParam");
         cv_it; cv_it = cv_it->next_sibling("cvParam") ) {
      char * acc = cv_it->first_attribute("accession")->value();
      switch( std::stoi( &acc[3] ) ) {
      case 1000523: 
        hold_fpe_ = 64;
        break;
      case 1000521: 
        hold_fpe_ = 32;
        break;  
      case 1000576: 
        hold_zlc_ = 0;
        break;
      case 1000574: 
        hold_zlc_ = 1;
        break;
      case 1000514:
        current_mz_ = ba_it;
        break;
      case 1000515:
        current_intns_ = ba_it;
        break;
      }
    }
    if (current_mz_ == ba_it) {
      current_mz_el_ = std::stoi(current_mz_->first_attribute("encodedLength")->value());
      current_m_z_fpe_ = hold_fpe_;
      current_m_z_zlc_ = hold_zlc_;
      current_m_z_data_ = current_mz_->first_node("binary")->value();
    }
    if (current_intns_ == ba_it) {
      current_intensity_el_ = std::stoi(current_intns_->first_attribute("encodedLength")->value());
      current_intensity_fpe_ = hold_fpe_;
      current_intensity_zlc_ = hold_zlc_;
      current_intensity_data_ = current_intns_->first_node("binary")->value();
    }
  }
}//specToChrom::set_mz_intns

void specToChrom::decoding_step(){
  post_base64_m_z_size_ =
    static_cast<long unsigned int>((current_m_z_data_[current_mz_el_ - 1] == '='
                                      ? (current_m_z_data_[current_mz_el_ - 2] == '='
                                           ? (3 * (current_mz_el_ - 4)) / 4 + 1
                                           : (3 * (current_mz_el_ - 4)) / 4 + 2)
                                      : (3 * current_mz_el_) / 4));
  post_base64_intensity_size_ =
    static_cast<long unsigned int>((current_intensity_data_[current_intensity_el_ - 1] == '='
                                      ? (current_intensity_data_[current_intensity_el_ - 2] == '='
                                           ? (3 * (current_intensity_el_ - 4)) / 4 + 1
                                           : (3 * (current_intensity_el_ - 4)) / 4 + 2)
                                      : (3 * current_intensity_el_) / 4));
  
  text_decoder_->decode_base64(&tmp_mz_[0], current_m_z_data_,
                              post_base64_m_z_size_);
  text_decoder_->decode_base64(&tmp_intns_[0], current_intensity_data_,
                              post_base64_intensity_size_);
}//specToChrom::decoding_step

void specToChrom::decompression_step(double * mz, double * intns){
  post_zlib_m_z_size_ = 1048576UL * sizeof(double);
  post_zlib_intensity_size_ = 1048576UL * sizeof(double);
  
  MZSTAT_ =
    uncompress((unsigned char *)&mz[0], &post_zlib_m_z_size_,
               (unsigned char *)&tmp_mz_[0], post_base64_m_z_size_);
  if (MZSTAT_ == Z_BUF_ERROR) {
    std::cout << "Z_BUF_ERROR. data not decoded" << std::endl;
  }

/*  tmp_mz_.resize(1048576);
    tmp_intns_.resize(1048576); */
  
  if (current_m_z_fpe_==32) {
    memcpy( &hold_floats_[0], &mz[0], post_zlib_m_z_size_ );
    for (size_t i = 0UL; i < (post_zlib_m_z_size_/sizeof(float)); i++ ) {
      mz[i] = static_cast<double>(hold_floats_[i]);
    }
    post_zlib_m_z_size_*=2;
  }
  
  INSTAT_ = uncompress(
    (unsigned char *)&intns[0], &post_zlib_intensity_size_,
    (unsigned char *)&tmp_intns_[0], post_base64_intensity_size_);
  if (INSTAT_ == Z_BUF_ERROR) {
    std::cout << "Z_BUF_ERROR. data not decoded" << std::endl;
  }
  
  if ( current_intensity_fpe_==32) {
    memcpy( &hold_floats_[0], &intns[0], post_zlib_intensity_size_ );
    for (size_t i = 0UL; i < (post_zlib_intensity_size_/sizeof(float)); i++ ) {
      intns[i] = static_cast<double>(hold_floats_[i]);
    }
    post_zlib_intensity_size_*=2;
  }
}//specToChrom::decompression_step

void specToChrom::advance_run() {
  if (current_spec_ == nullptr) {
    current_spec_ = lcms_DOC_.first_node()
                        ->first_node()
                        ->first_node("run")
                        ->first_node("spectrumList")
                        ->first_node("spectrum");
  } else if (std::stoi(current_spec_->first_attribute("index")->value()) ==
             (spec_count_ - 1)) {
    current_spec_ = nullptr;
    return;
  } else {
    current_spec_ = current_spec_->next_sibling();
  }
}

void specToChrom::specRunToFront() {
  current_spec_ = nullptr;
}

void specToChrom::findChromatograms(){
  calc_windows();
  /* count the points for each window in each spectrum */
  account_for_points();
  /* This is the function where we allocate space for points */
  fill_windows();
}

void specToChrom::calc_windows(){
  /* first count the number of windows */
  double current_wall = smallest_mz_;
  int wall_count = 1;
  double bot_factor = 1 / (1 - (2*width_));
  std::cout << "bot_factor is " << bot_factor << std::endl;
  while (current_wall < biggest_mz_){
    current_wall  = current_wall * bot_factor;
    wall_count++;
    //std::cout << current_wall << std::endl;
  }

  windows_.resize(wall_count);
  windows_[0] = smallest_mz_;

  for (int i = 1; i < wall_count; i++) {
    windows_[i] = windows_[i-1] * bot_factor;
  }
  for (int i = 1; i < wall_count; i++) {
    if (windows_[i] < 70.0) {
      std::cout << "i is " << i << "windows_[i] is " << windows_[i] << std::endl;
    }
  }
}

void specToChrom::account_for_points() {
  for (int i = 0; i < windows_.size(); i++) {
    if (windows_[i] < 70.0) {
      std::cout << windows_[i] << std::endl;
    }
  }
  int wind_ind = 0;
  int n_wind = 0;
  int point_ofs = 0;
  std::vector<int> window_counts;
  /* the last element will ALWAYS be 0*/
  big_points_=0;
  window_counts.resize(windows_.size(),0);
  //std::cout << "wondow_counts.size() is " << window_counts.size() << std::endl;
  double * spec_i_mzs_;
  double * spec_i_intns_;

  for (int i = 0; i < spectra_.size(); i++) {
    spec_i_mzs_ = spectra_[i].get_mzs();
    spec_i_intns_ = spectra_[i].get_intns();
    for (int j = 0; j < spectra_[i].get_n_pts(); j++) {
      if (spec_i_intns_[j] < minimum_intensity_) {
        continue;
      }
      /* climb up to the right window. the for loop overshoots by 1,
         so count it back down */
      while ( spec_i_mzs_[j] > windows_[wind_ind+1] ) {wind_ind++;}
      
      window_counts[wind_ind]++;
    }
    wind_ind = 0;
  }
  std::cout << "before counting windows" << std::endl;
  for (int i = 0; i < window_counts.size(); i++) {
    if (window_counts[i]!=0) {
      n_wind++;
      big_points_+=window_counts[i];
    }
  }
  mz_windows_.resize(n_wind);
  win_to_mzwin_.resize(n_wind);
  wind_ind = 0;
  for (int i = 0; i < window_counts.size(); i++) {
    if (window_counts[i]!=0) {
      mz_windows_[wind_ind].min_mz_ = windows_[i];
      mz_windows_[wind_ind].min_ind_ = i;
      mz_windows_[wind_ind].population_ = window_counts[i];
      mz_windows_[wind_ind].start_ = &point_windows_[point_ofs];

      win_to_mzwin_[wind_ind] = i;

      wind_ind++;
      point_ofs += window_counts[i];
    }
  }

}

void specToChrom::fill_windows(){
  /* how many points have we touched yet? */
  std::cout << "inside specToChrom::fill_windows" << std::endl;
  double * spec_i_mzs_;
  double * spec_i_intns_;
  int pw_address = 0;
  int specId;
  int wind_ind;
  point_windows_.resize(big_points_);

  for (int i = 0; i < spectra_.size(); i++) {
    wind_ind=0;
    specId = spectra_[i].get_id();
    for (int j = 0; j < spectra_[i].get_n_pts(); j++) {
      if (spec_i_intns_[j] < minimum_intensity_) { continue; }
      while ( spec_i_mzs_[j] > windows_[wind_ind+1] ) {wind_ind++;}
      point_windows_[pw_address].mz_ = spec_i_mzs_[j];
      point_windows_[pw_address].intensity_ = spec_i_intns_[j];
      point_windows_[pw_address].spec_id_ = specId;
      point_windows_[pw_address].chrom_id_ = -1;
      if ( &spec_i_mzs_[j] > &mzs_[ mzs_.size()] ) {
        std::cout << "mzs_ overwritten" << std::endl;
      }
      pw_address++;
      if (pw_address > big_points_) {
        std::cout << "pw violates" << std::endl;
      }
    }
  }
}

void specToChrom::writeChromatograms(std::string fname){
  point * chek_ptr;
  chek_ptr = mz_windows_[0].start_;
  std::cout << chek_ptr[mz_windows_[0].population_ - 1].mz_ << " " << chek_ptr[mz_windows_[0].population_ - 1].spec_id_  << std::endl;
  chek_ptr = mz_windows_[1].start_;
  std::cout << chek_ptr[0].mz_ << chek_ptr[0].spec_id_ << std::endl;
}

void specToChrom::print_filename () {
  std::cout << readFilename_.c_str() << std::endl;
}

void specToChrom::reset() {
  lcms_DOC_.clear();
  chrom_doc_.clear();
  parsedMzML_.clear();
  spectra_.clear();
  current_spec_ = nullptr;
}

} // namespace asaristc
