#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>
#include <zlib.h>

#include "base64.h"
#include "specToChrom.h"


#ifndef UNASSIGNED
#define UNASSIGNED -1
#endif

#ifndef SHORT_PEAK
#define SHORT_PEAK -2
#endif



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

  std::cout << "spec_count_ is " << spec_count_ << std::endl;
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
  /* This is the function where we allocate space for points 
   *   sorts points by intensity */
  fill_windows();
  /* form chromatograms */
  determine_chromatograms();
}

void specToChrom::calc_windows(){
  std::cout << "entering specToChrom::calc_windows(" << std::endl;
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
  rt_quai_.resize(wall_count);
  in_quai_.resize(wall_count);
  mz_quai_.resize(wall_count);
  windows_[0] = smallest_mz_;
  std::cout << "wall_count is " << wall_count << std::endl;
  chroms_.resize(wall_count);

  for (int i = 1; i < wall_count; i++) {
    windows_[i] = windows_[i-1] * bot_factor;
  }
  win_to_mzwin_.resize(wall_count, -1);
  std::cout << "leaving specToChrom::calc_windows(" << std::endl;
}

void specToChrom::account_for_points() {
std::cout << "entering specToChrom::account_for_points" << std::endl;
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

std::cout << "before counting big_points_" << std::endl;

  for (int i = 0; i < window_counts.size(); i++) {
    if (window_counts[i]!=0) {
      win_to_mzwin_[i]=n_wind;
      n_wind++;
      big_points_+=window_counts[i];
      max_window_pop_ = std::max(max_window_pop_, window_counts[i]);
    }
  }

std::cout << "after counting big_points_" << std::endl;

  point_windows_.resize(big_points_);
  window_addresses_.resize(big_points_);

  mz_windows_.resize(n_wind);
  win_to_mzwin_.resize(n_wind);
  mzwin_to_win_.resize(n_wind);
  wind_ind = 0;
  int if_counter=0;

std::cout << "before assigning window params" << std::endl;

  for (int i = 0; i < window_counts.size(); i++) {
    if (window_counts[i]!=0) {

      mz_windows_[wind_ind].min_mz_ = windows_[i];
      //std::cout << "windows_[i]" << std::endl;
      mz_windows_[wind_ind].min_ind_ = i;
      //std::cout << "i" << std::endl;
      mz_windows_[wind_ind].population_ = window_counts[i];
      //std::cout << "window_counts[i]" << std::endl;
      mz_windows_[wind_ind].start_ = &point_windows_[point_ofs];

      //std::cout << "&point_windows_[point_ofs]" << std::endl;

      mzwin_to_win_[wind_ind] = i;

      //std::cout << " mzwin_to_win_[wind_ind] = i;" << std::endl;

      wind_ind++;
      point_ofs += window_counts[i];
    }
  }

std::cout << "after assigning window params" << std::endl;
std::cout << "leaving specToChrom::account_for_points" << std::endl;
}

void specToChrom::fill_windows(){
  std::cout << "entering specToChrom::fill_windows" << std::endl;
  /* how many points have we touched yet? */
  double * spec_i_mzs_;
  double * spec_i_intns_;
  int pw_address = 0;
  int specId;
  int wind_ind;
  point * destination;

  std::vector<int> wind_pops(mz_windows_.size());
  std::fill(wind_pops.begin(),wind_pops.end(),0);

  std::cout << "before windows" << std::endl;

  int mz_wind_ind;

  for (int i = 0; i < spectra_.size(); i++) {
    wind_ind=0;
    specId = spectra_[i].get_id();
    spec_i_mzs_ = spectra_[i].get_mzs();
    spec_i_intns_ = spectra_[i].get_intns();
    for (int j = 0; j < spectra_[i].get_n_pts(); j++) {
      if (spec_i_intns_[j] < minimum_intensity_) { continue; }
      while ( spec_i_mzs_[j] > windows_[wind_ind+1] ) {wind_ind++;}

      mz_wind_ind = win_to_mzwin_[wind_ind];

      if (mz_wind_ind==-1) {std::cout << "mz_wind_ind ==0" << std::endl; double saboteur = * (double *)nullptr; }

      destination = &(mz_windows_[ mz_wind_ind  ].start_[ wind_pops[ mz_wind_ind ] ]);

      (*destination).mz_ = spec_i_mzs_[j];
      (*destination).intensity_ = spec_i_intns_[j];
      (*destination).spec_id_ = specId;
      (*destination).chrom_id_ = UNASSIGNED;
      (*destination).window_id_ = wind_ind;

      wind_pops[ mz_wind_ind ]++;

      /*
      point_windows_[pw_address].mz_ = spec_i_mzs_[j];
      point_windows_[pw_address].intensity_ = spec_i_intns_[j];
      point_windows_[pw_address].spec_id_ = specId;
      point_windows_[pw_address].chrom_id_ = UNASSIGNED;
      point_windows_[pw_address].window_id_ = wind_ind;
*/
      //pw_address++;

    }
  } 

  {
  std::ofstream out_putf("/home/joseph/Practice/chromatograms/tryfile.txt");
  for (int i = 0; i < point_windows_.size(); i++) {
    out_putf << point_windows_[i].mz_ << " " << point_windows_[i].intensity_ << std::endl;
  }
  out_putf.close();
  }


  for (int i = 0; i < window_addresses_.size(); i++) {
    window_addresses_[i] = i;
  }
  
  std::sort(window_addresses_.begin(), window_addresses_.end(), [&](int a, int b){
    return point_windows_[a].intensity_ > point_windows_[b].intensity_;
  });


  /*for (int i = 0; i < point_windows_.size(); i++) {
    std::cout << point_windows_[window_addresses_[i]].intensity_ << std::endl;
  }*/

  std::cout << "leaving specToChrom::fill_windows" << std::endl;
}

/* figure out if the neighboring windows are actually neighboring m/z regions */
void specToChrom::neighbor_tables(){
  double bot_factor = 1 / (1 - (2*width_));
  bigger_NT_.resize(mz_windows_.size(), false);
  smaller_NT_.resize(mz_windows_.size(), false);
 
  smaller_NT_[0] = false;
  bigger_NT_[mz_windows_.size()-1] = false;
 
  for (int i = 1; i < mz_windows_.size()-1; i++) {
    smaller_NT_[i] = 
!(std::abs(
(mz_windows_[i-1].min_mz_ * width_* 2 * bot_factor) - mz_windows_[i].min_mz_
) > 1e-11) ;
    bigger_NT_[i] = 
!(std::abs(
(mz_windows_[i].min_mz_ * width_* 2 * bot_factor) - mz_windows_[i+1].min_mz_
)>1e-11);
  }

}

void specToChrom::determine_chromatograms(){
  //std::cout << "entering specToChrom::form_chromatograms" << std::endl;
  double dist;
  double wait;
  double now;
  neighbor_tables();
  int max_points = 0;
  int local_points;
  int mz_window_id;
  int chrom_id = 0;
  mz_window * lower_w = nullptr;
  mz_window * same_w = nullptr;
  mz_window * high_w = nullptr;
  //chrom_count_ = 0;
  /* list of points in the current */
  int assignees_count;
  std::cout << "max_window_pop_ is " << max_window_pop_ << std::endl;
  std::vector<point*> assignees_bag(max_window_pop_ * 3);
  for (int i =0; i < window_addresses_.size(); i++) {
    local_points=0;
    assignees_count=0;
    if ( point_windows_[ window_addresses_[i] ].chrom_id_ !=UNASSIGNED ) {
      continue;
    }

    mz_window_id = win_to_mzwin_[ point_windows_[ window_addresses_[i] ] .window_id_] ;

    same_w = &mz_windows_[ mz_window_id ];

    if (!smaller_NT_[ win_to_mzwin_[point_windows_[window_addresses_[i]].window_id_]] ) {
      lower_w = nullptr;
    } else {
      lower_w = &mz_windows_[ mz_window_id -1 ];
    }

    if (!bigger_NT_[ win_to_mzwin_[point_windows_[window_addresses_[i]].window_id_  ]] ) {
      high_w = nullptr;
    } else {
      high_w = &mz_windows_[ mz_window_id +1 ];
    }

    if (chrom_id > windows_.size()) {
      chroms_.resize(chroms_.size()+windows_.size());
    }

    check_point( point_windows_[ window_addresses_[i] ], local_points, 
                 chrom_id, dist, wait, now, lower_w, same_w, high_w,
                 &assignees_bag[0], chroms_[chrom_id] );

    //std::cout << "chrom_id is " << chrom_id << std::endl;

    max_points = std::max(max_points, local_points);
    if (local_points != SHORT_PEAK) {chrom_id++;}
  }

  std::cout << "chrom_id is " << chrom_id << std::endl;
  print_membership();
  //std::cout << "leaving specToChrom::form_chromatograms" << std::endl;
}

void specToChrom::check_point(point &candidate, int &n_points, int &chrom_id, 
                              double &chrom_dist, double &wait, 
                              double &check_time, mz_window * lower_w,
                              mz_window * same_w, mz_window * high_w,
                              point ** assignees, asari_point &chrom) {
  if ( lower_w == nullptr && high_w == nullptr ) {
    check_point_1w(candidate, n_points, chrom_id, chrom_dist, wait, check_time,
                   same_w, assignees, chrom);
  } else if (lower_w == nullptr || high_w == nullptr) {
    mz_window * other_w = (lower_w == nullptr ? high_w : lower_w);
    check_point_2w(candidate, n_points, chrom_id, chrom_dist, wait, check_time,
                   same_w, other_w,assignees, chrom, 
                   (lower_w == nullptr ? chrom.h_count_ : chrom.l_count_ ),
                   (lower_w == nullptr ? &(chrom.h_start_) : &(chrom.l_start_)));
  } else {
    check_point_3w(candidate, n_points, chrom_id, chrom_dist, wait, check_time,
                   lower_w, same_w, high_w, assignees, chrom);
  }
}

void specToChrom::print_membership() {
  //for (int i = 0; i < point_windows_.size(); i++) {
   // if (point_windows_[i].chrom_id_ < 100 && point_windows_[i].chrom_id_ > -1) {
      //std::cout << "point_windows_[i].chrom_id_ is " << point_windows_[i].chrom_id_ << '\n';
    //}
  //}
  //std::cout << std::endl;
}


void specToChrom::check_point_3w(point &candidate, int &n_points,
                                 int &chrom_id, double &chrom_dist, 
                                 double &wait, double &check_time, 
                                 mz_window * lower_w, mz_window * same_w,
                                 mz_window * high_w, point ** assignees,
                                 asari_point &chrom) {
  
  chrom.l_start_ = nullptr;
  chrom.h_start_ = nullptr;
  chrom.l_count_ = 0;
  chrom.s_count_ = 0;
  chrom.h_count_ = 0;
  chrom.start_id_ = spectra_.size();
  assignees[0] = &candidate;
  n_points = 1;
  int ind;
  int l_ind;
  int h_ind;

  int pop = same_w->population_;
  int l_pop = lower_w->population_;
  int h_pop = high_w->population_;

  int walk_spec = candidate.spec_id_;

  int h_top_id;
  int l_top_id;
  int top_id;

  int h_bot_id;
  int l_bot_id;
  int bot_id;

  int top_spec = candidate.spec_id_;
  int bot_spec = candidate.spec_id_;

  bool walk_c;
  bool l_walk_c;
  bool h_walk_c;

  /* we can't do out hacky idea from before. we have to walk all the way up
   * on both windows */
  for (int i = 0; i < same_w->population_; i++ ) {
    if (std::abs( same_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2) {
      ind = i;
      break;
    }
  }

  top_id = ind;
  bot_id = ind;

  for (int i = 0; i < lower_w->population_; i++ ) {
    if ( std::abs( lower_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2 ) {
      l_ind = i;
      break;
    }
  }

  l_top_id = l_ind;
  l_bot_id = l_ind;

  for (int i = 0; i < high_w->population_; i++ ) {
    if ( std::abs( high_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2 ) {
      h_ind = i;
      break;
    }
  }

  h_top_id = h_ind;
  h_bot_id = h_ind;

  walk_c  = true;
  l_walk_c = true;
  h_walk_c = true;

  while (walk_c || l_walk_c || h_walk_c) {
    walk_c  = false;
    l_walk_c = false;
    h_walk_c = false;
    for (int i = top_id; i < pop; i++) {
      if (same_w->start_[i].spec_id_ == top_spec+1 ) {
        top_id = i;
        break;
      }
    }
    for (int i = l_top_id; i < l_pop; i++) {
      if (lower_w->start_[i].spec_id_ == top_spec+1 ) {
        l_top_id = i;
        break;
      }
    }
    for (int i = h_top_id; i < h_pop; i++) {
      if (high_w->start_[i].spec_id_ == top_spec+1 ) {
        h_top_id = i;
        break;
      }
    }
    while ( top_id < pop && same_w->start_[top_id].spec_id_ - walk_spec < 2) {
      if (std::abs(same_w->start_[top_id].mz_ - candidate.mz_ ) < width_ * candidate.mz_) {
        top_spec = same_w->start_[top_id].spec_id_;
        //same_w->start_[top_id].chrom_id_ = chrom_id;
        assignees[n_points]= &(same_w->start_[top_id]);
        chrom.s_count_++;
        n_points++;
        walk_c = true;
      }
      top_id++;
    }
    while (l_top_id < l_pop&&lower_w->start_[l_top_id].spec_id_ - walk_spec < 2) {
      if (std::abs(lower_w->start_[l_top_id].mz_ - candidate.mz_ ) < width_ * candidate.mz_) {
        top_spec = lower_w->start_[l_top_id].spec_id_;
        //lower_w->start_[l_top_id].chrom_id_ = chrom_id;
        assignees[n_points] = &(lower_w->start_[l_top_id]);
        n_points++;
        chrom.l_count_++;
        l_walk_c = true;
      }
      l_top_id++;
    }
    while (h_top_id < h_pop&&high_w->start_[h_top_id].spec_id_ - walk_spec < 2) {
      if (std::abs(high_w->start_[h_top_id].mz_ - candidate.mz_ ) < width_ * candidate.mz_) {
        top_spec = high_w->start_[h_top_id].spec_id_;
        //high_w->start_[h_top_id].chrom_id_ = chrom_id;
        assignees[n_points] = &(high_w->start_[h_top_id]);
        chrom.h_count_++;
        n_points++;
        h_walk_c = true;
      }
      h_top_id++;
    }
    walk_spec = top_spec;
  }

  walk_c  = true;
  l_walk_c = true;
  h_walk_c = true;

  walk_spec = bot_spec;

  // loop over the bot inds
  while (walk_c || l_walk_c || h_walk_c) {
    walk_c  = false;
    l_walk_c = false;
    h_walk_c = false;
    for (int i = bot_id; i > 0; i--) {
      if (same_w->start_[i].spec_id_ == bot_spec -1) {
        bot_id = i;
        break;
      }
    }
    for (int i = l_bot_id; i > 0; i--) {
      if (lower_w->start_[i].spec_id_ == bot_spec -1) {
        l_bot_id = i;
        break;
      }
    }
    for (int i = h_bot_id; i > 0; i--) {
      if (high_w->start_[i].spec_id_ == bot_spec -1) {
        h_bot_id = i;
        break;
      }
    }
    while (bot_id > 0 && walk_spec - same_w->start_[bot_id].spec_id_ < 2 ) {
      if (std::abs( same_w->start_[bot_id].mz_ - candidate.mz_ )< width_ * candidate.mz_) {
        bot_spec = same_w->start_[bot_id].spec_id_;
        //same_w->start_[bot_id].chrom_id_ = chrom_id;
        chrom.s_count_++;
        chrom.s_start_ = &(same_w->start_[bot_id]);
        assignees[n_points] = &(same_w->start_[bot_id]);
        n_points++;
      }
      bot_id--;
    }
    while (l_bot_id > 0 && walk_spec - lower_w->start_[l_bot_id].spec_id_ < 2) {
      if (std::abs(lower_w->start_[l_bot_id].mz_ - candidate.mz_ )< width_ * candidate.mz_) {
        bot_spec = lower_w->start_[l_bot_id].spec_id_;
        chrom.l_count_++;
        //lower_w->start_[l_bot_id].chrom_id_ = chrom_id;
        chrom.l_start_ = &(lower_w->start_[l_bot_id]);
        assignees[n_points] = &(lower_w->start_[l_bot_id]);
        n_points++;
      }
      l_bot_id--;
    }
    while (h_bot_id > 0 && walk_spec - high_w->start_[h_bot_id].spec_id_ < 2) {
      if (std::abs(high_w->start_[h_bot_id].mz_ - candidate.mz_ )< width_ * candidate.mz_) {
        bot_spec = high_w->start_[h_bot_id].spec_id_;
        //high_w->start_[h_bot_id].chrom_id_ = chrom_id;
        chrom.h_start_ = &(high_w->start_[h_bot_id]);
        chrom.h_count_++;
        assignees[n_points] =&(high_w->start_[h_bot_id]);
        n_points++;
      }
      h_bot_id--;
    }
    walk_spec = bot_spec;
  }
  if (top_spec - bot_spec < min_steps_) {
    for (int i = 0; i < n_points; i++) {
      assignees[i]->chrom_id_ = SHORT_PEAK;
    }
    n_points=SHORT_PEAK;
  } else {
    for (int i = 0; i < n_points; i++) {
      assignees[i]->chrom_id_ = chrom_id;
    }
    chrom.chrom_len_ = top_spec - bot_spec;
  }
}//specToChrom::check_point_3w

void specToChrom::check_point_2w(point &candidate, int &n_points,
                                 int &chrom_id, double &chrom_dist, 
                                 double &wait, double &check_time,
                                 mz_window * same_w, mz_window * other_w,
                                 point ** assignees, asari_point &chrom,
                                 int &o_count, point ** o_start) {
  chrom.l_start_ = nullptr;
  chrom.h_start_ = nullptr;
  chrom.l_count_ = 0;
  chrom.s_count_ = 0;
  chrom.h_count_ = 0;
  chrom.start_id_ = spectra_.size();
  assignees[0] = &candidate;
  n_points = 1;
  int ind;
  int o_ind;

  int pop = same_w->population_;
  int o_pop = other_w->population_;

  int walk_spec = candidate.spec_id_;

  int o_top_id;
  int o_bot_id;
  int top_id;
  int bot_id;

  int top_spec = candidate.spec_id_;//int o_top_spec;
  int bot_spec = candidate.spec_id_;//int o_bot_spec;
  /* walk control variables */
  bool walk_c;
  bool o_walk_c;

  /* we can't do out hacky idea from before. we have to walk all the way up
   * on both windows */
  for (int i = 0; i < same_w->population_; i++ ) {
    if (std::abs( same_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2) {
      ind = i;
      break;
    }
  }

  top_id = ind;
  bot_id = ind;

  for (int i = 0; i < other_w->population_; i++ ) {
    if ( std::abs( other_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2 ) {
      o_ind = i;
      break;
    }
  }

  o_top_id = o_ind;
  o_bot_id = o_ind;
  walk_c  = true;
  o_walk_c = true;
  // loop over the top spectra
  // Every entry point to t
  while (walk_c || o_walk_c) {
    walk_c  = false;
    o_walk_c = false;
    for (int i = top_id; i < pop; i++) {
      if (same_w->start_[i].spec_id_ == top_spec+1 ) {
        top_id = i;
        break;
      }
    }
    for (int i = o_top_id; i < o_pop; i++) {
      if (other_w->start_[i].spec_id_ == top_spec+1 ) {
        o_top_id = i;
        break;
      }
    }
    while ( top_id < pop && same_w->start_[top_id].spec_id_ - walk_spec < 2) {
      if (std::abs(same_w->start_[top_id].mz_ - candidate.mz_ ) < width_ * candidate.mz_) {
        top_spec = same_w->start_[top_id].spec_id_;
        //same_w->start_[top_id].chrom_id_ = chrom_id;
        chrom.s_count_++;
        assignees[n_points] = &(same_w->start_[top_id]);
        n_points++;
        walk_c = true;
      }
      top_id++;
    }
    while (o_top_id < o_pop&&other_w->start_[o_top_id].spec_id_ - walk_spec < 2) {
      if (std::abs(other_w->start_[o_top_id].mz_ - candidate.mz_ ) < width_ * candidate.mz_) {
        top_spec = other_w->start_[o_top_id].spec_id_;
        //other_w->start_[o_top_id].chrom_id_ = chrom_id;
        assignees[n_points] = &(other_w->start_[o_top_id]);
        n_points++;
        o_count++;
        o_walk_c = true;
      }
      o_top_id++;
    }
    walk_spec = top_spec;
  }

  walk_c  = true;
  o_walk_c = true;

  walk_spec = bot_spec;

  // loop over the bot inds
  while (walk_c || o_walk_c) {
    walk_c  = false;
    o_walk_c = false;
    for (int i = bot_id; i > 0; i--) {
      if (same_w->start_[i].spec_id_ == bot_spec -1) {
        bot_id = i;
        break;
      }
    }
    for (int i = o_bot_id; i > 0; i--) {
      if (other_w->start_[i].spec_id_ == bot_spec -1) {
        o_bot_id = i;
        break;
      }
    }
    while (bot_id > 0 && walk_spec - same_w->start_[bot_id].spec_id_ < 2 ) {
      if (std::abs( same_w->start_[bot_id].mz_ - candidate.mz_ )< width_ * candidate.mz_) {
        bot_spec = same_w->start_[bot_id].spec_id_;
        //same_w->start_[bot_id].chrom_id_ = chrom_id;
        assignees[n_points] = &(same_w->start_[bot_id]);
        chrom.s_start_ = &(same_w->start_[bot_id]);
        chrom.s_count_++;
        n_points++;
      }
      bot_id--;
    }
    while (o_bot_id > 0 && walk_spec - other_w->start_[o_bot_id].spec_id_ < 2) {
      if (std::abs(other_w->start_[o_bot_id].mz_ - candidate.mz_ )< width_ * candidate.mz_) {
        bot_spec = other_w->start_[o_bot_id].spec_id_;
        //other_w->start_[o_bot_id].chrom_id_ = chrom_id;
        assignees[n_points] = &(other_w->start_[o_bot_id]);
        *o_start = &(other_w->start_[o_bot_id]);
        o_count++;
        n_points++;
      }
      bot_id--;
    }
    walk_spec = bot_spec;
  }
  if (top_spec - bot_spec < min_steps_) {
    for (int i = 0; i < n_points; i++) {
      assignees[i]->chrom_id_ = SHORT_PEAK;
    }
    n_points=SHORT_PEAK;
    //std::cout << "Assigning SHORT_PEAK for length: " << top_id - bot_id << std::endl;
  } else {
    for (int i = 0; i < n_points; i++) {
      assignees[i]->chrom_id_ = chrom_id;
    }
    chrom.chrom_len_ = top_spec - bot_spec;
    //std::cout << "Assigning chrom_id for length: " << top_id - bot_id << std::endl;
  }
}// specToChrom::check_point_2w

void specToChrom::check_point_1w(point &candidate, int &n_points,
                                 int &chrom_id, /* gets reset with the id for 
                                                   the chromatogram */
                                 double &chrom_dist, /* can't see use */
                                 double &wait,
                                 double &check_time, 
                                 mz_window * same_w,
                                 point ** assignees,
                                 asari_point &chrom /*,
                                 mz_window * lower_w, 
                                 mz_window * high_w*/ ) {
  chrom.l_start_ = nullptr;
  chrom.h_start_ = nullptr;
  chrom.l_count_ = 0;
  chrom.s_count_ = 0;
  chrom.h_count_ = 0;
  chrom.start_id_ = spectra_.size();


  assignees[0] = &candidate;
  n_points = 1;
  int top_id;
  int bot_id;
  for (int i = 0; i < same_w->population_; i++) {
    if ( std::abs(  same_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2 ) {
      top_id = same_w->start_[i].spec_id_;
      bot_id = same_w->start_[i].spec_id_;
      /* walk up */
      for (int j = i; j < same_w->population_; j++ ) {
        if ( std::abs(top_id - same_w->start_[j].spec_id_ ) > 1 ) {break;}
        if (std::abs(same_w->start_[j].mz_ - candidate.mz_) > width_ * candidate.mz_) {
          continue;
        } else {
          assignees[n_points] = &(same_w->start_[j]);
          n_points++;
          chrom.s_count_++;
          //same_w->start_[j].chrom_id_ = chrom_id;
          top_id = same_w->start_[j].spec_id_;
        }
      }
      /* walk down */
      for (int j = i; j >=0; j--) {
        if ( std::abs(bot_id - same_w->start_[j].spec_id_ ) > 1  ) {break;}
        if (std::abs(same_w->start_[j].mz_ - candidate.mz_) > width_ * candidate.mz_) {
          continue;
        } else {
          assignees[n_points] = &(same_w->start_[j]);
          chrom.s_count_++;
          chrom.start_id_ = std::min(same_w->start_[j].spec_id_, 
                                     chrom.start_id_);
          n_points++;
          //same_w->start_[j].chrom_id_ = chrom_id;
          bot_id = same_w->start_[j].spec_id_;
          chrom.s_start_ = &(same_w->start_[j]);
        }
      }
      break;
    }
  }
  if (top_id - bot_id < min_steps_) {
    for (int i = 0; i < n_points; i++) {
      assignees[i]->chrom_id_ = SHORT_PEAK;
    }
    n_points = SHORT_PEAK;
    //std::cout << "assigning SHORT_PEAK for length: " << top_id - bot_id << std::endl;
  } else {
    for (int i = 0; i < n_points; i++) {
      assignees[i]->chrom_id_ = chrom_id;
    }
    chrom.chrom_len_ = top_id - bot_id;
    //std::cout << "assigning chom_id for length: " << top_id - bot_id << std::endl;
  }
}// specToChrom::check_point_1w


/* bit of a misnomer. This function looks in the area around a point to decide
 * which points around it belong in its chromatogram. 
 * crucially, it sets those chromatograms to  */
void specToChrom::othercheck_point(point &candidate, int &n_points, 
                              int &chrom_id, double &chrom_dist, double &wait,
                              double &check_time, mz_window * lower_w,
                              mz_window * same_w, mz_window * high_w) {
std::cout << "entering specToChrom::check_point" << std::endl;
  n_points = 1;

  double current_length_l;
  int added_per_spec;
  /* spec_id's we're going to walk along while examining the lower neighbor */
  int curr_ind_l = 0;  int top_ind_l;  int bot_ind_l;

  double current_length_h;  int top_ind_h;  int bot_ind_h;

  int top_ind;  int bot_ind;

  check_time = spectra_[candidate.spec_id_].get_rt();
  if (lower_w != nullptr) {
    for (int i = 0; i < lower_w->population_; i++) {
      if ( std::abs(  lower_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2 ) {
        bot_ind_l = lower_w->start_[i].spec_id_;
        top_ind_l = lower_w->start_[i].spec_id_;
 
	/* walk along preceding points to find chromatogram members */
        for ( int j = i; j >= 0; j-- ) {
	  /* end search if not in a neighboring spectrum */
          if ( (bot_ind_l - lower_w->start_[j].spec_id_  ) > 1 ) {
            break;
} else if ( std::abs(lower_w->start_[j].mz_ - candidate.mz_) > width_ ) {
	    continue;
          } else {
            n_points++;
            lower_w->start_[j].chrom_id_ = chrom_id;
            bot_ind_l = lower_w->start_[j].spec_id_;
	      }
        }

        /* walk along following points to find chromatogram members */
        for (int j = i; j < lower_w->population_; j++) {
	      if ( (lower_w->start_[j].spec_id_ - top_ind_l) > 1 ) {
	        break;
} else if ( std::abs(lower_w->start_[j].mz_ - candidate.mz_) > width_  ) {
            continue;
	      } else {
            n_points++;
            lower_w->start_[j].chrom_id_ = chrom_id;
            top_ind_l = lower_w->start_[j].spec_id_;
          }
	    }
        break;
      }
    }
  } 
  if (high_w != nullptr) {
    for (int i = 0; i < high_w->population_; i++) {
      if ( std::abs(  high_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2 ) {
        bot_ind_h = high_w->start_[i].spec_id_;
        top_ind_h = high_w->start_[i].spec_id_;
 
	/* walk along preceding points to find chromatogram members */
        for ( int j = i; j >= 0; j-- ) {
	  /* end search if not in a neighboring spectrum */
          if ( (bot_ind_h - high_w->start_[j].spec_id_  ) > 1 ) {
            break;
} else if ( std::abs(high_w->start_[j].mz_ - candidate.mz_) > width_ ) {
	    continue;
          } else {
            n_points++;
            high_w->start_[j].chrom_id_ = chrom_id;
            bot_ind_h = high_w->start_[j].spec_id_;
	      }
        }

        /* walk along following points to find chromatogram members */
        for (int j = i; j < high_w->population_; j++) {
	      if ( (high_w->start_[j].spec_id_ - top_ind_l) > 1 ) {
	        break;
} else if ( std::abs(high_w->start_[j].mz_ - candidate.mz_) > width_  ) {
            continue;
	      } else {
            n_points++;
            high_w->start_[j].chrom_id_ = chrom_id;
            top_ind_h = high_w->start_[j].spec_id_;
          }
	    }
        break;
      }
    }

  }
  

  for ( int i = 0; i < same_w->population_; i++ ) {
    if (std::abs( same_w->start_[i].spec_id_ - candidate.spec_id_ ) < 2 ) {

      top_ind = same_w->start_[i].spec_id_;
      bot_ind = same_w->start_[i].spec_id_;

      for (int j = i; j >= 0; j--) {
        if ( (bot_ind - same_w->start_[j].spec_id_ ) > 1 ) {
          break;
} else if ( std::abs(same_w->start_[j].mz_ - candidate.mz_ ) > width_ ) {
          continue;
        } else {
          n_points++;
          same_w->start_[j].chrom_id_ = chrom_id;
          bot_ind = same_w->start_[j].spec_id_;
        }
      }

      for ( int j = i; j < same_w->population_; j++ ) {
        if ( same_w->start_[j].spec_id_  - top_ind ) {
          break;
} else if (std::abs(same_w->start_[j].mz_ - candidate.mz_ ) > width_) {
          continue;
        } else {
          n_points++;
          same_w->start_[j].chrom_id_ = chrom_id;
          top_ind = same_w->start_[j].spec_id_;
        }
      }

      break;
    }
  }
 
}

void specToChrom::writeChromatograms(std::string fname){
  
}

void specToChrom::print_filename () {
  std::cout << readFilename_.c_str() << std::endl;
}

void specToChrom::reset() {
  lcms_DOC_.clear();
  chrom_doc_.clear();
  parsedMzML_.clear();
  spectra_.clear();
  
  windows_.clear();
  mz_windows_.clear();
  point_windows_.clear();
  window_addresses_.clear();

  win_to_mzwin_.clear();
  mzwin_to_win_.clear();

  bigger_NT_.clear();
  smaller_NT_.clear();

  chroms_.clear();

  rt_quai_.clear();
  in_quai_.clear();
  mz_quai_.clear();

  current_spec_ = nullptr;
}

} // namespace asaristc
