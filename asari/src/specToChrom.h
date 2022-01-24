#ifndef __ASARI_SPECTOCHROM__
#define __ASARI_SPECTOCHROM__

#define UNASSIGNED -1

#include <vector>
#include <cstring>
#include <string>
#include <memory>

#include "rapidxml/rapidxml.hpp"

#include "base64.h"

namespace asaristc {

struct point {
  double mz_;
  double intensity_;
  int spec_id_;
  int chrom_id_;
  int window_id_;
};

struct asari_point {
  double mz_;
  double rt_;
  double intn_;
};

struct mz_window {
  double min_mz_;
  /* index of min_mz_ in windows_ */
  int min_ind_;
  int population_;
  point * start_;  
};

class spectrum {
public:
  spectrum(int &n_pts, double &RT, int &id, int &offset);
  void copy_values(int &count, double* mzs, double * intns);
  //sets the classes data pointers to the address study_start + offset_
  void pointers_to_offset(double* mz_start, double * intns_start);

  int get_id() {return id_;}

  int get_n_pts() {return n_pts_;}

  double get_rt() {return rt_;}

  double* get_mzs() {return mzs_;}
  double* get_intns() {return intns_;}

protected:
  /* number of points */
  int n_pts_;
  /* spectrum id. used for indexing later */
  int id_;
  /* index of the 0th rt or mz in the array for the entire run*/
  int offset_;
  /* spectrum retention time*/
  double rt_;
  /* location of this spectrum's m/z values in the global array */
  double* mzs_;
  /* location of this spectrum's intensity values in the global array */
  double* intns_;
};

class chromatogram {
public:
  chromatogram(int &count, point &peak, double * mzs, double * rts, double * intns);
protected:
  double mz_;
  std::vector<asari_point> data;
};

// reads spectra from an mzML file then writes them to a chromatogram file
class specToChrom {
public:
  // constructor
  specToChrom();
  /* set the name of the mzml file used for data input. 
   * EXPORTED TO PYTHON */
  void set_filename(std::string filename);

  void set_minimum_intensity(double minimum_intensity) {
    minimum_intensity_ = minimum_intensity;
  }
  /* Print the name of the mzml file to screen
   * EXPORTED TO PYTHON */
  void print_filename();
  /* Free memory in containers and set pointers to NULL
   * Always call before opening a new mzML file
   * EXPORTED TO PYTHON */
  void reset();
  bool is_set() {return true;}
  /* Read mzML data from disk
   * convert data from xml document to regular datatypes
   * EXPORTED TO PYTHON */
  void readSpectra();
  /* process data from spectrum array to mass traces */
  void findChromatograms();
  /* write chromatograms to a rapidxml file*/
  void writeChromatograms(std::string fname);
protected:
  void parse_xml();
  void initiate_spectra();
  void convert_spectra();
  void advance_run();

  /*functions to convert m/z and intensity arrays from base64 
   *spec data from b_64 puts spectra data from the mzml file into the 
   *spectra objects. it uses 
   *set_mz_intns **gets metadata from the spectrum to aid with text to binary **
   *decoding_step **does text to binary conversion**
   *decompression_step **does zlib to decompressed conversion** */
  void spec_data_from_b64(double* mz, double* intns);

  void set_mz_intns();
  void decoding_step();
  void decompression_step(double* mz, double* intns);
  void check_spec(int spec_id);

  /* determine the necessary window boundaries for this lcms run */
  void calc_windows();
  void account_for_points();
  void fill_windows();
  void determine_chromatograms();

  void print_membership();

  void neighbor_tables();


  void check_point_2w(point &candidate, int &n_points,
                                 int &chrom_id, double &chrom_dist, 
                                 double &wait, double &check_time,
                                 mz_window * same_w, mz_window * other_w,
                                 point ** assignees);

  void check_point_3w(point &candidate, int &n_points,
                                 int &chrom_id, double &chrom_dist,
                                 double &wait, double &check_time,
                                 mz_window * lower_w, mz_window * same_w,
                                 mz_window * high_w, point ** assignees);


  void check_point_1w(point &candidate, int &n_points, int &chrom_id, 
                      double &chrom_dist, double &wait, double &check_time, 
                      mz_window * same_w, point ** assignees);

  void othercheck_point(point &candidate, int &n_points, int &chrom_id, 
                   double &chrom_dist, double &wait, double &check_time,
                   mz_window * lower_w, mz_window * same_w, mz_window * high_w);

  void check_point(point &candidate, int &n_points, int &chrom_id, 
                   double &chrom_dist, double &wait, double &check_time,
                   mz_window * lower_w, mz_window * same_w, mz_window * high_w,
                   point ** assignees);

  /* reset linked list traversal to front */
  void specRunToFront();

  std::shared_ptr<b64_decoder> text_decoder_;
  rapidxml::xml_node<> * spec_paramlist_;

  /* extract the retention time for a spectrum from the rapidxml document
   * then convert from text to binary */
  inline void xmlRtntnTime(double &rt) {
    for (spec_paramlist_ = current_spec_->first_node("scanList")
                                   ->first_node("scan")
                                   ->first_node("cvParam");
         spec_paramlist_;spec_paramlist_ = spec_paramlist_->next_sibling("cvParam")) {
      if ( std::stoi(&(spec_paramlist_->first_attribute("accession")->value()[3]))==1000016) {
        if (!strcmp(spec_paramlist_->first_attribute("unitName")->value(), "minute")) {
          rt = std::stod(spec_paramlist_->first_attribute("value")->value()) * 60;
        } else {
          rt = std::stod(spec_paramlist_->first_attribute("value")->value());
        }
      }
    }
  }


  /* extract the number of points in the curren spectrum from the rapidxml
   * document and convert from text to binary */
  inline void xmlSpecPts(int &pts) {
    pts = std::stoi(current_spec_->first_attribute("defaultArrayLength")->value());
  }

  /* we need to know what the biggest and smallest mz values */
  inline void min_max_mz(){
    for (spec_paramlist_ = current_spec_->first_node("cvParam");
         spec_paramlist_; spec_paramlist_ = spec_paramlist_
                                            ->next_sibling("cvParam") ) {  
      char * acc = spec_paramlist_->first_attribute("accession")->value();
      switch ( std::stoi( &acc[3]) ) {
        case 1000528:
          smallest_mz_ = std::min(smallest_mz_,
                                  std::stod(spec_paramlist_
                                  ->first_attribute("value")
                                  ->value()));
          break;
        case 1000527:
          biggest_mz_ = std::max(biggest_mz_, 
                                 std::stod(spec_paramlist_
                                 ->first_attribute("value")
                                 ->value()));
          break;
      }
    }
  }
  //double biggest_mz_=-1.0;
  //double smallest_mz_ = 2e29;

  /* rapidxml document used for parsing mzml data */
  rapidxml::xml_document<> lcms_DOC_;

  /* rapidxml document used for writing mzml data */
  rapidxml::xml_document<> chrom_doc_;


  /* 3 xml nodes used for holding data used during each step of mzml data
   * traversal
   * current_spec is the MS spectrum of "the current point of traversal"
   * current_mz_ is the m/z "binaryDataArray" for "" 
   * current_intns is the intensity "" for ""*/
  rapidxml::xml_node<> * current_spec_;
  rapidxml::xml_node<> * current_mz_;
  rapidxml::xml_node<> * current_intns_;

  double tolerance_ = 5.0;
  double width_ = tolerance_ * 1e-6;

  double min_len_ = 2.0;
  int min_steps_ = 5;

  /* double time_len_ = ;
  double neg_time_len_ = ;*/

  double minimum_intensity_;

  int hold_fpe_;
  int hold_zlc_;

  int current_m_z_fpe_;
  int current_intensity_fpe_;

  int current_m_z_zlc_;
  int current_intensity_zlc_;

  int current_mz_el_;
  int current_intensity_el_;

  char * current_m_z_data_;
  char * current_intensity_data_;

  long unsigned int post_base64_m_z_size_;
  long unsigned int post_base64_intensity_size_;

  long unsigned int post_zlib_m_z_size_;
  long unsigned int post_zlib_intensity_size_;

  int spec_count_;

  int big_points_ = 0;

  int chrom_count_;

  /* point count of the biggest window */
  int max_window_pop_ = 0;

  std::string readFilename_;
  std::vector<char> parsedMzML_;

  std::vector<char> tmp_mz_;
  std::vector<char> tmp_intns_;
  std::vector<char> tmp_cmprs_;

  std::vector<char> hold_floats_;

  std::vector<double> mzs_;
  std::vector<double> intns_;

  std::vector<double> windows_;
  std::vector<mz_window> mz_windows_;
  std::vector<point> point_windows_;

  std::vector<int> win_to_mzwin_;
  std::vector<int> mzwin_to_win_;;

  std::vector<bool> bigger_NT_;
  std::vector<bool> smaller_NT_;

  /* addresses of points in point_windows */
  /* this list will be sorted so that 
   * point_windows_[windows_addresses_[0]] is the highest intensity point */
  std::vector<int> window_addresses_;


  std::vector<spectrum> spectra_;

  int INSTAT_;
  int MZSTAT_;

  double biggest_mz_=-1.0;
  double smallest_mz_ = 2e29;

  std::vector<std::shared_ptr<chromatogram>> chromatograms;

}; // class specToChrom

} // namespace asaristc

#else
#endif
