#ifndef __SPECTOCHROM__
#define __SPECTOCHROM__

#include "rapidxml/rapidxml.hpp"
#include <vector>
#include <cstring>
#include <string>

namespace asaristc {

class spectrum {
public:
  spectrum(int &n_pts, double &RT, int &id, int &offset);
  void copy_values(int &count, double* mzs, double * intns);
  //sets the classes data pointers to the address study_start + offset_
  void pointers_to_offset(double* study_start);
protected:
  int n_pts_;
  double rt_;
  int id_;
  int offset_;
  double* mzs_;
  double* intns_;
};

// reads spectra from an mzML file then writes them to a chromatogram file
class specToChrom {
public:
  specToChrom();
  void set_filename(std::string filename);
  void print_filename();
  void reset();
  bool is_set() {return true;}
  void readSpectra();
  void findChromatograms();
  void writeChromatograms(std::string fname);
protected:
  void parse_xml();
  void initiate_spectra();
  void advance_run();
  rapidxml::xml_node<> * spec_paramlist_;

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

  inline void xmlSpecPts(int &pts) {
    spec_paramlist_ = current_spec_->first_node("binaryDataArrayList")
                                  ->first_node("binaryDataArray");
    pts = std::stoi(spec_paramlist_->first_attribute("encodedLength")->value());
  }

  rapidxml::xml_document<> lcms_DOC_;

  rapidxml::xml_document<> chrom_doc_;

  rapidxml::xml_node<> * current_spec_;

  int spec_count_;

  std::string readFilename_;
  std::vector<char> parsedMzML_;

  std::vector<char> tmp_mz_;
  std::vector<char> tmp_intns_;
  std::vector<char> tmp_cmprs_;

  std::vector<double> rts_;
  std::vector<double> intns_;

  std::vector<spectrum> spectra_;

}; // class specToChrom

} // namespace asaristc

#else
#endif
