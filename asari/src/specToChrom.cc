#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include "specToChrom.h"

namespace asaristc {

specToChrom::specToChrom() {
  current_spec_ = nullptr;

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
    if (i==0) {std::cout << npts << std::endl;}
    spectra_.push_back(spectrum(npts, rt, i, ofst));
    ofst += npts;
  }
  rts_.resize(ofst);
  intns_.resize(ofst);
}

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

void specToChrom::findChromatograms(){
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
  current_spec_ = nullptr;
}

} // namespace asaristc
