#ifndef __SPECTOCHROM__
#define __SPECTOCHROM__

#include "rapidxml/rapidxml.hpp"
#include <vector>
#include <cstring>
#include <string>

namespace asaristc {
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

  rapidxml::xml_document<> lcms_DOC_;

  rapidxml::xml_document<> chrom_doc_;


  int spec_count_;

  std::string readFilename_;
  std::vector<char> parsedMzML_;
}; // class specToChrom

class spectrum {
public:
protected:
};

} // namespace asaristc

#else
#endif
