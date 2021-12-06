#ifndef __SPECTOCHROM__
#define __SPECTOCHROM__

#include <cstring>
#include <string>

namespace asaristc {
  class specToChrom {
public:
  specToChrom();
  void set_filename(std::string filename);
  void print_filename();
protected:
  std::string filename_;
}; // class specToChrom
} // namespace asaristc

#else
#endif
