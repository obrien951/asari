#include <cstring>
#include <string>
#include <iostream>

#include "specToChrom.h"

namespace asaristc {

specToChrom::specToChrom() {}

void specToChrom::set_filename(std::string filename) {
  filename_ = filename;
}

void specToChrom::print_filename () {
  std::cout << filename_.c_str() << std::endl;
}

} // namespace asaristc
