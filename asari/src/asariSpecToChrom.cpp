#include <pybind11/pybind11.h>
//#include <pybind11/embed.h>

#include "specToChrom.h"

namespace py = pybind11;

PYBIND11_MODULE(asariSpecToChrom, m) {
  py::class_<asaristc::specToChrom>(m, "specToChrom")
    .def(py::init<>())
    .def("set_filename", &asaristc::specToChrom::set_filename)
    .def("print_filename", &asaristc::specToChrom::print_filename);
}
