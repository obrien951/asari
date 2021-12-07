#include <pybind11/pybind11.h>
//#include <pybind11/embed.h>

#include "specToChrom.h"

namespace py = pybind11;

PYBIND11_MODULE(asariSpecToChrom, m) {
  py::class_<asaristc::specToChrom>(m, "specToChrom")
    .def(py::init<>())
    .def("set_filename", &asaristc::specToChrom::set_filename)
    .def("print_filename", &asaristc::specToChrom::print_filename)
    .def("is_set", &asaristc::specToChrom::is_set)
    .def("reset", &asaristc::specToChrom::reset)
    .def("readSpectra", &asaristc::specToChrom::readSpectra)
    .def("findChromatograms", &asaristc::specToChrom::findChromatograms)
    .def("writeChromatograms", &asaristc::specToChrom::writeChromatograms);
}
