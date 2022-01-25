#include <pybind11/pybind11.h>
//#include <pybind11/embed.h>

#include "specToChrom.h"

namespace py = pybind11;

PYBIND11_MODULE(asariSpecToChrom, m) {
  py::class_<asaristc::specToChrom>(m, "specToChrom")
    .def(py::init<>())
    .def("set_filename", &asaristc::specToChrom::set_filename)
    .def("set_minimum_intensity", &asaristc::specToChrom::set_minimum_intensity)
    .def("print_filename", &asaristc::specToChrom::print_filename)
    .def("is_set", &asaristc::specToChrom::is_set)
    .def("reset", &asaristc::specToChrom::reset)
    .def("readSpectra", &asaristc::specToChrom::readSpectra)
    .def("findChromatograms", &asaristc::specToChrom::findChromatograms)
    .def("writeChromatograms", &asaristc::specToChrom::writeChromatograms)
    .def("get_chrom", &asaristc::specToChrom::get_chrom)
    .def("get_nchrom", &asaristc::specToChrom::get_nchrom)
    .def("get_chrom_len", &asaristc::specToChrom::get_chrom_len);
}
