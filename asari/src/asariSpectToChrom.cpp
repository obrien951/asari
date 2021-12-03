#include <pybind11/pybind11.h>
#include "specToChrom.h"

namespace py = pybind11;


PYBIND11MODULE(example, m) {
  py::class_<asaristc::specToChrom>(m, "specToChrom")
    .def(py::init<>())
    .def("", &asaristc::specToChrom::set_filename);
}
