#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


#include "det.h"


void bind_det(py::module &m) {
    m.def("det", &det, "det");
    m.def("idd", &idd, "idd");
    
}

PYBIND11_MODULE(mimo, m) {
    bind_det(m);
}