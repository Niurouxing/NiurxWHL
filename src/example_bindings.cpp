#include <pybind11/pybind11.h>
// include numpy header for usage of array_t
#include <pybind11/numpy.h>

#include "utils.h"


namespace py = pybind11;

PYBIND11_MODULE(mimo, m) {
    bind_Detection(m);
}
