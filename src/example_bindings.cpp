#include <pybind11/pybind11.h>
// include numpy header for usage of array_t
#include <pybind11/numpy.h>

#include "utils.h"


namespace py = pybind11;

void bind_Detection(py::module &m) {
    m.def(
        "Detection",
        &Detection,
        R"pbdoc(
            Example function performing OLS.

            Parameters
            ----------
            arr : np.array
                input array

            Returns
            -------
            coeffs: np.ndarray
                coefficients
            std_err : np.ndarray
                standard error on the coefficients
        )pbdoc"
        // py::arg("TxAntNum"),
        // py::arg("RxAntNum"),
        // py::arg("SNRdB"),
        // py::arg("ModType"),
        // py::arg("sample")
    );
}

PYBIND11_MODULE(mimo, m) {
    bind_Detection(m);
}
