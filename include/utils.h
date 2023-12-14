#include <carma>
#include <armadillo>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>

arma::mat GenerateChannelMatrix(int TxAntNum, int RxAntNum);

 
std::tuple<arma::vec, arma::mat> GenerateCons(int ModType);

std::tuple<arma::uvec, arma::vec> GenerateTxSignals(int TxAntNum, arma::vec Cons);

std::tuple<arma::vec, double> GenerateRxSignals(arma::mat H, arma::vec TxSignals, double SNRdB, int ModType);

std::tuple<arma::uvec, arma::vec> MMSE(arma::mat H,arma::vec RxSignals,arma::vec Cons,double Nv);

int Detection(int TxAntNum, int RxAntNum, double SNRdB, int ModType, int sample);

void bind_Detection(py::module &m);
