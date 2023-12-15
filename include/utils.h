#include <armadillo>

arma::mat GenerateChannelMatrix(int TxAntNum, int RxAntNum);

 
std::tuple<arma::vec, arma::mat> GenerateCons(int ModType);

arma::mat GenerateChannelMatrix(int TxAntNum, int RxAntNum);

std::tuple<arma::uvec, arma::vec> GenerateTxSignals(int TxAntNum, arma::vec Cons);

std::tuple<arma::vec, double> GenerateRxSignals(arma::mat H, arma::vec TxSignals, double SNRdB, int ModType);

std::tuple<arma::uvec, arma::vec> MMSE(arma::mat H,arma::vec RxSignals,arma::vec Cons,double Nv);

std::tuple<arma::uvec, arma::vec> EP(arma::mat H, arma::vec RxSignals, arma::vec Cons, double Nv, int iter);

int Detection(int TxAntNum, int RxAntNum, double SNRdB, int ModType, int sample);

 
