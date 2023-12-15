#include "utils.h"

int main(int argc, char** argv) {
    int TxAntNum = 4;
    int RxAntNum = 4;
    double SNRdB = 10;
    int ModType = 4;
    int sample = 1000;

    // Generate Cons
    auto [Cons, bitCons] = GenerateCons(ModType);
    std::cout << "Cons:" << std::endl;
    std::cout << Cons << std::endl;
    std::cout << "bitCons:" << std::endl;
    std::cout << bitCons << std::endl;

    // // Generate Channel Matrix
    arma::mat H = GenerateChannelMatrix(TxAntNum, RxAntNum);
    std::cout << "H:" << std::endl;
    std::cout << H << std::endl;

    // // Generate TxSignals
    auto [TxIndice,TxSignals] = GenerateTxSignals(TxAntNum, Cons);
    std::cout << "TxIndice:" << std::endl;
    std::cout << TxIndice << std::endl;
    std::cout << "TxSignals:" << std::endl;
    std::cout << TxSignals << std::endl;

    // // Generate RxSignals
    auto [RxSignals,Nv] = GenerateRxSignals(H,TxSignals,SNRdB,ModType);
    std::cout << "RxSignals:" << std::endl;
    std::cout << RxSignals << std::endl;

    // // MMSE
    auto [RxIndice,RxSymbols] = MMSE(H,RxSignals,Cons,Nv);
    std::cout << "RxIndice:" << std::endl;
    std::cout << RxIndice << std::endl;
    std::cout << "RxSymbols:" << std::endl;
    std::cout << RxSymbols << std::endl;
     
 

  return 0;
}