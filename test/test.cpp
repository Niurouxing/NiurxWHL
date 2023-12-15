#include "utils.h"

int main(int argc, char** argv) {
    int TxAntNum = 4;
    int RxAntNum = 4;
    double SNRdB = 10;
    int ModType = 4;
    int sample = 1000;

    // 创建两个相同形状的矩阵
    arma::mat A = arma::ones<arma::mat>(3, 3) * 9; // 3x3 矩阵，所有元素为 9
    arma::mat B = arma::ones<arma::mat>(3, 3) * 3; // 3x3 矩阵，所有元素为 3

    // 执行元素级除法，C的每个元素都是A中对应元素除以B中对应元素
    arma::mat C = A /B;

    // 输出结果
    C.print("Matrix C is:");

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
    // auto [RxIndice,RxSymbols] = MMSE(H,RxSignals,Cons,Nv);
    // std::cout << "RxIndice:" << std::endl;
    // std::cout << RxIndice << std::endl;
    // std::cout << "RxSymbols:" << std::endl;
    // std::cout << RxSymbols << std::endl;

    // EP
    auto [RxIndice,RxSymbols] = EP(H,RxSignals,Cons,Nv,4);
    std::cout << "RxIndice:" << std::endl;
    std::cout << RxIndice << std::endl;
    std::cout << "RxSymbols:" << std::endl;
    std::cout << RxSymbols << std::endl;

 

  return 0;
}