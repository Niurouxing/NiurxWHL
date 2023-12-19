#include "MMSE.h"
#include "utils.h"
#include "Mimo.h"

#include <iostream>




int main(){

    Mimo::createMimo(16,32,4,20);

    std::cout << "Cons : " << std::endl;
    printVector(Mimo::getMimo()->Cons, Mimo::getMimo()->ConSize);

    std::cout << "bitCons : " << std::endl;
    printVector(Mimo::getMimo()->bitCons, Mimo::getMimo()->ConSize * Mimo::getMimo()->bitLength);

    Algorithm * mmse = new MMSE();

    Mimo::getMimo()->reset();

    std::cout << "H : " << std::endl;
    printMatrix(Mimo::getMimo()->H, Mimo::getMimo()->RxAntNum2, Mimo::getMimo()->TxAntNum2);

    std::cout << "RxSymbols : " << std::endl;
    printVector(Mimo::getMimo()->RxSymbols, Mimo::getMimo()->RxAntNum2);

    std::cout << "TxSymbols : " << std::endl;
    printVector(Mimo::getMimo()->TxSymbols, Mimo::getMimo()->TxAntNum2);

    std::cout << "TxBits : " << std::endl;
    printMatrix(Mimo::getMimo()->TxBits, Mimo::getMimo()->TxAntNum2, Mimo::getMimo()->bitLength);

    mmse->execute();

    std::cout << "TxBitsEst : " << std::endl;
    printMatrix(Mimo::getMimo()->TxBitsEst, Mimo::getMimo()->TxAntNum2, Mimo::getMimo()->bitLength);

    for (int i = 0; i < 100; i++) {
        Mimo::getMimo()->reset();
        mmse->execute();
        mmse->check();
    }

    std::cout << "Error Bits : " << mmse->getErrorBits() << std::endl;
    std::cout << "Error Frames : " << mmse->getErrorFrames() << std::endl;

    return 0;
}