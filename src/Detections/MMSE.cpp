
#include "CholeskyInv.h"
#include "MMSE.h"  

#include "utils.h"

MMSE::MMSE(): DetectionAlgorithmRD() {
    HtH = nullptr;
    HtHInv = nullptr;
    HtR = nullptr;
    TxSymbolsEst = nullptr;
}

void MMSE::bind(Detection* detection) {
    DetectionAlgorithmRD::bind(detection);

    HtH = new double[TxAntNum2 * TxAntNum2];
    HtHInv = new double[TxAntNum2 * TxAntNum2];
    HtR = new double[TxAntNum2];
    TxSymbolsEst = new double[TxAntNum2];
}


void MMSE::execute() {

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum2, TxAntNum2, HtR);

    MatrixTransposeMultiplySelf(H, RxAntNum2, TxAntNum2, HtH);
 
    for (int i = 0; i < TxAntNum2; i++) {
        HtH[i * TxAntNum2 + i] += Nv;
    }

    solveHermitianPositiveDefiniteSystem(HtH, HtR, TxAntNum2);

    symbolsToBits(HtR);
}

MMSE::~MMSE() {
    delete[] HtH;
    delete[] HtHInv;
    delete[] HtR;
    delete[] TxSymbolsEst;
}

MMSECD::MMSECD(): DetectionAlgorithmCD() {
    HtH = nullptr;
    HtHInv = nullptr;
    HtR = nullptr;
    TxSymbolsEst = nullptr;
}

void MMSECD::bind(Detection* detection) {
    DetectionAlgorithmCD::bind(detection);

    HtH= new std::complex<double> [TxAntNum * TxAntNum];
    HtHInv= new std::complex<double> [TxAntNum * TxAntNum];

    HtR = new std::complex<double>[TxAntNum];
    TxSymbolsEst = new std::complex<double>[TxAntNum];
}

void MMSECD::execute() {

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum, TxAntNum, HtR);

    MatrixTransposeMultiplySelf(H, RxAntNum, TxAntNum, HtH);
 
    for (int i = 0; i < TxAntNum; i++) {
        HtH[i * TxAntNum + i] += Nv;
    }

    solveHermitianPositiveDefiniteSystem(HtH, HtR, TxAntNum);


    symbolsToBits(HtR);
}

MMSECD::~MMSECD() {
    
    delete[] HtH;
    delete[] HtHInv;
    delete[] HtR;
    delete[] TxSymbolsEst;
}

