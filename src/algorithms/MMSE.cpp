
#include "CholeskyInv.h"
#include "MMSE.h"  

#include "utils.h"

MMSE::MMSE(): DetectionAlgorithmRD() {
    HtH = nullptr;
    HtHInv = nullptr;
    HtR = nullptr;
    TxSymbolsEst = nullptr;
    choleskyInv = nullptr;
}

void MMSE::bind(Detection* detection) {
    DetectionAlgorithmRD::bind(detection);

    HtH = new double[TxAntNum2 * TxAntNum2];
    HtHInv = new double[TxAntNum2 * TxAntNum2];
    HtR = new double[TxAntNum2];
    TxSymbolsEst = new double[TxAntNum2];
    choleskyInv = new CholeskyInv(TxAntNum2);
}


void MMSE::execute() {

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum2, TxAntNum2, HtR);

    MatrixTransposeMultiplySelf(H, RxAntNum2, TxAntNum2, HtH);
 
    for (int i = 0; i < TxAntNum2; i++) {
        HtH[i * TxAntNum2 + i] += Nv;
    }

    choleskyInv->execute(HtH, HtHInv);
     
    MatrixMultiplyVector(HtHInv, HtR, TxAntNum2, TxAntNum2, TxSymbolsEst);

    symbolsToBits(TxSymbolsEst);
}

MMSE::~MMSE() {
    delete[] HtH;
    delete[] HtHInv;
    delete[] HtR;
    delete[] TxSymbolsEst;
    delete choleskyInv;
}

MMSECD::MMSECD(): DetectionAlgorithmCD() {
    HtH = nullptr;
    HtHInv = nullptr;
    HtR = nullptr;
    TxSymbolsEst = nullptr;
    choleskyInv = nullptr;
}

void MMSECD::bind(Detection* detection) {
    DetectionAlgorithmCD::bind(detection);

    HtH= new std::complex<double> [TxAntNum * TxAntNum];
    HtHInv= new std::complex<double> [TxAntNum * TxAntNum];

    HtR = new std::complex<double>[TxAntNum];
    TxSymbolsEst = new std::complex<double>[TxAntNum];
    choleskyInv = new CholeskyInv(TxAntNum, true);
}

void MMSECD::execute() {

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum, TxAntNum, HtR);

    MatrixTransposeMultiplyMatrix(H, H, RxAntNum, TxAntNum, TxAntNum, HtH);
 
    for (int i = 0; i < TxAntNum; i++) {
        HtH[i * TxAntNum + i] += Nv;
    }

    choleskyInv->execute(HtH, HtHInv);

    MatrixMultiplyVector(HtHInv, HtR, TxAntNum, TxAntNum, TxSymbolsEst);

    symbolsToBits(TxSymbolsEst);
}

MMSECD::~MMSECD() {
    
    delete[] HtH;
    delete[] HtHInv;
    delete[] HtR;
    delete[] TxSymbolsEst;
    delete choleskyInv;
}

