#include "MMSE.h"  

#include "utils.h"

MMSE::MMSE(): Algorithm() {
    HtH = new double[TxAntNum2 * TxAntNum2];
    HtHInv = new double[TxAntNum2 * TxAntNum2];
    HtR = new double[TxAntNum2];

    TxSymbolsEst = new double[TxAntNum2];

    choleskyInv = new CholeskyInv(TxAntNum2);

}

MMSE::~MMSE() {
    delete[] HtH;
    delete[] HtHInv;
    delete[] HtR;

    delete[] TxSymbolsEst;
    delete choleskyInv;
}

void MMSE::execute() {

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum2, TxAntNum2, HtR);

    MatrixTransposeMultiplyMatrix(H, H, RxAntNum2, TxAntNum2, TxAntNum2, HtH);
 
    for (int i = 0; i < TxAntNum2; i++) {
        HtH[i * TxAntNum2 + i] += Nv;
    }

    choleskyInv->execute(HtH, HtHInv);
     
    MatrixMultiplyVector(HtHInv, HtR, TxAntNum2, TxAntNum2, TxSymbolsEst);

    symbolsToBits(TxSymbolsEst);
}

