#include <cstring>
#include <cmath>

#include "MMSE.h"  

#include "utils.h"

MMSE::MMSE(): Algorithm() {

    H=mimo->H;
    RxSymbols=mimo->RxSymbols;
    Cons=mimo->Cons;
    bitCons=mimo->bitCons;
    Nv=mimo->Nv;

    TxAntNum2=mimo->TxAntNum2;
    RxAntNum2=mimo->RxAntNum2;
    ConSize=mimo->ConSize;
    bitLength=mimo->bitLength;

    HtH = new double[TxAntNum2 * TxAntNum2];
    HtHInv = new double[TxAntNum2 * TxAntNum2];

    HtR = new double[TxAntNum2];

    L = new double[TxAntNum2 * TxAntNum2];
    y = new double[TxAntNum2];

    TxSymbolsEst = new double[TxAntNum2];
    TxBitsEst = mimo->TxBitsEst;

}

MMSE::~MMSE() {
    delete[] HtH;
    delete[] HtHInv;
    delete[] HtR;
    delete[] L;
    delete[] y;
    delete[] TxSymbolsEst;
}

void MMSE::execute() {

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum2, TxAntNum2, HtR);

    MatrixTransposeMultiplyMatrix(H, H, RxAntNum2, TxAntNum2, TxAntNum2, HtH);

    for (int i = 0; i < TxAntNum2; i++) {
        HtH[i * TxAntNum2 + i] += Nv;
    }

    std::memset(L, 0, TxAntNum2 * TxAntNum2 * sizeof(double));
    


    // Cholesky 分解
    for (int i = 0; i < TxAntNum2; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += L[i * TxAntNum2 + k] * L[j * TxAntNum2 + k];
            }
            L[i * TxAntNum2 + j] = (i == j) ? sqrt(HtH[i * TxAntNum2 + i] - sum) : (1.0 / L[j * TxAntNum2 + j] * (HtH[i * TxAntNum2 + j] - sum));
        }
    }

    // 求逆
    for (int i = 0; i < TxAntNum2; ++i) {
        memset(y, 0, TxAntNum2 * sizeof(double));

        // 前向替换求解 L * y = e_i
        for (int j = 0; j < TxAntNum2; ++j) {
            y[j] = (j == i) ? 1.0 : 0.0;
            for (int k = 0; k < j; ++k) {
                y[j] -= L[j * TxAntNum2 + k] * y[k];
            }
            y[j] /= L[j * TxAntNum2 + j];
        }

        // 后向替换求解 L^T * x = y
        for (int j = TxAntNum2 - 1; j >= 0; --j) {
            HtHInv[j * TxAntNum2 + i] = y[j];
            for (int k = j + 1; k < TxAntNum2; ++k) {
                HtHInv[j * TxAntNum2 + i] -= L[k * TxAntNum2 + j] * HtHInv[k * TxAntNum2 + i];
            }
            HtHInv[j * TxAntNum2 + i] /= L[j * TxAntNum2 + j];
        }
    }

    MatrixMultiplyVector(HtHInv, HtR, TxAntNum2, TxAntNum2, TxSymbolsEst);

    for(int i=0;i<TxAntNum2;i++){
        double minDistance = 100000000;
        int minIndex = 0;

        for(int j=0;j<ConSize;j++){
            double distance = 0;
            distance = abs(TxSymbolsEst[i] - Cons[j]);

            if(distance < minDistance){
                minDistance = distance;
                minIndex = j;
            }
        }

        for(int j=0;j<bitLength;j++){
            TxBitsEst[i * bitLength + j] = bitCons[minIndex * bitLength + j];
        }
    }


}

