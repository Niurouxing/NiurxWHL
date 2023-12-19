
#include "utils.h"
#include "Algorithm.h"


Algorithm::Algorithm(){
    mimo = Mimo::getMimo();
    errorBits = 0;
    errorFrames = 0;

    H=mimo->H;
    RxSymbols=mimo->RxSymbols;
    Cons=mimo->Cons;
    bitCons=mimo->bitCons;
    Nv=mimo->Nv;

    TxAntNum2=mimo->TxAntNum2;
    RxAntNum2=mimo->RxAntNum2;
    ConSize=mimo->ConSize;
    bitLength=mimo->bitLength;

    TxBitsEst = new int[TxAntNum2 * bitLength];
}

Algorithm::~Algorithm(){
    delete[] TxBitsEst;
}


void Algorithm::check(){
    int currentErrorBits = 0;
    for (int i = 0; i < mimo->TxAntNum2 * mimo->bitLength; i++) {
        if (mimo->TxBits[i] != TxBitsEst[i]) {
            currentErrorBits++;
        }
    }
    if (currentErrorBits > 0) {
        errorFrames++;
        errorBits += currentErrorBits;
    }
}

void Algorithm::symbolsToBits(double * TxSymbolsEst){
    for(int i=0;i<TxAntNum2;i++){
        double minDistance = 100000000;
        int minIndex = 0;

        for(int j=0;j<ConSize;j++){
            double distance = 0;
            distance = std::abs(TxSymbolsEst[i] - Cons[j]);

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
