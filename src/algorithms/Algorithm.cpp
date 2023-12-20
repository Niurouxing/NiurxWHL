
#include "utils.h"
#include "Algorithm.h"


Algorithm::Algorithm(){
    detection = Detection::getDetection();
    errorBits = 0;
    errorFrames = 0;

    H=detection->H;
    RxSymbols=detection->RxSymbols;
    Cons=detection->Cons;
    bitCons=detection->bitCons;
    Nv=detection->Nv;

    TxAntNum2=detection->TxAntNum2;
    RxAntNum2=detection->RxAntNum2;
    ConSize=detection->ConSize;
    bitLength=detection->bitLength;

    TxBitsEst = new int[TxAntNum2 * bitLength];
}


void Algorithm::check(){
    int currentErrorBits = 0;
    for (int i = 0; i < detection->TxAntNum2 * detection->bitLength; i++) {
        if (detection->TxBits[i] != TxBitsEst[i]) {
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
