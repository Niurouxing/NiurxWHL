
#include "utils.h"
#include "Algorithm.h"


Algorithm::Algorithm(){
    mimo = Mimo::getMimo();
    errorBits = 0;
    errorFrames = 0;
}


void Algorithm::check(){
    int currentErrorBits = 0;
    for (int i = 0; i < mimo->TxAntNum2 * mimo->bitLength; i++) {
        if (mimo->TxBits[i] != mimo->TxBitsEst[i]) {
            currentErrorBits++;
        }
    }
    if (currentErrorBits > 0) {
        errorFrames++;
        errorBits += currentErrorBits;
    }
}
