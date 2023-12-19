#include <tuple>
#include "MMSE.h"
#include "utils.h"
#include "Mimo.h"


std::tuple<int,int> detection(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){
    Mimo::createMimo(TxAntNum, RxAntNum, ModType, SNRdB);
    MMSE mmse;
    for (int i = 0; i < sample; i++) {
        Mimo::getMimo()->reset();
        mmse.execute();
        mmse.check();
    }
    return std::make_tuple(mmse.getErrorBits(), mmse.getErrorFrames());
}