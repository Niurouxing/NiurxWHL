#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Mimo.h"


std::tuple<int,int> detection(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){
    Mimo::createMimo(TxAntNum, RxAntNum, ModType, SNRdB);
    Algorithm * alg = new EP(5,0.9);
    for (int i = 0; i < sample; i++) {
        Mimo::getMimo()->reset();
        alg->execute();
        alg->check();
    }
    return std::make_tuple(alg->getErrorBits(), alg->getErrorFrames());
}