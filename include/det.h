#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"


std::tuple<int,int> det(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){
    Detection::createDetection(TxAntNum, RxAntNum, ModType, SNRdB);
    Algorithm * alg = new EP(5,0.9);
    for (int i = 0; i < sample; i++) {
        Detection::getDetection()->generate();
        alg->execute();
        alg->check();
    }
    return std::make_tuple(alg->getErrorBits(), alg->getErrorFrames());
}