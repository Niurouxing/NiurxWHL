#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"


std::tuple<int,int> det(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){

    DetectionAlgorithm * alg = new EP(5,0.9);

    return std::make_tuple(alg->getErrorBits(), alg->getErrorFrames());
}