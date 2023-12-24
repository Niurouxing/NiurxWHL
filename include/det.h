#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"
#include "ExBsP.h"


std::tuple<int,int> det(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){
    auto det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    auto alg = new EP(5, 0.9);
    alg->bind(det);

    for (int i = 0; i < sample; ++i) {
        det->generate();
        alg->execute();
        alg->check();
    }


    return std::make_tuple(alg->getErrorBits(), alg->getErrorFrames());
}
 