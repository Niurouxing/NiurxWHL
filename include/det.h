#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"


std::tuple<int,int> det(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){

    Detection * det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    DetectionAlgorithm * alg = new MMSE();

    alg->bind(det);

    auto start = std::chrono::high_resolution_clock::now();

    for(int i=0;i<sample;i++){
        det->generate();
        alg->execute();
        alg->check();
    }
    
    return std::make_tuple(alg->getErrorBits(), alg->getErrorFrames());
}

