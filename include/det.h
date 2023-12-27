#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"
#include "ExBsP.h"


std::tuple<int,int> det(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){
    openblas_set_num_threads(1);
    Detection * det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    DetectionAlgorithm * alg = new EP(5,0.9);

    alg->bind(det);

    for(int i=0;i<sample;i++){
        det->generate();
        alg->execute();
        alg->check();
    }


    auto errorBits = alg->getErrorBits();
    auto errorFrames = alg->getErrorFrames();

    delete det;
    delete alg;


    return std::make_tuple(errorBits, errorFrames);
}
 