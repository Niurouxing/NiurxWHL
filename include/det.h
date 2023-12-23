#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"
#include "ExBsP.h"



std::tuple<int,int> det(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample){

    Detection * det = new DetectionCD(TxAntNum, RxAntNum, ModType, SNRdB);
    DetectionAlgorithm * alg = new ExBsPCD(3, 8);


    alg->bind(det);

 

    for(int i=0;i<sample;i++){
        det->generate();
        alg->execute();
        alg->check();
    }
    
    return std::make_tuple(alg->getErrorBits(), alg->getErrorFrames());
}

