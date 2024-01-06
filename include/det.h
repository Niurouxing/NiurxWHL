#pragma once
#include <tuple>
#include "MMSE.h"
#include "EP.h"
#include "utils.h"
#include "Detection.h"
#include "ExBsP.h"
#include "ExBsP_NB.h"
#include "MIMO.h"
#include "NBLDPC.h"
#include <iostream>

std::tuple<int, int> det(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample)
{
    openblas_set_num_threads(1);
    Detection *det = new DetectionCD(TxAntNum, RxAntNum, ModType, SNRdB);
    DetectionAlgorithm *alg = new ExBsPCD(2, 10);

    alg->bind(det);

    for (int i = 0; i < sample; i++)
    {
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

std::tuple<int, int> idd(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample)
{

    openblas_set_num_threads(1);
    static NBLDPC *nbldpc = new NBLDPC(96, 48, 64, 20, 20);

    auto mimo = MIMO::getMIMO();
    mimo->addCode(nbldpc);
    mimo->addDetection(true, TxAntNum, RxAntNum, ModType, SNRdB);

    static ExBsP_NB *exbsp = new ExBsP_NB(5, 2, 2, 100, 0.3);

    for (int i = 0; i < sample; i++)
    {

        mimo->generate();

        exbsp->execute();
    }

    auto errorBits = exbsp->getCode()->getErrorBits();
    auto errorFrames = exbsp->getCode()->getErrorFrames();

    return std::make_tuple(errorBits, errorFrames);
}