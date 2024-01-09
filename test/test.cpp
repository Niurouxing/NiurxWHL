#include "MMSE.h"
#include "utils.h"
#include "Detection.h"
#include "EP.h"
#include "EPAwNSA.h"
#include "ExBsP.h"

#include <iostream>
#include <chrono>

#include "ExBsP_NB.h"
#include "MIMO.h"
#include "NBLDPC.h"

int main()
{
    int TxAntNum = 64;
    int RxAntNum = 128;
    int ModType = 4;
    double SNRdB = 10;
    int sample = 100;
    Detection *det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    DetectionAlgorithm *alg = new EPAwNSA(0.9,0.5,13, 7);

    alg->bind(det);

    for (int i = 0; i < sample; i++)
    {
        det->generate();
        alg->execute();
        alg->check();
    }

    auto errorBits = alg->getErrorBits();
    auto errorFrames = alg->getErrorFrames();

    std::cout << "EP: " << errorBits << " " << errorFrames << std::endl;

    delete det;
    delete alg;
}