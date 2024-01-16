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
    int TxAntNum = 32;
    int RxAntNum = 64;
    int ModType = 4;
    double SNRdB = 12;
    int sample = 100000;
    Detection *det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    DetectionAlgorithm *alg = new EPAwNSA(0.9, 17, 8);

    alg->bind(det);

    for (int i = 0; i < sample; i++)
    {
        det->generate();
        alg->execute();
        alg->check();
    }

    double errorBits = static_cast<double>(alg->getErrorBits());
    double totalBits = static_cast<double>(TxAntNum * ModType * sample);
    double ber = errorBits / totalBits;

    std::cout << "BER: " << ber << std::endl;

    delete det;
    delete alg;
}