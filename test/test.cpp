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
    int ModType = 8;
    double SNRdB = 20;
    int sample = 1000;
    Detection *det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    EPAwNSA *alg = new EPAwNSA(0.9, 40, 7);

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