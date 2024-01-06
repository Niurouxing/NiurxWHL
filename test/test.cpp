#include "MMSE.h"
#include "utils.h"
#include "Detection.h"
#include "EP.h"
#include "ExBsP.h"

#include <iostream>
#include <chrono>

#include "ExBsP_NB.h"
#include "MIMO.h"
#include "NBLDPC.h"

int main()
{
    int TxAntNum = 4;
    int RxAntNum = 8;
    int ModType = 6;
    double SNRdB = 10;
    int sample = 100;
    static NBLDPC *nbldpc = new NBLDPC(96, 48, 64, 20, 20);

    auto mimo = MIMO::getMIMO();
    mimo->addCode(nbldpc);
    mimo->addDetection(true, TxAntNum, RxAntNum, ModType, SNRdB);

    static ExBsP_NB *exbsp = new ExBsP_NB(5, 1, 1, 100, 0.3);
    for (int i = 0; i < 3; i++)
    {
        std::cout << "i: " << i << std::endl;
        for (int i = 0; i < sample; i++)
        {
            mimo->generate();
            exbsp->execute();
        }

        auto errorBits = exbsp->getCode()->getErrorBits();
        auto errorFrames = exbsp->getCode()->getErrorFrames();

        std::cout << "errorBits: " << errorBits << std::endl;
        std::cout << "errorFrames: " << errorFrames << std::endl;
    }
    // delete nbldpc;
    // delete mimo;
    return 0;
}