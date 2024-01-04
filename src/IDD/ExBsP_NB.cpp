#include "ExBsP_NB.h"
#include "ExBsP.h"
#include "NBLDPC.h"
#include "MIMO.h"
#include "Detection.h"
#include "BinGF.h"


void ExBsP_NB::initialize()
{
    int gf = dynamic_cast<NBLDPC *>(code)->code->GF;

    switch (gf)
    {
    case 16:
        BinGF_aux = BinGF_16_aux;
        break;
    case 32:
        BinGF_aux = BinGF_32_aux;
        break;
    case 64:
        BinGF_aux = BinGF_64_aux;
        break;
    case 256:
        BinGF_aux = BinGF_256_aux;
        break;
    default:
        throw std::runtime_error("GF size not supported");
        break;
    }
 
}


void ExBsP_NB::execute()
{

    auto mimo = MIMO::getMIMO();

    auto intrinsic_LLR = code -> decoder ->intrinsic_LLR;
    int TxAntNum = mimo->detections[0]->TxAntNum;
    int ConSize = mimo->detections[0]->ConSize;

    // let each block do pre-processing
    for (auto detectionAlg : detectionAlgorithms)
    {
        detectionAlg->preProcess();
    }

    // main loop
    for (int mloop = 0; mloop < mainLoop; mloop++)
    {

        // detection loop
        for (auto detectionAlg : detectionAlgorithms)
        {
            detectionAlg->mainLoop(detectionLoop);
        }

        for (int block = 0; block < blockNum; block++)
        {
            for (int Nt = 0; Nt < TxAntNum; Nt++)
            {
                for (int sym = 0; sym < ConSize; sym++)
                {
                    int symGF = BinGF_aux[sym];
                    intrinsic_LLR[block * TxAntNum + Nt][sym] =  2;
                }
            }
        }
    }
}
