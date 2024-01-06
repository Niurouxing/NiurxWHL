#include "MIMO.h"

#include "Detection.h"
#include "baseCode.h"

MIMO *MIMO::getMIMO()
{
    static MIMO *mimo = new MIMO();
    return mimo;
}

void MIMO::addCode(BaseCode *baseCode)
{
    code = baseCode;
}

void MIMO::addDetection(bool isComplex, int TxAntNum, int RxAntNum, int ModType, double SNRdB)
{
    int bitsPerDetection = TxAntNum * ModType;
    if (code == nullptr)
    {
        throw std::runtime_error("No code is added to MIMO");
    }

    blockNum = (code->codeLength + bitsPerDetection - 1) / bitsPerDetection;
    for (int i = 0; i < blockNum; i++)
    {
        Detection *det = nullptr;
        if (isComplex)
        {
            det = new DetectionCD(TxAntNum, RxAntNum, ModType, SNRdB);
        }
        else
        {
            det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
        }
        detections.push_back(det);
    }
}

void MIMO::generate()
{
    code->encode();

    for (int i = 0; i < detections.size(); i++)
    {
        detections[i]->generate(code->codedBits + i * detections[i]->TxAntNum * detections[i]->ModType);
    }
}
