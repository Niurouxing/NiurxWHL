#include "MIMO.h"

#include "Detection.h"
#include "baseCode.h"

#include <iostream>

MIMO *MIMO::getMIMO()
{
    static MIMO *mimo = new MIMO();
    return mimo;
}

void MIMO::addCode(BaseCode *baseCode)
{
    if(code != nullptr)
    {
        code->resetError();
        return;
    }
    code = baseCode;
}

void MIMO::addDetection(bool isComplex, int TxAntNum, int RxAntNum, int ModType, double SNRdB)
{
    if (detections.size() != 0)
    {
        return;
    }
    int bitsPerDetection = TxAntNum * ModType;
    if (code == nullptr)
    {
        throw std::runtime_error("No code is added to MIMO");
    }
    // check if the code length is a multiple of bitsPerDetection
    if (code->codeLength % bitsPerDetection != 0)
    {
        throw std::runtime_error("Code length is not a multiple of bitsPerDetection");
    }

    blockNum = (code->codeLength + bitsPerDetection - 1) / bitsPerDetection;
    detections.clear();
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

MIMO::~MIMO()
{
    for (int i = 0; i < detections.size(); i++)
    {
        delete detections[i];
    }
}