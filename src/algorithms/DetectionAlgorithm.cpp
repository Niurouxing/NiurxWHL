
#include "utils.h"
#include "DetectionAlgorithm.h"
#include "Detection.h"

DetectionAlgorithm::DetectionAlgorithm()
{
    detection = nullptr;
    errorBits = 0;
    errorFrames = 0;
}

void DetectionAlgorithmRD::bind(Detection *detection)
{
    //  transfer the pointer of detection to the pointer of DetectionRD, use dynamic_cast
    DetectionRD *rd = dynamic_cast<DetectionRD *>(detection);
    if (rd == nullptr)
    {
        throw std::runtime_error("only DetectionRD can be binded to DetectionAlgorithmRD");
    }

    this->detection = detection;
    this->detectionRD = rd;


    Nv = rd->Nv;
    TxAntNum = rd->TxAntNum;
    RxAntNum = rd->RxAntNum;
    TxAntNum2 = rd->TxAntNum2;
    RxAntNum2 = rd->RxAntNum2;
    ConSize = rd->ConSize;
    bitLength = rd->bitLength;

    H = rd->H;
    RxSymbols = rd->RxSymbols;
    Cons = rd->Cons;
    bitCons = rd->bitCons;

    TxBitsEst = new int[TxAntNum2 * bitLength];
}

DetectionAlgorithmRD::DetectionAlgorithmRD() : DetectionAlgorithm()
{
    H = nullptr;
    RxSymbols = nullptr;
    Cons = nullptr;
    bitCons = nullptr;
}

void DetectionAlgorithmRD::check()
{
    int currentErrorBits = 0;
    for (int i = 0; i < detection->TxAntNum2 * detection->bitLength; i++)
    {
        if (detection->TxBits[i] != TxBitsEst[i])
        {
            currentErrorBits++;
        }
    }
    if (currentErrorBits > 0)
    {
        errorFrames++;
        errorBits += currentErrorBits;
    }
}

void DetectionAlgorithmRD::symbolsToBits(double *TxSymbolsEst)
{
    for (int i = 0; i < TxAntNum2; i++)
    {
        double minDistance = 100000000;
        int minIndex = 0;

        for (int j = 0; j < ConSize; j++)
        {
            double distance = 0;
            distance = std::abs(TxSymbolsEst[i] - Cons[j]);

            if (distance < minDistance)
            {
                minDistance = distance;
                minIndex = j;
            }
        }

        for (int j = 0; j < bitLength; j++)
        {
            TxBitsEst[i * bitLength + j] = bitCons[minIndex * bitLength + j];
        }
    }
}
