#pragma once

#include <vector>
class BaseCode;
class Detection;

// the class governing the MIMO system, including the detections and decodings
class MIMO
{
public:
    std::vector<Detection *> detections;
    BaseCode * code;

    int blockNum;

    static MIMO * getMIMO();

    void addCode(BaseCode * baseCode);

    void addDetection(bool isComplex, int TxAntNum, int RxAntNum, int ModType, double SNRdB);

    void generate();

    ~MIMO();
}; 