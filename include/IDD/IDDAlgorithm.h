#pragma once

#include <vector>
#include <stdexcept>
#include "MIMO.h"

#include <utility>

class MIMO;

template <typename DetectionAlgType, typename CodeType>
class IDDAlgorithm
{
protected:
    std::vector<DetectionAlgType *> detectionAlgorithms;
    CodeType *code;
    int mainLoop, detectionLoop, decodingLoop;

    int blockNum;

public:
    template <typename... Args>
    IDDAlgorithm(int mainLoop, int detectionLoop, int decodingLoop, Args... args);

    virtual void execute() = 0;
    CodeType *getCode() { return code; }
};

template <typename DetectionAlgType, typename CodeType>
template <typename... Args>
IDDAlgorithm<DetectionAlgType, CodeType>::IDDAlgorithm(int mainLoop, int detectionLoop, int decodingLoop, Args... args)
{

    this->mainLoop = mainLoop;
    this->detectionLoop = detectionLoop;
    this->decodingLoop = decodingLoop;

    auto mimo = MIMO::getMIMO();
    blockNum = mimo->blockNum;

    if (blockNum == 0)
    {
        throw std::runtime_error("add code and detection first");
    }

    for (int i = 0; i < blockNum; i++)
    {
        DetectionAlgType *detectionAlg = new DetectionAlgType(std::forward<Args>(args)...);
        detectionAlgorithms.push_back(detectionAlg);
    }

    // check if mimo->code is the same type as CodeType
    code = dynamic_cast<CodeType *>(mimo->code);
    if (code == nullptr)
    {
        throw std::runtime_error("Code type mismatch");
    }

    for (int i = 0; i < mimo->blockNum; i++)
    {
        detectionAlgorithms[i]->bind(mimo->detections[i]);
    }
}
