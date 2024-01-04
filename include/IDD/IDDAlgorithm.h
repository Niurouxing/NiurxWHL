#pragma once 

#include <vector>

#include "MIMO.h"

#include <utility>

class MIMO;

template<typename DetectionAlgType, typename CodeType>
class IDDAlgorithm{
protected:
    std::vector<DetectionAlgType *> detectionAlgorithms;
    CodeType * code;
public:
    template<typename... Args>
    IDDAlgorithm(Args... args);
    void bind(MIMO * mimo);
};  

template < typename DetectionAlgType, typename CodeType>
void IDDAlgorithm< DetectionAlgType, CodeType>::bind(MIMO *mimo)
{
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

template<typename DetectionAlgType, typename CodeType>
template<typename... Args>
IDDAlgorithm<DetectionAlgType, CodeType>::IDDAlgorithm(Args... args) {

    int blockNum = MIMO::getMIMO()->blockNum;
    for(int i=0;i<blockNum;i++){
        DetectionAlgType * detectionAlg = new DetectionAlgType(std::forward<Args>(args)...);
        detectionAlgorithms.push_back(detectionAlg);
    }
}