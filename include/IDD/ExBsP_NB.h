#pragma once
#include "IDDAlgorithm.h"
#include <vector>
 

class ExBsPCD;
class NBLDPC;
class MIMO;

class ExBsP_NB : public IDDAlgorithm<ExBsPCD, NBLDPC>
{
private:
    const int *BinGF_aux;

public:
    template <typename... Args>
    ExBsP_NB(int mainLoop, int detectionLoop, int decodingLoop, Args... args)
        : IDDAlgorithm<ExBsPCD, NBLDPC>(mainLoop, detectionLoop, decodingLoop, std::forward<Args>(args)...)
    {
        // check NBLDPC GF size
        initialize();
    }

 
    void initialize();
 

    void execute() override;
};
