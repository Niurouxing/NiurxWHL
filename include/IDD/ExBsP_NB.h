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
    const int *BinGF_real;
    int NbOper;
    double offset;
 

public:
    template <typename... Args>
    ExBsP_NB(int mainLoop, int detectionLoop, int decodingLoop, int NbOper, double offset, Args... args)
        : IDDAlgorithm<ExBsPCD, NBLDPC>(mainLoop, detectionLoop, decodingLoop, std::forward<Args>(args)...)
    {
        this->NbOper = NbOper;
        this->offset = offset;
        // check NBLDPC GF size
        initialize();
    }

 
    void initialize();
 

    void execute() override;

    ~ExBsP_NB();

};
