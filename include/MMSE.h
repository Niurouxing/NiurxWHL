#pragma once

#include "Algorithm.h"
#include "CholeskyInv.h"

class MMSE: public Algorithm {
    private:
        double * HtH;
        double * HtHInv;
        double * HtR;

        double * TxSymbolsEst;

        CholeskyInv * choleskyInv;
 
    public:
        MMSE();
        void execute() override;

        ~MMSE();

};