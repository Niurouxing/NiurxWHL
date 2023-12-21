#pragma once

#include "DetectionAlgorithm.h"

class CholeskyInv;

class MMSE: public DetectionAlgorithmRD {
    private:
        double * HtH;
        double * HtHInv;
        double * HtR;

        double * TxSymbolsEst;

        CholeskyInv * choleskyInv;
 
    public:
        MMSE();
        void bind(Detection* detection) override;
        void execute() override;

        ~MMSE()=default;

};