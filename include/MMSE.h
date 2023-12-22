#pragma once

#include "DetectionAlgorithm.h"
#include <complex>

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

class MMSECD: public DetectionAlgorithmCD {
    private:
        std::complex<double> ** HtH;
        std::complex<double> ** HtHInv;
        std::complex<double> * HtR;

        std::complex<double> * TxSymbolsEst;

        CholeskyInv * choleskyInv;
 
    public:
        MMSECD();
        void bind(Detection* detection) override;
        void execute() override;

        ~MMSECD()=default;

};