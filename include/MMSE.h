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

 
    public:
        MMSE();
        void bind(Detection* detection) override;
        void execute() override;

        ~MMSE() override;

};

class MMSECD: public DetectionAlgorithmCD {
    private:
        std::complex<double> * HtH;
        std::complex<double> * HtHInv;
        std::complex<double> * HtR;

        std::complex<double> * TxSymbolsEst;

 
    public:
        MMSECD();
        void bind(Detection* detection) override;
        void execute() override;

        ~MMSECD() override;

};