#pragma once

#include "Algorithm.h"


class MMSE: public Algorithm {
    private:
        double ** H;
        double * RxSymbols;
        const double * Cons;
        const int * bitCons;
        double  Nv;

        int TxAntNum2, RxAntNum2, ConSize, bitLength;

        double * HtH;
        double * HtHInv;

        double * HtR;

        double * L;
        double * y;

        double * TxSymbolsEst;
        int * TxBitsEst;

    public:
        MMSE();
        void execute() override;

        ~MMSE();

};