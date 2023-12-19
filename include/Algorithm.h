#pragma once
#include "Mimo.h"


// class Algorithm, a virtual base class for all Mimo detection algorithms
class Algorithm {
    protected:
        Mimo * mimo;
        int errorBits;
        int errorFrames;

        double ** H;
        double * RxSymbols;
        const double * Cons;
        const int * bitCons;
        double  Nv;

        int TxAntNum2, RxAntNum2, ConSize, bitLength;

        int * TxBitsEst;
    
    public:
        Algorithm();
        virtual void execute() = 0;
        virtual ~Algorithm()=default;
        void check();
        void symbolsToBits(double * TxSymbolsEst);

        int * getTxBitsEst() const { return TxBitsEst; }
        int getErrorBits() const { return errorBits; }
        int getErrorFrames() const { return errorFrames; }
};
        


