#pragma once

#include "DetectionAlgorithm.h"
#include <complex>

class CholeskyInv;


class ExBsPCD: public DetectionAlgorithmCD {
    private:
        int iter, dm;
        std::complex<double> ** HtH;
        std::complex<double> ** HtHInv;
        std::complex<double> * HtR;

        CholeskyInv * choleskyInv;


        double ** gamma;

        // alpha 3D Matrix with shape [TxAntNum][RxAntNum][ConSize]
        double *** alpha;

        // beta 3D Matrix with shape [RxAntNum][TxAntNum][ConSize]
        double *** beta;

        // Px 3D Matrix with shape [TxAntNum][RxAntNum][ConSize]
        double *** Px;

        // sIndex 3D Matrix with shape [TxAntNum][RxAntNum][dm]
        // represents the index of the dm most possible symbols, decided by MMSE
        int *** sIndex;

        std::complex<double> * sMean;
        // double * sVar;



        // intermediate variables here
        double * distList; // shape [ConSize], temp array used to store the distance between MMSEEst and constellation for each Tx antenna
        int * minkRes; // shape [dm], temp array used to store the index of the dm smallest elements in distList

        // double * alpha_ems; 
        // int * idx_ems;

        double * expAlpha;

        std::complex<double> * precomputedHCons;

    public:
        ExBsPCD(int iter, int dm);  
        void bind(Detection* detection) override;
        void execute() override;

        ~ExBsPCD()=default;

};