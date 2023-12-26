#pragma once

#include "DetectionAlgorithm.h"
#include "CholeskyInv.h"

class EP: public DetectionAlgorithmRD {
    private:
        int iter;
        double delta;
        double NvInv;

        const double * Cons2;


        double * Alpha;
        double * Gamma;

        double * AlphaInit;

        double * Alpha_new;
        double * Gamma_new;



        double * sig;
        double * h2;
        double * t;

        double * prob;
        double * sigma2_p;

        double * mu_p;

        double * HtH;
        double * HtHMod;

        double * HtR;

        double * Sigma_q;

        double * Mu_q;


 
        CholeskyInv * choleskyInv;

        // 中间变量

        double * invOneMinusSigAlpha;
        double * HtRAddGamma;

 
    public:
        EP(int iter, double delta);
        void bind(Detection* detection) override;
        void execute() override;

        ~EP() override;
};