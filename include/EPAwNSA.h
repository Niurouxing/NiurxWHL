#pragma once

#include "DetectionAlgorithm.h"
#include <vector>



class EPAwNSA: public DetectionAlgorithmRD {
    private:
    double delta;
    double alpha;
    int NSAiter;
    int iter;

    double * A;
    double * bhat;
    double * Alpha;
    double * Gamma;
    double * W;
    double * b;
    double * DInv;
    double * ps;
    double * mu;
    double * mu_new;
    double * Dinvb;
    double * t;
    double * eta;
    double * m;

 
    public:
        EPAwNSA(double delta=0.9, double alpha=0.5, int NSAiter=5,int iter=5);
        void bind(Detection* detection) override;
        void execute() override;


~EPAwNSA() override;
};