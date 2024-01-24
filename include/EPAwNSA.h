#pragma once

#include "DetectionAlgorithm.h"
#include <vector>



class EPAwNSA: public DetectionAlgorithmRD {
    private:
    double delta;
    int NSAiter;
    int iter;

    double * bhat;
 
    double * A;
 
    double * DInv;
    double * mu;
    double * mu_new;
    double * t;
    double * eta;
    double * m;
    double * momentum;
    double * grad;

    std::vector<double> alphaVec;
    std::vector<double> accuVec;

 
    public:
        EPAwNSA(double delta=0.9, int NSAiter=5,int iter=5);
        void bind(Detection* detection) override;
        void execute() override;
        void setAlphaVec(std::vector<double> alphaVec);
        void setAccuVec(std::vector<double> accuVec);
 


~EPAwNSA() override;
};