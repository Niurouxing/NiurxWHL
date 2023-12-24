#pragma once

#include <complex>

class CholeskyInv
{
    private:
        bool isComplex;
        double * lReal;
        double * yReal;
        std::complex<double> * lComplex;
        std::complex<double> * yComplex;
        int size;

    public:
        CholeskyInv(int size, bool isComplex = false);
        void execute(double * A, double * AInv);
        void execute(double ** A, double ** AInv);
        void execute(std::complex<double> ** A, std::complex<double> ** AInv);
        ~CholeskyInv();
};