#include <cstring>
#include "EPAwNSA.h"
#include "utils.h"
#include <cmath>
#include <Detection.h>
#include <fstream>
#include <vector>

EPAwNSA::EPAwNSA(double delta, int NSAiter, int iter) : DetectionAlgorithmRD()
{
    this->iter = iter;
    this->NSAiter = NSAiter;
    this->delta = delta;

    bhat = nullptr;

    A = nullptr;

    DInv = nullptr;
    ps = nullptr;
    mu = nullptr;
    mu_new = nullptr;
    mu_0 = nullptr;
    t = nullptr;
    eta = nullptr;
    m = nullptr;

    alphaVec = std::vector<double>(TxAntNum2, 0.5);
    accuVec = std::vector<double>(NSAiter, 1.0);

 
}

void EPAwNSA::setAlphaVec(std::vector<double> alphaVec)
{   
    // check if the size of alphaVec is equal to TxAntNum2
    if(alphaVec.size() != TxAntNum2){
        std::cout << "The size of alphaVec is not equal to TxAntNum2" << std::endl;
        std::cout << "alphaVec.size() = " << alphaVec.size() << std::endl;
        std::cout << "TxAntNum2 = " << TxAntNum2 << std::endl;
        std::cout << "The input alphaVec is: " << std::endl;
        for(int i = 0; i < alphaVec.size(); i++){
            std::cout << alphaVec[i] << " ";
        }
        exit(1);
    }
    this->alphaVec = alphaVec;
}

void EPAwNSA::setAccuVec(std::vector<double> accuVec)
{
    // check if the size of accuVec is equal to NSAiter
    if(accuVec.size() != NSAiter){
        std::cout << "The size of accuVec is not equal to NSAiter" << std::endl;
        std::cout << "accuVec.size() = " << accuVec.size() << std::endl;
        std::cout << "NSAiter = " << NSAiter << std::endl;
        std::cout << "The input accuVec is: " << std::endl;
        for(int i = 0; i < accuVec.size(); i++){
            std::cout << accuVec[i] << " ";
        }
        exit(1);
    }
    this->accuVec = accuVec;
}

void EPAwNSA::bind(Detection *detection)
{
    DetectionAlgorithmRD::bind(detection);
    bhat = new double[TxAntNum2];

    A = new double[TxAntNum2 * TxAntNum2];

    DInv = new double[TxAntNum2];
    ps = new double[TxAntNum2 * TxAntNum2];
    mu = new double[TxAntNum2];
    mu_new = new double[TxAntNum2];
    mu_0 = new double[TxAntNum2];
    t = new double[TxAntNum2];
    eta = new double[TxAntNum2];
    m = new double[TxAntNum2];
}

EPAwNSA::~EPAwNSA()
{
    delete[] bhat;

    delete[] A;

    delete[] DInv;
    delete[] ps;
    delete[] mu;
    delete[] mu_new;
    delete[] mu_0;
    delete[] t;
    delete[] eta;
    delete[] m;
}

void EPAwNSA::execute()
{

    double Es = 2;
    // A = H' * H /Nv
    MatrixTransposeMultiplySelf(H, RxAntNum2, TxAntNum2, A, NvInv);

    // bhat = H' * RxSymbols / Nv
    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum2, TxAntNum2, bhat, NvInv);

    // W = A + diag(Alpha)
    // DInv = 1./diag(W)
    for (int i = 0; i < TxAntNum2; i++)
    {
        DInv[i] = 1.0 / (A[i * TxAntNum2 + i] + Es);
    }

    // ps = I - DInv * W
    for (int i = 0; i < TxAntNum2; i++)
    {
        for (int j = 0; j < TxAntNum2; j++)
        {
            if (i != j)
            {
                ps[i * TxAntNum2 + j] = -alphaVec[i] * DInv[i] * A[i * TxAntNum2 + j];
            }
        }
        ps[i * TxAntNum2 + i] = 1 - alphaVec[i];
    }

    // mu_0 = alpha * DInv * b
    for (int i = 0; i < TxAntNum2; i++)
    {
        mu_0[i] = accuVec[0] * alphaVec[i] * DInv[i] * bhat[i];
    }
    memcpy(mu, mu_0, sizeof(double) * TxAntNum2);

    for (int k = 1; k < NSAiter; k++)
    {
        double *input = (k % 2 == 0) ? mu_new : mu_0;
        double *output = (k % 2 == 0) ? mu_0 : mu_new;

        MatrixMultiplyVector(ps, input, TxAntNum2, TxAntNum2, output);

        cblas_daxpy(TxAntNum2, accuVec[k], output, 1, mu, 1);

    }

 

    // t_i = mu_i / (1 - DInv_i)
    for (int i = 0; i < TxAntNum2; i++)
    {
        t[i] = mu[i] / (1 - DInv[i] * Es);
    }

    // main loop
    for (int loop = 0; loop < iter; loop++)
    {
        // eta 为距离t最近的Cons
        for (int i = 0; i < TxAntNum2; i++)
        {
            double minDist = std::abs(t[i] - Cons[0]);
            eta[i] = Cons[0];
            for (int j = 1; j < ConSize; j++)
            {
                double dist = std::abs(t[i] - Cons[j]);
                if (dist < minDist)
                {
                    minDist = dist;
                    eta[i] = Cons[j];
                }
            }
        }

        // m = bhat - A * eta
        memcpy(m, bhat, sizeof(double) * TxAntNum2);
        cblas_dsymv(CblasRowMajor, CblasUpper,
                    TxAntNum2,
                    -1.0, A, TxAntNum2,
                    eta, 1,
                    1.0, m, 1);

        // t_i = m_i * W_{ii} + eta_i
        for (int i = 0; i < TxAntNum2; i++)
        {
            double t_new = m[i] / A[i * TxAntNum2 + i] + eta[i];
            t[i] = delta * t_new + (1 - delta) * t[i];
        }
    }
    symbolsToBits(eta);
}
