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
    momentum = nullptr;
    grad = nullptr;
    DInv = nullptr;
    mu = nullptr;
    mu_new = nullptr;
    t = nullptr;
    eta = nullptr;
    m = nullptr;



 
}

void EPAwNSA::setAlphaVec(std::vector<double> alphaVec)
{   
    // check if the size of alphaVec is equal to TxAntNum2
    if(alphaVec.size() != NSAiter){
        std::cout << "The size of alphaVec is not equal to NSAiter" << std::endl;
        std::cout << "alphaVec.size() = " << alphaVec.size() << std::endl;
        std::cout << "NSAiter = " << NSAiter << std::endl;
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
    momentum = new double[TxAntNum2 * TxAntNum2];
    grad = new double[TxAntNum2 * TxAntNum2];

    DInv = new double[TxAntNum2];
    mu = new double[TxAntNum2];
    mu_new = new double[TxAntNum2];
    t = new double[TxAntNum2];
    eta = new double[TxAntNum2];
    m = new double[TxAntNum2];

    alphaVec = std::vector<double>(TxAntNum2, 0.5);
    accuVec = std::vector<double>(NSAiter, 1.0);
}

EPAwNSA::~EPAwNSA()
{
    delete[] bhat;

    delete[] A;
    delete[] momentum;
    delete[] grad;
    delete[] DInv;
    delete[] mu;
    delete[] mu_new;
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

    for (int i = 0; i < TxAntNum2; i++)
    {
        mu[i] =   accuVec[0] * alphaVec[0] * DInv[i] * bhat[i];
    }
    memset(momentum, 0, sizeof(double) * TxAntNum2 * TxAntNum2);

    // Nestrerov Accelerated Gradient Descent. Dynamic step size in alphaVec[k] and dynamic momentum shrinkage in accuVec[k]
    for (int n = 1; n < NSAiter; n++)
    {
        // mu_new = mu - accuVec[n] * momentum
        for(int i = 0; i < TxAntNum2; i++){
            mu_new[i] = mu[i] - accuVec[n] * DInv[i] * momentum[i];
        }

        // grad = -A * mu_new + bhat, grad:size(TxAntNum2,1)
        memcpy(grad, bhat, sizeof(double) * TxAntNum2);
        cblas_dsymv(CblasRowMajor, CblasUpper,
                    TxAntNum2,
                    -1.0, A, TxAntNum2,
                    mu_new, 1,
                    1.0, grad, 1);
        
        // momentum = accuVec[n] * momentum - alphaVec[n] * DInv .* grad
        for (int i = 0; i < TxAntNum2; i++)
        {
            momentum[i] = accuVec[n] * DInv[i] *momentum[i] - alphaVec[n] * DInv[i] * grad[i];
        }

        // mu = mu_new - momentum
        for (int i = 0; i < TxAntNum2; i++)
        {
            mu[i] = mu_new[i] - momentum[i];
        }
 
 

        
 
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
