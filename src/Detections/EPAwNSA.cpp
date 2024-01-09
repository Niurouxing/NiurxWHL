#include <cstring>
#include "EPAwNSA.h"
#include "utils.h"
#include <cmath>
#include <Detection.h>

EPAwNSA::EPAwNSA(double delta, double alpha, int NSAiter, int iter) : DetectionAlgorithmRD()
{
    this->iter = iter;
    this->NSAiter = NSAiter;
    this->delta = delta;
    this->alpha = alpha;

    bhat = nullptr;
    Alpha = nullptr;
 
    W = nullptr;
 
    DInv = nullptr;
    ps = nullptr;
    mu = nullptr;
    mu_new = nullptr;
    Dinvb = nullptr;
    t = nullptr;
    eta = nullptr;
    m = nullptr;
}

void EPAwNSA::bind(Detection *detection)
{
    DetectionAlgorithmRD::bind(detection);
    bhat = new double[TxAntNum2];
    Alpha = new double[TxAntNum2];
 
    W = new double[TxAntNum2 * TxAntNum2];
 
    DInv = new double[TxAntNum2];
    ps = new double[TxAntNum2 * TxAntNum2];
    mu = new double[TxAntNum2];
    mu_new = new double[TxAntNum2];
    Dinvb = new double[TxAntNum2];
    t = new double[TxAntNum2];
    eta = new double[TxAntNum2];
    m = new double[TxAntNum2];
}

EPAwNSA::~EPAwNSA()
{
    delete[] bhat;
    delete[] Alpha;
 
    delete[] W;
 
    delete[] DInv;
    delete[] ps;
    delete[] mu;
    delete[] mu_new;
    delete[] Dinvb;
    delete[] t;
    delete[] eta;
    delete[] m;
}

void EPAwNSA::execute()
{

    // A = H' * H /Nv
    MatrixTransposeMultiplySelf(H, RxAntNum2, TxAntNum2, W, NvInv);

    // bhat = H' * RxSymbols / Nv
    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum2, TxAntNum2, bhat, NvInv);

    for (int i = 0; i < TxAntNum2; i++)
    {
        Alpha[i] = 2;
    }
 

    // W = A + diag(Alpha)
    for (int i = 0; i < TxAntNum2; i++)
    {
        W[i * TxAntNum2 + i] += 2;
    }

 

    // DInv = 1./diag(W)
    for (int i = 0; i < TxAntNum2; i++)
    {
        DInv[i] = 1.0 / W[i * TxAntNum2 + i];
    }

    // ps = I - DInv * W
    for (int i = 0; i < TxAntNum2; i++)
    {
        for (int j = 0; j < TxAntNum2; j++)
        {
            ps[i * TxAntNum2 + j] = - alpha * DInv[i] * W[i * TxAntNum2 + j];
        }
        ps[i * TxAntNum2 + i] += 1;
    }

    // mu = DInv * b
    for (int i = 0; i < TxAntNum2; i++)
    {
        Dinvb[i] = alpha * DInv[i] * bhat[i];
    }
    memcpy(mu, Dinvb, sizeof(double) * TxAntNum2);

    for (int k = 0; k < NSAiter; k++)
    {
        
        SymMatrixMultiplyVector(ps, mu, TxAntNum2, mu_new);
        for (int i = 0; i < TxAntNum2; i++)
        {
            mu_new[i] += Dinvb[i];
        }
        memcpy(mu, mu_new, sizeof(double) * TxAntNum2);
    }

    // t_i = mu_i / (1 - DInv_i)
    for (int i = 0; i < TxAntNum2; i++)
    {
        t[i] = mu[i] / (1 - DInv[i]);
    }

    // main loop
    for (int loop = 0; loop < iter; loop++)
    {
        // eta 为距离t最近的Cons
        for (int i = 0; i < TxAntNum2; i++)
        {
            eta[i] = Cons[0];
            double minDis = std::abs(t[i] - Cons[0]);
            for (int j = 1; j < ConSize; j++)
            {
                double dis = std::abs(t[i] - Cons[j]);
                if (dis < minDis)
                {
                    minDis = dis;
                    eta[i] = Cons[j];
                }
            }
        }

        // m_i = bhat_i - \Sum_{j=1}^{TxAntNum2} A_{ij} * eta_j
        memcpy(m, bhat, sizeof(double) * TxAntNum2);
        for (int i = 0; i < TxAntNum2; i++)
        {
            for (int j = 0; j < TxAntNum2; j++)
            {
                double A = W[i * TxAntNum2 + j];
                if (i == j) A-=2;
                m[i] -= A * eta[j];
            }
        }

        // t_i = m_i * W_{ii} + eta_i
        for (int i = 0; i < TxAntNum2; i++)
        {
            double t_new =  m[i] / W[i * TxAntNum2 + i] + eta[i];
            t[i] = delta * t_new + (1 - delta) * t[i];
        }
    }
    symbolsToBits(eta);
}
