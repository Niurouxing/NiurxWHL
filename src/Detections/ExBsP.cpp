
#include "ExBsP.h"
#include "CholeskyInv.h"
#include "utils.h"
#include "Detection.h"
#include <cstring>

ExBsPCD::ExBsPCD(int iter, int dm) : DetectionAlgorithmCD()
{
    this->iter = iter;
    this->dm = dm;
    HtH = nullptr;
    HtHInv = nullptr;
    HtR = nullptr;
    choleskyInv = nullptr;

    gamma = nullptr;
    alpha = nullptr;
    beta = nullptr;
    Px = nullptr;
    sIndex = nullptr;
    sMean = nullptr;
    // sVar = nullptr;

    distList = nullptr;
    minkRes = nullptr;

    // alpha_ems = nullptr;
    // idx_ems = nullptr;

    expAlpha = nullptr;

    precomputedHCons = nullptr;
}

void ExBsPCD::bind(Detection *detection)
{
    DetectionAlgorithmCD::bind(detection);

    HtH = new std::complex<double> [TxAntNum * TxAntNum];


    HtHInv = new std::complex<double> [TxAntNum * TxAntNum];

    HtR = new std::complex<double>[TxAntNum];

    choleskyInv = new CholeskyInv(TxAntNum, true);

    gamma = new double [TxAntNum * ConSize];


    alpha = new double [TxAntNum * RxAntNum * ConSize];

    beta = new double [RxAntNum * TxAntNum * ConSize];


    Px = new double [TxAntNum * RxAntNum * ConSize];


    sIndex = new int [TxAntNum * RxAntNum * dm];


    sMean = new std::complex<double>[TxAntNum];
    // sVar = new double[TxAntNum];

    distList = new double[ConSize];
    minkRes = new int[dm];

    // alpha_ems = new double[ConSize];
    // idx_ems = new int[dm];

    expAlpha = new double[dm];

    precomputedHCons = new std::complex<double>[RxAntNum * TxAntNum * ConSize];
}

ExBsPCD::~ExBsPCD()
{

    delete[] HtH;
    delete[] HtHInv;
    delete[] HtR;
    delete choleskyInv;
    delete[] gamma;

    delete[] alpha;

    delete[] beta;

    delete[] Px;

    delete[] sIndex;

    delete[] sMean;
    // delete[] sVar;

    delete[] distList;
    delete[] minkRes;

    // delete[] alpha_ems;
    // delete[] idx_ems;

    delete[] expAlpha;

    delete[] precomputedHCons;
}


void ExBsPCD::execute()
{

    MatrixTransposeMultiplySelf(H, RxAntNum, TxAntNum, HtH);

    for (int i = 0; i < TxAntNum; i++)
    {
        HtH[i * TxAntNum + i] += Nv;
    }

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum, TxAntNum, HtR);

    choleskyInv->execute(HtH, HtHInv);

    // HtHInv * HtR, result stored in gamma
    // meanwhile, initialize gamma, alpha, Px, sIndex

    // usually, use i for TxAntNum, j for RxAntNum, k for ConSize

    memset(Px, 0, TxAntNum * RxAntNum * ConSize * sizeof(double));
    for (int i = 0; i < TxAntNum; i++)
    {
        std::complex<double> MMSEEst = 0;
        for (int j = 0; j < TxAntNum; j++)
        {
            MMSEEst += HtHInv[i * TxAntNum + j] * HtR[j];
        }

        double minDist = 1e10;
        int bestIndex = 0;
        for (int k = 0; k < ConSize; k++)
        {
            double dist = std::abs(ConsComplex[k] - MMSEEst);
            distList[k] = dist;
            gamma[i * ConSize + k] = -dist;

            if (dist < minDist)
            {
                minDist = dist;
                bestIndex = k;
            }
        }

        for (int j = 0; j < RxAntNum; j++)
        {
            std::copy_n(&gamma[i * ConSize], ConSize, &alpha[(i * RxAntNum + j) * ConSize]);
            
            Px[(i * RxAntNum + j) * ConSize + bestIndex] = 1;
        }

        mink(distList, ConSize, minkRes, dm);
        for (int j = 0; j < RxAntNum; j++)
        {
            std::copy_n(minkRes, dm, &sIndex[(i * RxAntNum + j) * dm]);
        }
    }

    // precompute H * ConsComplex
    for (int j = 0; j < RxAntNum; j++)
    {
        for (int i = 0; i < TxAntNum; i++)
        {
            for (int k = 0; k < ConSize; k++)
            {
                precomputedHCons[j * TxAntNum * ConSize + i * ConSize + k] = H[j * TxAntNum + i] * ConsComplex[k];
            }
        }
    }

    // iteration
    for (int L = 0; L < iter; L++)
    {
        for (int j = 0; j < RxAntNum; j++)
        {
            std::complex<double> mean_incoming_all = 0; // use the name in zwy's code
            // double var_incoming_all = 0;  // turn off variance according to wendy's code
            for (int i = 0; i < TxAntNum; i++)
            {
                std::complex<double> sGAI = 0;
                // double sGAI_var = 0;
                for (int k = 0; k < dm; k++)
                {
                    int index = sIndex[(i * RxAntNum + j) * dm + k];
                    sGAI += Px[(i * RxAntNum + j) * ConSize + index] * ConsComplex[index];
                    // sGAI_var += std::abs(Px[i][j][index] * ConsComplex[index]) * std::abs(Px[i][j][index] * ConsComplex[index]);
                }

                sMean[i] = H[j * TxAntNum + i] * sGAI;
                // sVar[i] = std::norm(H[j][i]) * (sGAI_var - std::norm(sGAI));

                mean_incoming_all += sMean[i];
                // var_incoming_all += sVar[i];
            }

            // Compute each beta message;
            for (int i = 0; i < TxAntNum; i++)
            {
                std::complex<double> sMean_in = mean_incoming_all - sMean[i];
                // double sVar_in = var_incoming_all - sVar[i] + Nv;
                std::complex<double> HS_0 = H[j * TxAntNum + i] * ConsComplex[0];

                for (int k = 0; k < ConSize; k++)
                {
                    std::complex<double> HS = precomputedHCons[j * TxAntNum * ConSize + i * ConSize + k];
                    // beta[j][i][k] =  -std::norm(RxSymbols[j] - sMean_in - HS) / sVar_in + std::norm(RxSymbols[j] - sMean_in - HS_0) / sVar_in;
                    beta[j * TxAntNum * ConSize + i * ConSize + k] = (-std::norm(RxSymbols[j] - sMean_in - HS) + std::norm(RxSymbols[j] - sMean_in - HS_0)) * NvInv * 0.5;
                }
            }
        }

        for (int i = 0; i < TxAntNum; i++)
        {
            for (int k = 0; k < ConSize; k++)
            {
                gamma[i * ConSize + k] = 0;
                for (int j = 0; j < RxAntNum; j++)
                {
                    gamma[i * ConSize + k] += beta[j * TxAntNum * ConSize + i * ConSize + k];
                }
            }
        }

        memset(Px, 0, TxAntNum * RxAntNum * ConSize * sizeof(double));
        for (int j = 0; j < RxAntNum; j++)
        {
            for (int i = 0; i < TxAntNum; i++)
            {
                for (int k = 0; k < ConSize; k++)
                {
                    // alpha[i][j][k] = gamma[i][k] - beta[j][i][k];
                    alpha[(i * RxAntNum + j) * ConSize + k] = gamma[i * ConSize + k] - beta[j * TxAntNum * ConSize + i * ConSize + k];
                }

                double expAlphaSum = 0.0;
                for (int k = 0; k < dm; k++)
                {
                    int idx = sIndex[(i * RxAntNum + j) * dm + k];
                    // expAlpha[k] = exp(alpha[i][j][idx] - alpha[i][j][sIndex[i][j][0]]);
                    expAlpha[k] = std::exp(alpha[(i * RxAntNum + j) * ConSize + idx] - alpha[(i * RxAntNum + j) * ConSize + sIndex[(i * RxAntNum + j) * dm + 0]]);
                    expAlphaSum += expAlpha[k];
                }

                // Update Px[i][j]
                for (int k = 0; k < dm; k++)
                {
                    int idx = sIndex[(i * RxAntNum + j) * dm + k];
                    Px[(i * RxAntNum + j) * ConSize + idx] = expAlpha[k] / expAlphaSum;
                }
            }
        }
    }

    // symbolsToBits
    for (int i = 0; i < TxAntNum; i++)
    {
        // find largest in gamma[i]
        double maxGamma = -1e10;
        int maxIndex = 0;
        for (int k = 0; k < ConSize; k++)
        {
            if (gamma[i * ConSize + k] > maxGamma)
            {
                maxGamma = gamma[i * ConSize + k];
                maxIndex = k;
            }
        }

        for (int j = 0; j < bitLength; j++)
        {
            TxBitsEst[i * bitLength + j] = bitConsComplex[maxIndex * bitLength + j];
        }
    }
}
