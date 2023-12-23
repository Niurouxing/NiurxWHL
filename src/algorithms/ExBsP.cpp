
#include "ExBsP.h"
#include "CholeskyInv.h"
#include "utils.h"
#include "Detection.h"

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
    sVar = nullptr;

    distList = nullptr;
    minkRes = nullptr;

    alpha_ems = nullptr;
    idx_ems = nullptr;

    expAlpha = nullptr;
}

void ExBsPCD::bind(Detection *detection)
{
    DetectionAlgorithmCD::bind(detection);

    HtH = new std::complex<double> *[TxAntNum];
    for (int i = 0; i < TxAntNum; i++)
    {
        HtH[i] = new std::complex<double>[TxAntNum];
    }

    HtHInv = new std::complex<double> *[TxAntNum];
    for (int i = 0; i < TxAntNum; i++)
    {
        HtHInv[i] = new std::complex<double>[TxAntNum];
    }

    HtR = new std::complex<double>[TxAntNum];

    choleskyInv = new CholeskyInv(TxAntNum,true);

    gamma = new double *[TxAntNum];
    for (int i = 0; i < TxAntNum; i++)
    {
        gamma[i] = new double[ConSize];
    }

    alpha = new double **[TxAntNum];
    for (int i = 0; i < TxAntNum; i++)
    {
        alpha[i] = new double *[RxAntNum];
        for (int j = 0; j < RxAntNum; j++)
        {
            alpha[i][j] = new double[ConSize];
        }
    }

    beta = new double **[RxAntNum];
    for (int i = 0; i < RxAntNum; i++)
    {
        beta[i] = new double *[TxAntNum];
        for (int j = 0; j < TxAntNum; j++)
        {
            beta[i][j] = new double[ConSize];
        }
    }

    Px = new double **[TxAntNum];
    for (int i = 0; i < TxAntNum; i++)
    {
        Px[i] = new double *[RxAntNum];
        for (int j = 0; j < RxAntNum; j++)
        {
            Px[i][j] = new double[ConSize];
        }
    }

    sIndex = new int **[TxAntNum];
    for (int i = 0; i < TxAntNum; i++)
    {
        sIndex[i] = new int *[RxAntNum];
        for (int j = 0; j < RxAntNum; j++)
        {
            sIndex[i][j] = new int[dm];
        }
    }

    sMean = new std::complex<double>[TxAntNum];
    sVar = new double[TxAntNum];

    distList = new double[ConSize];
    minkRes = new int[dm];

    alpha_ems = new double[ConSize];
    idx_ems = new int[dm];

    expAlpha = new double[dm];
}
void ExBsPCD::execute()
{

    MatrixTransposeMultiplyMatrix(H, H, RxAntNum, TxAntNum, TxAntNum, HtH);

    for (int i = 0; i < TxAntNum; i++)
    {
        HtH[i][i] += Nv;
    }

    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum, TxAntNum, HtR);

    choleskyInv->execute(HtH, HtHInv);

    // HtHInv * HtR, result stored in gamma
    // meanwhile, initialize gamma, alpha, Px, sIndex

    // usually, use i for TxAntNum, j for RxAntNum, k for ConSize
    for (int i = 0; i < TxAntNum; i++)
    {
        std::complex<double> MMSEEst = 0;
        for (int j = 0; j < TxAntNum; j++)
        { // exception, j iterate TxAntNum here
            MMSEEst += HtHInv[i][j] * HtR[j];
        }

        double minDist = 1e10;
        int bestIndex = 0;
        // gamma[i][k] = - abs ( Cons[k] - MMSEres )
        for (int k = 0; k < ConSize; k++)
        {
            double dist = std::abs(ConsComplex[k] - MMSEEst);
            distList[k] = dist;
            gamma[i][k] = -dist;

            // initialize alpha
            for (int j = 0; j < RxAntNum; j++)
            {
                alpha[i][j][k] = -dist;
            }

            if (dist < minDist)
            {
                minDist = dist;
                bestIndex = k;
            }
        }

        // initialize Px, the probability of the best symbol is 1, others are 0
        // the 0 is gruanteed by memset
        for (int j = 0; j < RxAntNum; j++)
        {
            memset(Px[i][j], 0, ConSize * sizeof(double));
            Px[i][j][bestIndex] = 1;
        }

        // initialize sIndex with the dm most least distance symbols
        mink(distList, ConSize, minkRes, dm);
        for (int j = 0; j < RxAntNum; j++)
        {
            for (int k = 0; k < dm; k++)
            {
                sIndex[i][j][k] = minkRes[k];
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
                    int index = sIndex[i][j][k];
                    sGAI += Px[i][j][index] * ConsComplex[index];
                    // sGAI_var += std::abs(Px[i][j][index] * ConsComplex[index]) * std::abs(Px[i][j][index] * ConsComplex[index]);
                }

                sMean[i] = H[j][i] * sGAI;
                // sVar[i] = std::norm(H[j][i]) * (sGAI_var - std::norm(sGAI));

                mean_incoming_all += sMean[i];
                // var_incoming_all += sVar[i];
            }

            // Compute each beta message;
            for (int i = 0; i < TxAntNum; i++)
            {
                std::complex<double> sMean_in = mean_incoming_all - sMean[i];
                // double sVar_in = var_incoming_all - sVar[i] + Nv;
                std::complex<double> HS_0 = H[j][i] * ConsComplex[0];

                for (int k = 0; k < ConSize; k++)
                {
                    std::complex<double> HS = H[j][i] * ConsComplex[k];
                    // beta[j][i][k] =  -std::norm(RxSymbols[j] - sMean_in - HS) / sVar_in + std::norm(RxSymbols[j] - sMean_in - HS_0) / sVar_in;
                    beta[j][i][k] =  (-std::norm(RxSymbols[j] - sMean_in - HS) + std::norm(RxSymbols[j] - sMean_in - HS_0)) / (2*Nv);

         
                }
            }
        }

        for (int i = 0; i < TxAntNum; i++)
        {
            for (int k = 0; k < ConSize; k++)
            {
                gamma[i][k] = 0;
                for (int j = 0; j < RxAntNum; j++)
                {
                    gamma[i][k] += beta[j][i][k];
                }
            }
        }

        // printMatrix(gamma, TxAntNum, ConSize, "gamma");

        for (int j = 0; j < RxAntNum; j++)
        {
            for (int i = 0; i < TxAntNum; i++)
            {
                for (int k = 0; k < ConSize; k++)
                {
                    alpha[i][j][k] = gamma[i][k] - beta[j][i][k];

                    alpha_ems[k] = alpha[i][j][k];
                }

                maxk(alpha_ems, ConSize, idx_ems, dm);
                // printVector(alpha_ems, ConSize, "alpha_ems");
                // printVector(idx_ems, dm, "idx_ems");

                double expAlphaSum = 0;
                for (int k = 0; k < dm; k++)
                {
                    expAlpha[k] = exp(alpha_ems[idx_ems[k]] - alpha_ems[idx_ems[0]]);
                    expAlphaSum += expAlpha[k];
                }

                memset(Px[i][j], 0, ConSize * sizeof(double));
                for (int k = 0; k < dm; k++)
                {
                    Px[i][j][idx_ems[k]] = expAlpha[k] / expAlphaSum;
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
            if (gamma[i][k] > maxGamma)
            {
                maxGamma = gamma[i][k];
                maxIndex = k;
            }
        }

        for (int j = 0; j < bitLength; j++)
        {
            TxBitsEst[i * bitLength + j] = bitConsComplex[maxIndex * bitLength + j];
        }
    }

}
