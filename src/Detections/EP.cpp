#include <cstring>
#include "EP.h"
#include "utils.h"
#include <cmath>
#include <Detection.h>


EP::EP(int iter, double delta): DetectionAlgorithmRD() {
    this -> iter = iter;
    this -> delta = delta;

    NvInv = 0;
    Cons2 = nullptr;
    Alpha = nullptr;
    Gamma = nullptr;
    AlphaInit = nullptr;
    Alpha_new = nullptr;
    Gamma_new = nullptr;
    sig = nullptr;
    h2 = nullptr;
    t = nullptr;
    prob = nullptr;
    sigma2_p = nullptr;
    mu_p = nullptr;
    HtH = nullptr;
    HtHMod = nullptr;
    HtR = nullptr;
    Sigma_q = nullptr;
    Mu_q = nullptr;
    invOneMinusSigAlpha = nullptr;
    HtRAddGamma = nullptr;
}

void EP::bind(Detection* detection) {
    DetectionAlgorithmRD::bind(detection);

    NvInv = detectionRD -> NvInv;
    Cons2 = detectionRD -> Cons2;

    Alpha = new double[TxAntNum2];
    Gamma = new double[TxAntNum2];

    AlphaInit = new double[TxAntNum2];
    for (int i = 0; i < TxAntNum2; i++) {
        AlphaInit[i] = 2.0;
    }

    Alpha_new = new double[TxAntNum2];
    Gamma_new = new double[TxAntNum2];


    sig = new double[TxAntNum2];
    h2 = new double[TxAntNum2];
    t = new double[TxAntNum2];

    prob = new double [TxAntNum2*ConSize];

    sigma2_p = new double[TxAntNum2];

    mu_p = new double[TxAntNum2];

    HtH = new double [TxAntNum2*TxAntNum2];

    HtHMod = new double [TxAntNum2*TxAntNum2];
    HtR = new double[TxAntNum2];

    Sigma_q = new double [TxAntNum2*TxAntNum2];

    Mu_q = new double[TxAntNum2];


    invOneMinusSigAlpha = new double[TxAntNum2];
    HtRAddGamma = new double[TxAntNum2];
}

EP::~EP() {
    delete[] Alpha;
    delete[] Gamma;
    delete[] AlphaInit;
    delete[] Alpha_new;
    delete[] Gamma_new;
    delete[] sig;
    delete[] h2;
    delete[] t;
    delete[] prob;
    delete[] sigma2_p;
    delete[] mu_p;
    delete[] HtH;
    delete[] HtHMod;
    delete[] HtR;
    delete[] Sigma_q;
    delete[] Mu_q;
    delete[] invOneMinusSigAlpha;
    delete[] HtRAddGamma;
}

void EP::execute(){
    std::memcpy(Alpha, AlphaInit, TxAntNum2 * sizeof(double));
    memset(Gamma, 0, TxAntNum2 * sizeof(double));
    memset(Alpha_new, 0, TxAntNum2 * sizeof(double));
    memset(Gamma_new, 0, TxAntNum2 * sizeof(double));



    MatrixTransposeMultiplySelf(H, RxAntNum2, TxAntNum2, HtH, NvInv);
 
    memcpy(HtHMod, HtH, TxAntNum2 * TxAntNum2 * sizeof(double));
    for (int i = 0; i < TxAntNum2; i++) {
        HtHMod[i * TxAntNum2 + i] += 2;
    }


    // HtR = H' * R / Nv
    MatrixTransposeMultiplyVector(H, RxSymbols, RxAntNum2, TxAntNum2, HtR, NvInv);
 

    memcpy(Sigma_q, HtHMod, TxAntNum2 * TxAntNum2 * sizeof(double));
    solveHermitianPositiveDefiniteInv(Sigma_q, TxAntNum2);

    SymMatrixMultiplyVector(Sigma_q, HtR, TxAntNum2,  Mu_q);


    for (int i = 0; i < iter; i++) {
        for (int j = 0; j < TxAntNum2; j++) {

            // sig = diag(Sigma_q)
            sig[j] = Sigma_q[j * TxAntNum2 + j];

            // 1 ./ (1 - sig * Alpha)
            invOneMinusSigAlpha[j] = 1.0 / (1.0 - sig[j] * Alpha[j]);

            // h2 = sig ./ (1 - sig * Alpha)
            h2[j] = sig[j] * invOneMinusSigAlpha[j];

            // t = (Mu_q - sig * Gamma) ./ (1 - sig * Alpha)
            t[j] = (Mu_q[j] - sig[j] * Gamma[j]) * invOneMinusSigAlpha[j];

            double probSum = 0;
            for(int k = 0; k < ConSize; k++) {
                // prob = exp( -( t - sym' ).^2./( 2 * h2 ));
                prob[j * ConSize + k] = std::exp(-(t[j] - Cons[k]) * (t[j] - Cons[k]) / (2 * h2[j]));

                probSum += prob[j * ConSize + k];
            }

            probSum = 1.0 / probSum;

            mu_p[j] = 0;
            for(int k = 0; k < ConSize; k++) {
                prob[j * ConSize + k] *= probSum;

                // mu_p = prob * sym;
                mu_p[j] += prob[j * ConSize + k] * Cons[k];
            }

            // sigma2_p=prob*(sym.^2)-mu_p.^2;
            sigma2_p[j] = -mu_p[j] * mu_p[j];

            for(int k = 0; k < ConSize; k++) {
                sigma2_p[j] += prob[j * ConSize + k] * Cons2[k];
            }

            sigma2_p[j] = sigma2_p[j] < 5e-7 ? 5e-7 : sigma2_p[j];

            // tempAlpha = 1./sigma2_p - 1./ h2;
            double tempAlpha = 1.0 / sigma2_p[j] - 1.0 / h2[j];

            if (tempAlpha > 5e-7) {

                Alpha_new[j] = tempAlpha;
                Gamma_new[j] = mu_p[j] / sigma2_p[j] - t[j] / h2[j];

            }

            Alpha[j] = Alpha_new[j] * delta + Alpha[j] * (1 - delta);
            Gamma[j] = Gamma_new[j] * delta + Gamma[j] * (1 - delta);
        }

        for(int j = 0; j < TxAntNum2; j++) {
            HtHMod[j * TxAntNum2 + j] = HtH[j * TxAntNum2 + j] + Alpha[j];
        }

        for (int j = 0; j < TxAntNum2; j++) {
            HtRAddGamma[j] = HtR[j] + Gamma[j];
        }

        memcpy(Sigma_q, HtHMod, TxAntNum2 * TxAntNum2 * sizeof(double));
        solveHermitianPositiveDefiniteInv(Sigma_q, TxAntNum2);

        SymMatrixMultiplyVector(Sigma_q, HtRAddGamma, TxAntNum2,  Mu_q);

    }

    symbolsToBits(Mu_q);

}

