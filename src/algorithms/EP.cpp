#include <cstring>
#include "EP.h"

#include "utils.h"


EP::EP(int iter, double delta): Algorithm() {
    this -> iter = iter;
    this -> delta = delta;

    choleskyInv = new CholeskyInv(TxAntNum2);

    NvInv = mimo -> NvInv;
    Cons2 = mimo -> Cons2;

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

    prob = new double*[TxAntNum2];
    for (int i = 0; i < TxAntNum2; i++) {
        prob[i] = new double[2];
    }

    sigma2_p = new double[TxAntNum2];

    mu_p = new double[TxAntNum2];

    HtH = new double*[TxAntNum2];
    for (int i = 0; i < TxAntNum2; i++) {
        HtH[i] = new double[TxAntNum2];
    }

    HtHMod = new double*[TxAntNum2];
    for (int i = 0; i < TxAntNum2; i++) {
        HtHMod[i] = new double[TxAntNum2];
    }

    HtR = new double[TxAntNum2];

    Sigma_q = new double*[TxAntNum2];
    for (int i = 0; i < TxAntNum2; i++) {
        Sigma_q[i] = new double[TxAntNum2];
    }

    Mu_q = new double[TxAntNum2];

    choleskyInv = new CholeskyInv(TxAntNum2);

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

    for (int i = 0; i < TxAntNum2; i++) {
        delete[] prob[i];
    }
    delete[] prob;

    delete[] sigma2_p;

    delete[] mu_p;

    for (int i = 0; i < TxAntNum2; i++) {
        delete[] HtH[i];
    }
    delete[] HtH;

    for (int i = 0; i < TxAntNum2; i++) {
        delete[] HtHMod[i];
    }
    delete[] HtHMod;

    delete[] HtR;

    for (int i = 0; i < TxAntNum2; i++) {
        delete[] Sigma_q[i];
    }
    delete[] Sigma_q;

    delete[] Mu_q;

    delete choleskyInv;

    delete[] HtRAddGamma;
}

void EP::execute(){
    memcpy(Alpha, AlphaInit, TxAntNum2 * sizeof(double));
    memset(Gamma, 0, TxAntNum2 * sizeof(double));
    memset(Alpha_new, 0, TxAntNum2 * sizeof(double));
    memset(Gamma_new, 0, TxAntNum2 * sizeof(double));

    // HtH = H' * H / Nv
    for (int i = 0; i < TxAntNum2; i++) {
        for (int j = 0; j < TxAntNum2; j++) {
            double sum = 0;
            for (int k = 0; k < RxAntNum2; k++) {
                sum += H[k][i] * H[k][j];
            }
            sum *= NvInv;
            HtH[i][j] = sum;
            HtHMod[i][j] = i == j ? sum + 2 : sum;
        }
    }

    // HtR = H' * R / Nv
    for (int i = 0; i < TxAntNum2; i++) {
        double sum = 0;
        for (int j = 0; j < RxAntNum2; j++) {
            sum += H[j][i] * RxSymbols[j];
        }
        HtR[i] = sum * NvInv;
    }


    // Sigma_q = inv(HtHMod)
    choleskyInv -> execute(HtHMod, Sigma_q);

    // Mu_q = Sigma_q * HtR
    MatrixMultiplyVector(Sigma_q, HtR, TxAntNum2, TxAntNum2, Mu_q);

    for (int i = 0; i < iter; i++) {
        for (int j = 0; j < TxAntNum2; j++) {

            // sig = diag(Sigma_q)
            sig[j] = Sigma_q[j][j];

            // h2 = sig ./ (1 - sig * Alpha)
            h2[j] = 1.0 / sig[j] - Alpha[j];

            // t = (Mu_q - sig * Gamma) ./ (1 - sig * Alpha)
            t[j] = (Mu_q[j] - sig[j] * Gamma[j]) / (1.0 - sig[j] * Alpha[j]);

            double probSum = 0;
            for(int k = 0; k < ConSize; k++) {
                // prob = exp( -( t - sym' ).^2./( 2 * h2 ));
                prob[j][k] = std::exp(-(t[j] - Cons[k]) * (t[j] - Cons[k]) * h2[j] * 0.5);

                probSum += prob[j][k];
            }

            probSum = 1.0 / probSum;

            mu_p[j] = 0;
            for(int k = 0; k < ConSize; k++) {
                prob[j][k] *= probSum;

                // mu_p = prob * sym;
                mu_p[j] += prob[j][k] * Cons[k];
            }

            // sigma2_p=prob*(sym.^2)-mu_p.^2;
            sigma2_p[j] = -mu_p[j] * mu_p[j];

            for(int k = 0; k < ConSize; k++) {
                sigma2_p[j] += prob[j][k] * Cons2[k];
            }

            sigma2_p[j] = sigma2_p[j] < 5e-7 ? 2e6 : 1.0/sigma2_p[j];

            // tempAlpha = 1./sigma2_p - 1./ h2;
            double tempAlpha = 1.0 * sigma2_p[j] - 1.0 * h2[j];

            if (tempAlpha > 5e-7) {

                Alpha_new[j] = tempAlpha;
                Gamma_new[j] = mu_p[j] * sigma2_p[j] - t[j] * h2[j];

            }

            Alpha[j] = Alpha_new[j] * delta + Alpha[j] * (1 - delta);
            Gamma[j] = Gamma_new[j] * delta + Gamma[j] * (1 - delta);
        }

        for(int j = 0; j < TxAntNum2; j++) {
            HtHMod[j][j] = HtH[j][j] + Alpha[j];
        }

        choleskyInv -> execute(HtHMod, Sigma_q);

        for (int j = 0; j < TxAntNum2; j++) {
            HtRAddGamma[j] = HtR[j] + Gamma[j];
        }

        MatrixMultiplyVector(Sigma_q, HtRAddGamma, TxAntNum2, TxAntNum2, Mu_q);
    }

    symbolsToBits(Mu_q);

}