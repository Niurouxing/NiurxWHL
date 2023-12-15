#include "utils.h"

arma::mat GenerateChannelMatrix(int TxAntNum, int RxAntNum)
{

    arma::mat HReal(RxAntNum, TxAntNum, arma::fill::randn);
    arma::mat HImag(RxAntNum, TxAntNum, arma::fill::randn);

    // multiply by 1/sqrt(2) to get unit variance
    HReal *= 1.0 / std::sqrt(2.0);
    HImag *= 1.0 / std::sqrt(2.0);

    // [HReal, -HImag; HImag, HReal]
    arma::mat H(2 * RxAntNum, 2 * TxAntNum, arma::fill::zeros);
    H.submat(0, 0, RxAntNum - 1, TxAntNum - 1) = HReal;
    H.submat(RxAntNum, TxAntNum, 2 * RxAntNum - 1, 2 * TxAntNum - 1) = HReal;
    H.submat(0, TxAntNum, RxAntNum - 1, 2 * TxAntNum - 1) = -HImag;
    H.submat(RxAntNum, 0, 2 * RxAntNum - 1, TxAntNum - 1) = HImag;

    return H;
}

std::tuple<arma::vec, arma::mat> GenerateCons(int ModType)
{
    double norm_f;
    int ConSize;
    int bitLength;

    if (ModType == 2)
    {
        norm_f = 1.0 / std::sqrt(2.0);
        ConSize = 2;
        bitLength = 1;
    }
    else if (ModType == 4)
    {
        norm_f = 1.0 / std::sqrt(10.0);
        ConSize = 4;
        bitLength = 2;
    }
    else if (ModType == 6)
    {
        norm_f = 1.0 / std::sqrt(42.0);
        ConSize = 8;
        bitLength = 3;
    }
    else if (ModType == 8)
    {
        norm_f = 1.0 / std::sqrt(170.0);
        ConSize = 16;
        bitLength = 4;
    }

    arma::vec Cons(ConSize, arma::fill::zeros);
    arma::mat bitCons(ConSize, bitLength, arma::fill::zeros);

    if (ModType == 2)
    {
        Cons = {-norm_f, norm_f};
        bitCons[1] = 1;
    }
    else if (ModType == 4)
    {
        Cons = {-norm_f, -3 * norm_f, norm_f, 3 * norm_f};
        bitCons = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    }
    else if (ModType == 6)
    {
        Cons = {-3 * norm_f, -norm_f, -5 * norm_f, -7 * norm_f, norm_f, 3 * norm_f, 5 * norm_f, 7 * norm_f};
        bitCons = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};
    }
    else if (ModType == 8)
    {
        Cons = {-5 * norm_f, -7 * norm_f, -3 * norm_f, -norm_f, -11 * norm_f, -9 * norm_f, -13 * norm_f, -15 * norm_f, 5 * norm_f, 7 * norm_f, 3 * norm_f, norm_f, 11 * norm_f, 9 * norm_f, 13 * norm_f, 15 * norm_f};
        bitCons = {{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 1, 1}, {0, 1, 0, 0}, {0, 1, 0, 1}, {0, 1, 1, 0}, {0, 1, 1, 1}, {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 1, 0}, {1, 0, 1, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 0}, {1, 1, 1, 1}};
    }

    return std::make_tuple(Cons, bitCons);
}

std::tuple<arma::uvec, arma::vec> GenerateTxSignals(int TxAntNum, arma::vec Cons)
{

    int ConSize = Cons.n_elem;

    auto TxIndice = arma::randi<arma::uvec>(2 * TxAntNum, arma::distr_param(0, ConSize - 1));

    auto TxSignals = Cons.elem(TxIndice);

    return std::make_tuple(TxIndice, TxSignals);
}

std::tuple<arma::vec, double> GenerateRxSignals(arma::mat H, arma::vec TxSignals, double SNRdB, int ModType)
{

    int RxAntNum = H.n_rows / 2;
    int TxAntNum = H.n_cols / 2;
    double SNR = std::pow(10, SNRdB / 10);
    double Nv = RxAntNum / (SNR * ModType);

    arma::vec Noise = arma::randn<arma::vec>(2 * RxAntNum);
    Noise *= std::sqrt(Nv / 2);

    arma::vec RxSignals = H * TxSignals + Noise;

    return std::make_tuple(RxSignals, Nv);
}

std::tuple<arma::uvec, arma::vec> MMSE(arma::mat H, arma::vec RxSignals, arma::vec Cons, double Nv)
{

    arma::mat HtH = H.t() * H;
    HtH.diag() += Nv;

    arma::vec TxEst = arma::inv(HtH) * H.t() * RxSignals;

    arma::mat dis = arma::abs(arma::repmat(TxEst, 1, Cons.n_elem) - arma::repmat(Cons, 1, TxEst.n_elem).t());

    arma::uvec RxIndice = arma::index_min(dis, 1);

    arma::vec RxSymbols = Cons.elem(RxIndice);

    return std::make_tuple(RxIndice, RxSymbols);
}

std::tuple<arma::uvec, arma::vec> EP(arma::mat H, arma::vec RxSignals, arma::vec Cons, double Nv, int iter)
{
    double delta=0.9;
    // initialize
    int RxAntNum2 = H.n_rows;
    int TxAntNum2 = H.n_cols;

    arma::vec Alpha = 2 * arma::vec(TxAntNum2, arma::fill::ones);
    arma::vec Gamma = arma::vec(TxAntNum2, arma::fill::zeros);

    arma::vec tempAlpha = arma::vec(TxAntNum2, arma::fill::zeros);
    arma::vec tempGamma = arma::vec(TxAntNum2, arma::fill::zeros);


    arma::vec Alpha_new = arma::vec(TxAntNum2, arma::fill::zeros);
    arma::vec Gamma_new = arma::vec(TxAntNum2, arma::fill::zeros);


    arma::vec sig = arma::vec(TxAntNum2, arma::fill::zeros);
    arma::vec h2 = arma::vec(TxAntNum2, arma::fill::zeros);
    arma::vec t = arma::vec(TxAntNum2, arma::fill::zeros);
    arma::vec sigma2_p = arma::vec(TxAntNum2, arma::fill::zeros);
    arma::vec mu_p = arma::vec(TxAntNum2, arma::fill::zeros);

    arma::mat prob = arma::mat(TxAntNum2, Cons.n_elem, arma::fill::zeros);

    arma::vec Cons2 = arma::square(Cons);

    arma::uvec in=arma::uvec(TxAntNum2, arma::fill::zeros);

    // perform MMSE
    arma::mat HtH = H.t() * H;
    HtH /= Nv;

    arma::mat HtHOne=HtH;
    HtHOne.diag() += 1;

    arma::mat Sigma_q = arma::inv(HtHOne);
    arma::vec Mu_q = Sigma_q * (H.t() * RxSignals / Nv);

    for (int i = 0; i < iter; i++)
    {
        sig = arma::diagvec(Sigma_q);
        h2 = sig / (1 - sig % Alpha);
        t = h2 % (Mu_q / sig - Gamma);

        prob = arma::exp(-arma::square(arma::repmat(t, 1, Cons.n_elem) - arma::repmat(Cons, 1, TxAntNum2).t()));
        prob.each_row() /=(2 * h2);
        prob.each_row() /= arma::sum(prob, 1);

        mu_p = arma::sum(prob.each_row() % Cons, 1);

        sigma2_p = arma::sum(prob.each_row() % Cons2, 1) - arma::square(mu_p);

        sigma2_p.elem(arma::find(sigma2_p < 5e-7)).fill(5e-7);

        tempAlpha=1/sigma2_p-1/h2;
        tempGamma=mu_p/sigma2_p-t/h2;

        in=arma::find(tempAlpha > 5e-7);
        Alpha_new.elem(in)=tempAlpha.elem(in);
        Gamma_new.elem(in)=tempGamma.elem(in);

        Alpha=delta*Alpha_new+(1-delta)*Alpha;
        Gamma=delta*Gamma_new+(1-delta)*Gamma;

        HtHOne=HtH;

        HtHOne.diag()+=Alpha;
        Sigma_q = arma::inv(HtHOne);
        Mu_q = Sigma_q * (H.t() * RxSignals / Nv);
    }

    arma::mat dis = arma::abs(arma::repmat(Mu_q, 1, Cons.n_elem) - arma::repmat(Cons, 1, TxAntNum2).t());

    arma::uvec RxIndice = arma::index_min(dis, 1);

    arma::vec RxSymbols = Cons.elem(RxIndice);

    return std::make_tuple(RxIndice, RxSymbols);

}

int Detection(int TxAntNum, int RxAntNum, double SNRdB, int ModType, int sample)
{

    int error = 0;
    auto [Cons, bitCons] = GenerateCons(ModType);

    for (int i = 0; i < sample; i++)
    {
        auto H = GenerateChannelMatrix(TxAntNum, RxAntNum);

        auto [TxIndice, TxSignals] = GenerateTxSignals(TxAntNum, Cons);

        auto [RxSignals, Nv] = GenerateRxSignals(H, TxSignals, SNRdB, ModType);

        auto [RxIndice, RxSymbols] = MMSE(H, RxSignals, Cons, Nv);

        arma::mat TxBits = bitCons.rows(TxIndice);
        arma::mat RxBits = bitCons.rows(RxIndice);

        error += arma::accu(TxBits != RxBits);
    }

    return error;
}
