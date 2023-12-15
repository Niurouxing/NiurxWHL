#include "utils.h"

arma::mat GenerateChannelMatrix(int TxAntNum, int RxAntNum)
{

    static arma::mat HReal(RxAntNum, TxAntNum);
    static arma::mat HImag(RxAntNum, TxAntNum);

    HReal.randn();  
    HImag.randn();

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
    static arma::uvec TxIndice(2 * TxAntNum);
    static arma::vec TxSignals(2 * TxAntNum);

    static int ConSize = Cons.n_elem;

    TxIndice = arma::randi<arma::uvec>(2 * TxAntNum, arma::distr_param(0, ConSize - 1));

    TxSignals = Cons.elem(TxIndice);

    return std::make_tuple(TxIndice, TxSignals);
}

std::tuple<arma::vec, double> GenerateRxSignals(arma::mat H, arma::vec TxSignals, double SNRdB, int ModType)
{

    static int RxAntNum = H.n_rows / 2;
    static int TxAntNum = H.n_cols / 2;
    static double SNR = std::pow(10, SNRdB / 10);
    static double Nv = RxAntNum / (SNR * ModType);
    static double sqrtNvd2 = std::sqrt(Nv / 2);

    static arma::vec Noise(RxAntNum * 2);
    Noise.randn();
    Noise *= sqrtNvd2;

    static arma::vec RxSignals(RxAntNum * 2);
    RxSignals = H * TxSignals + Noise;

    return std::make_tuple(RxSignals, Nv);
}

std::tuple<arma::uvec, arma::vec> MMSE(arma::mat H, arma::vec RxSignals, arma::vec Cons, double Nv)
{

    arma::mat HtH = H.t() * H;
    HtH.diag() += Nv;

    arma::vec TxEst = arma::inv(HtH) * H.t() * RxSignals;

    arma::mat dis = arma::abs(arma::repmat(TxEst, 1, Cons.n_elem) - arma::repmat(Cons, 1, TxEst.n_elem).t());

    arma::uvec TxBitsEst = arma::index_min(dis, 1);

    arma::vec TxSymbolsEst = Cons.elem(TxBitsEst);

    return std::make_tuple(TxBitsEst, TxSymbolsEst);
}

std::tuple<arma::uvec, arma::vec> EP(arma::mat H, arma::vec RxSignals, arma::vec Cons, double Nv, int iter)
{
    static double delta = 0.9;

    static int RxAntNum2 = H.n_rows;
    static int TxAntNum2 = H.n_cols;

    static arma::vec Alpha(TxAntNum2); 
    Alpha = 2 * arma::ones<arma::vec>(TxAntNum2);
    static arma::vec Gamma(TxAntNum2);
    Gamma.zeros();

    static arma::vec Alpha_new(TxAntNum2);
    static arma::vec Gamma_new(TxAntNum2);

    static arma::vec sig(TxAntNum2);
    static arma::vec h2(TxAntNum2);
    static arma::vec t(TxAntNum2);

    static arma::vec Cons2 = arma::square(Cons);
    static arma::mat ConsTRep = arma::repmat(Cons.t(), TxAntNum2, 1);

    static arma::mat prob(TxAntNum2, Cons.n_elem);
    static arma::vec sigma2_p(TxAntNum2);

    static arma::uvec in;

    static arma::mat HtH(TxAntNum2, TxAntNum2);
    HtH = H.t() * H;
    HtH /= Nv;

    static arma::mat HtHMod(TxAntNum2, TxAntNum2);
    HtHMod = HtH;
    HtHMod.diag() += 1;

    static arma::vec HtRxSignals(TxAntNum2);
    HtRxSignals = H.t() * RxSignals / Nv;

    static arma::vec Mu_q(TxAntNum2);
    Mu_q = arma::solve(HtHMod, HtRxSignals, arma::solve_opts::fast);

    for (int i = 0; i < iter; ++i) {
        sig = arma::diagvec(HtHMod - HtH);
        h2 = 1.0 / (1 / sig - Alpha);
        t = h2 % (Mu_q / sig - Gamma);

        prob = arma::exp(-arma::square(t * arma::ones<arma::rowvec>(Cons.n_elem) - ConsTRep));
        prob.each_col() /= (2 * h2);
        prob.each_col() /= arma::sum(prob, 1);

        sigma2_p = arma::sum(prob % arma::square(ConsTRep), 1) - arma::square(arma::sum(prob % ConsTRep, 1));
        sigma2_p.elem(arma::find(sigma2_p < 5e-7)).fill(5e-7);

        Alpha_new = 1 / sigma2_p - 1 / h2;
        Gamma_new = arma::sum(prob % ConsTRep, 1) / sigma2_p - t / h2;

        in = arma::find(Alpha_new > 5e-7);
        Alpha.elem(in) = delta * Alpha_new.elem(in) + (1 - delta) * Alpha.elem(in);
        Gamma.elem(in) = delta * Gamma_new.elem(in) + (1 - delta) * Gamma.elem(in);

        HtHMod.diag() = HtH.diag() + Alpha;
        Mu_q = arma::solve(HtHMod, HtRxSignals, arma::solve_opts::fast);
    }

    static arma::mat dis(TxAntNum2, Cons.n_elem);
    dis = arma::abs(Mu_q * arma::ones<arma::rowvec>(Cons.n_elem) - ConsTRep);

    static arma::uvec TxIndiceEst(TxAntNum2);
    TxIndiceEst = arma::index_min(dis, 1);

    static arma::vec TxSymbolsEst(TxAntNum2);
    TxSymbolsEst = Cons.elem(TxIndiceEst);

    return std::make_tuple(TxIndiceEst, TxSymbolsEst);
}


int Detection(int TxAntNum, int RxAntNum, double SNRdB, int ModType, int sample)
{
    int error = 0;
    arma::vec Cons;
    arma::mat bitCons;
    std::tie(Cons, bitCons) = GenerateCons(ModType);

    arma::mat H, TxBits, TxBitsEst;
    arma::uvec TxIndice, TxIndiceEst;
    arma::vec TxSignals, RxSignals;
    double Nv;

    for (int i = 0; i < sample; i++)
    {
        H = GenerateChannelMatrix(TxAntNum, RxAntNum);

        std::tie(TxIndice, TxSignals) = GenerateTxSignals(TxAntNum, Cons);

        std::tie(RxSignals, Nv) = GenerateRxSignals(H, TxSignals, SNRdB, ModType);

        std::tie(TxIndiceEst, std::ignore) = EP(H, RxSignals, Cons, Nv, 5);

        TxBits = bitCons.rows(TxIndice);
        TxBitsEst = bitCons.rows(TxIndiceEst);

        error += arma::accu(TxBits != TxBitsEst);
    }

    return error;
}

