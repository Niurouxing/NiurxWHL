#include <cmath>

#include <chrono>
#include <random>

#include "Cons.h"
#include "utils.h"
#include "Detection.h"

std::mt19937 seedInit()
{
    std::random_device rd;
    auto seed = rd() ^ std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);
    return rng;
}

static std::mt19937 rng1 = seedInit();
static std::mt19937 rng2 = seedInit();
static std::normal_distribution<double> distribution(0.0, 1.0);
static std::normal_distribution<double> distribution_divSqrt2(0.0, std::sqrt(0.5));

inline double randomGaussian()
{
    return distribution(rng1);
}

inline double randomGaussian_divSqrt2()
{
    return distribution_divSqrt2(rng2);
}


Detection::Detection(int TxAntNum, int RxAntNum, int ModType, double SNRdB)
{
    this->TxAntNum = TxAntNum;
    this->RxAntNum = RxAntNum;
    this->TxAntNum2 = 2 * TxAntNum;
    this->RxAntNum2 = 2 * RxAntNum;
    this->ModType = ModType;
    this->SNRdB = SNRdB;

    this->TxBits = new int[TxAntNum * ModType];

    this->Nv = TxAntNum * RxAntNum / (pow(10, SNRdB / 10) * ModType * TxAntNum);
    this->sqrtNvDiv2 = std::sqrt(Nv / 2);
    this->NvInv = 1 / Nv;
}

void Detection::generate()
{
    generateChannel();
    generateTxSignals();
    generateRxSignalsWithNoise();
}

Detection::~Detection()
{
    delete[] TxBits;
}

DetectionRD::DetectionRD(int TxAntNum, int RxAntNum, int ModType, double SNRdB) : Detection(TxAntNum, RxAntNum, ModType, SNRdB)
{
    switch (ModType)
    {
    case 2:
        this->ConSize = 2;
        this->bitLength = 1;

        this->Cons = realConsMod2;
        this->Cons2 = realCons2Mod2;
        this->bitCons = realBitConsMod2;
        break;
    case 4:
        this->ConSize = 4;
        this->bitLength = 2;

        this->Cons = realConsMod4;
        this->Cons2 = realCons2Mod4;
        this->bitCons = realBitConsMod4;
        break;
    case 6:
        this->ConSize = 8;
        this->bitLength = 3;

        this->Cons = realConsMod6;
        this->Cons2 = realCons2Mod6;
        this->bitCons = realBitConsMod6;
        break;
    case 8:
        this->ConSize = 16;
        this->bitLength = 4;

        this->Cons = realConsMod8;
        this->Cons2 = realCons2Mod8;
        this->bitCons = realBitConsMod8;
        break;
    }

    this->TxSymbols = new double[TxAntNum2];
    this->RxSymbols = new double[RxAntNum2];

    this->H = new double *[RxAntNum2];
    for (int i = 0; i < RxAntNum2; i++)
    {
        this->H[i] = new double[TxAntNum2];
    }

    this->randomInt = std::uniform_int_distribution<int>(0, ConSize - 1);
}

void DetectionRD::generateChannel()
{

    static double randomTemp;

    for (int r = 0; r < RxAntNum; r++)
    {
        for (int t = 0; t < TxAntNum; t++)
        {
            randomTemp = randomGaussian_divSqrt2();
            H[r][t] = randomTemp;
            H[r + RxAntNum][t + TxAntNum] = randomTemp;

            randomTemp = randomGaussian_divSqrt2();
            H[r][t + TxAntNum] = -randomTemp;
            H[r + RxAntNum][t] = randomTemp;
        }
    }
}

void DetectionRD::generateTxSignals()
{

    static int index;

    for (int i = 0; i < TxAntNum2; i++)
    {
        index = randomInt(rng1);
        TxSymbols[i] = Cons[index];
        for (int b = 0; b < bitLength; b++)
        {
            TxBits[i * bitLength + b] = bitCons[index * bitLength + b];
        }
    }
}

void DetectionRD::generateRxSignalsWithNoise()
{

    for (int r = 0; r < RxAntNum2; r++)
    {
        RxSymbols[r] = 0;
        for (int t = 0; t < TxAntNum2; t++)
        {
            RxSymbols[r] += H[r][t] * TxSymbols[t];
        }
        RxSymbols[r] += randomGaussian() * sqrtNvDiv2;
    }
}

DetectionRD::~DetectionRD()
{
    for (int i = 0; i < RxAntNum2; i++)
    {
        delete[] H[i];
    }
    delete[] H;

    delete[] TxSymbols;
    delete[] RxSymbols;
}


DetectionCD::DetectionCD(int TxAntNum, int RxAntNum, int ModType, double SNRdB) : Detection(TxAntNum, RxAntNum, ModType, SNRdB)
{
    switch (ModType)
    {
    case 2:
        this->ConSize = 4;
        this->bitLength = 2;

        this->ConsReal = realConsMod2;
        this->bitConsReal = realBitConsMod2;
        this->ConsComplex = complexConsMod2;
        this->bitConsComplex = complexBitConsMod2;
        break;
    case 4:
        this->ConSize = 16;
        this->bitLength = 4;

        this->ConsReal = realConsMod4;
        this->bitConsReal = realBitConsMod4;
        this->ConsComplex = complexConsMod4;
        this->bitConsComplex = complexBitConsMod4;
        break;
    case 6:
        this->ConSize = 64;
        this->bitLength = 6;

        this->ConsReal = realConsMod6;
        this->bitConsReal = realBitConsMod6;
        this->ConsComplex = complexConsMod6;
        this->bitConsComplex = complexBitConsMod6;
        break;
    case 8:
        this->ConSize = 256;
        this->bitLength = 8;

        this->ConsReal = realConsMod8;
        this->bitConsReal = realBitConsMod8;
        this->ConsComplex = complexConsMod8;
        this->bitConsComplex = complexBitConsMod8;
        break;
    }

    this->TxIndiceCD = new int[TxAntNum];
    this->TxSymbols = new std::complex<double>[TxAntNum];
    this->RxSymbols = new std::complex<double>[RxAntNum];

    this->H = new std::complex<double> *[RxAntNum];
    for (int i = 0; i < RxAntNum; i++)
    {
        this->H[i] = new std::complex<double>[TxAntNum];
    }

    this->randomInt = std::uniform_int_distribution<int>(0, ConSize - 1);
}

void DetectionCD::generateChannel()
{

    static double randomTemp1, randomTemp2;

    for (int r = 0; r < RxAntNum; r++)
    {
        for (int t = 0; t < TxAntNum; t++)
        {
            randomTemp1 = randomGaussian_divSqrt2();
            randomTemp2 = randomGaussian_divSqrt2();
            H[r][t] = std::complex<double>(randomTemp1, randomTemp2); 
        }
    }
}

void DetectionCD::generateTxSignals()
{

    static int index;

    for (int i = 0; i < TxAntNum; i++)
    {
        index = randomInt(rng1);
        TxSymbols[i] = ConsComplex[index];
        TxIndiceCD[i] = index;
        for (int b = 0; b < bitLength; b++)
        {
            TxBits[i * bitLength + b] = bitConsComplex[index * bitLength + b];
        }
    }
}

void DetectionCD::generateRxSignalsWithNoise()
{

    for (int r = 0; r < RxAntNum; r++)
    {
        RxSymbols[r] = 0;
        for (int t = 0; t < TxAntNum; t++)
        {
            RxSymbols[r] += H[r][t] * TxSymbols[t];
        }
        RxSymbols[r] += std::complex<double>(randomGaussian(), randomGaussian()) * sqrtNvDiv2;
    }

}

DetectionCD::~DetectionCD()
{
    for (int i = 0; i < RxAntNum; i++)
    {
        delete[] H[i];
    }
    delete[] H;

    delete[] TxIndiceCD;
    delete[] TxSymbols;
    delete[] RxSymbols;
}