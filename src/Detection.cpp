#include <cmath>

#include <chrono>
 

#include "utils.h"
#include "Detection.h"

std::mt19937 seedInit() {
    std::random_device rd;
    auto seed = rd() ^ std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);
    return rng;
}


static std::mt19937 rng1 = seedInit();
static std::mt19937 rng2 = seedInit();
static std::normal_distribution<double> distribution(0.0, 1.0);  
static std::normal_distribution<double> distribution_divSqrt2(0.0, std::sqrt(0.5));  

inline double randomGaussian() {
    return distribution(rng1);
}

inline double randomGaussian_divSqrt2() {
    return distribution_divSqrt2(rng2);
}

inline int randomInt(int max) {
    return rand() % max;
}

static constexpr double ConsMod2[2] = {0.70710678118, -0.70710678118};
static constexpr double ConsMod4[4] = {-0.31622776601683794, -0.9486832980505138, 0.31622776601683794, 0.9486832980505138};
static constexpr double ConsMod6[8] = {-0.4629100498862757, -0.1543033499620919, -0.7715167498104595, -1.0801234497346432, 0.1543033499620919, 0.4629100498862757, 0.7715167498104595, 1.0801234497346432};
static constexpr double ConsMod8[16] = {-0.3834824944236852, -0.5368754921931592, -0.2300894966542111, -0.07669649888473704, -0.8436614877321074, -0.6902684899626333, -0.9970544855015815, -1.1504474832710556, 0.3834824944236852, 0.5368754921931592, 0.2300894966542111, 0.07669649888473704, 0.8436614877321074, 0.6902684899626333, 0.9970544855015815, 1.1504474832710556};

static constexpr double Cons2Mod2[2] = {0.5, 0.5};
static constexpr double Cons2Mod4[4] = {0.1, 0.9, 0.1, 0.9};
static constexpr double Cons2Mod6[8] = {0.214285714285714,   0.023809523809524,   0.595238095238095,   1.166666666666666,   0.023809523809524,  0.214285714285714,   0.595238095238095,   1.166666666666666};
static constexpr double Cons2Mod8[16] = {0.147058823529412,   0.288235294117647,   0.052941176470588,   0.005882352941176,   0.711764705882353,  0.476470588235294,   0.994117647058823,   1.323529411764706,   0.147058823529412,   0.288235294117647,   0.052941176470588,   0.005882352941176,   0.711764705882353,   0.476470588235294,   0.994117647058823,   1.323529411764706};


static constexpr int bitConsMod2[2] = {0, 1};
static constexpr int bitConsMod4[8] = {0, 0, 
                                       0, 1,
                                       1, 0,
                                       1, 1};
static constexpr int bitConsMod6[24] = {0, 0, 0,
                                        0, 0, 1,
                                        0, 1, 0,
                                        0, 1, 1,
                                        1, 0, 0,
                                        1, 0, 1,
                                        1, 1, 0,
                                        1, 1, 1};
static constexpr int bitConsMod8[64] = {0, 0, 0, 0,
                                        0, 0, 0, 1,
                                        0, 0, 1, 0,
                                        0, 0, 1, 1,
                                        0, 1, 0, 0,
                                        0, 1, 0, 1,
                                        0, 1, 1, 0,
                                        0, 1, 1, 1,
                                        1, 0, 0, 0,
                                        1, 0, 0, 1,
                                        1, 0, 1, 0,
                                        1, 0, 1, 1,
                                        1, 1, 0, 0,
                                        1, 1, 0, 1,
                                        1, 1, 1, 0,
                                        1, 1, 1, 1};

Detection* Detection::detection = nullptr;

void Detection::createDetection(int TxAntNum, int RxAntNum, int ModType, double SNRdB){
    if (detection == nullptr) {
        detection = new Detection(TxAntNum, RxAntNum, ModType, SNRdB);
    }
}

Detection * Detection::getDetection(){
    return detection;
}


Detection::Detection(int TxAntNum, int RxAntNum, int ModType, double SNRdB){
    switch (ModType) {
        case 2:
            this->ConSize=2;
            this->bitLength=1;

            this->Cons = ConsMod2;
            this->Cons2 = Cons2Mod2;
            this->bitCons = bitConsMod2;
            break;
        case 4:
            this->ConSize=4;
            this->bitLength=2;

            this->Cons = ConsMod4;
            this->Cons2 = Cons2Mod4;
            this->bitCons = bitConsMod4;
            break;
        case 6:
            this->ConSize=8;
            this->bitLength=3;

            this->Cons = ConsMod6;
            this->Cons2 = Cons2Mod6;
            this->bitCons = bitConsMod6;
            break;
        case 8:
            this->ConSize=16;
            this->bitLength=4;

            this->Cons = ConsMod8;
            this->Cons2 = Cons2Mod8;
            this->bitCons = bitConsMod8;
            break;
    }
    this->TxAntNum = TxAntNum;
    this->RxAntNum = RxAntNum;
    this->TxAntNum2 = 2 * TxAntNum;
    this->RxAntNum2 = 2 * RxAntNum;
    this->ModType = ModType;
    this->SNRdB = SNRdB;

    this->TxSymbols = new double[TxAntNum2];
    this->RxSymbols = new double[RxAntNum2];

    this->TxBits = new int[TxAntNum2 * bitLength];

    this->H = new double*[RxAntNum2];
    for (int i = 0; i < RxAntNum2; i++) {
        this->H[i] = new double[TxAntNum2];
    }

    this->Nv = TxAntNum * RxAntNum / (pow(10, SNRdB / 10) * ModType * TxAntNum);
    this->sqrtNvDiv2 = std::sqrt(Nv / 2);
    this->NvInv = 1 / Nv;
}



void Detection::generateChannel(){

    static double randomTemp;

    for (int r = 0; r < RxAntNum; r++) {
        for (int t = 0; t < TxAntNum; t++) {
            randomTemp = randomGaussian_divSqrt2();
            H[r][t] = randomTemp;
            H[r+RxAntNum][t+TxAntNum] = randomTemp;

            randomTemp = randomGaussian_divSqrt2();
            H[r][t+TxAntNum] = -randomTemp;
            H[r+RxAntNum][t] = randomTemp;
        }
    }
}

void Detection::generateTxSignals(){

    static int index;

    for (int i = 0; i < TxAntNum2; i++) {
        index = randomInt(ConSize);
        TxSymbols[i] = Cons[index];
        for (int b = 0; b < bitLength; b++) {
            TxBits[i * bitLength + b] = bitCons[index * bitLength + b];
        }
    }
}

void Detection::generateRxSignalsWithNoise(){

    for (int r = 0; r < RxAntNum2; r++) {
        RxSymbols[r] = 0;
        for (int t = 0; t < TxAntNum2; t++) {
            RxSymbols[r] += H[r][t] * TxSymbols[t];
        }
        RxSymbols[r] += randomGaussian() * sqrtNvDiv2;
    }
}

void Detection::reset(){
    generateChannel();
    generateTxSignals();
    generateRxSignalsWithNoise();
}