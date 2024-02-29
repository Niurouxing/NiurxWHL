#include "MMSE.h"
#include "utils.h"
#include "Detection.h"
#include "EP.h"
#include "EPAwNSA.h"
#include "ExBsP.h"

#include <iostream>
#include <chrono>

#include "ExBsP_NB.h"
#include "MIMO.h"
#include "NBLDPC.h"
#include <tuple>
#include "det.h"

 

int main()
{
    int TxAntNum = 64;
    int RxAntNum = 128;
    int ModType = 8;
    double SNRdB = 22;
    int sample = 10000;

    std::vector<double> alphaVec = {0.32635631700904194,0.34260822294608356,0.36463710293423696,0.43494135528449085,0.741037241639934,1.65820034989151,6.655754002682825,0.4695309686920808,0.7765440601519576,0.42918308014071743};
    std::vector<double> accuVec = {9.957297158885503,-22.130643039168127,-13.101448937902674,-0.16183823042648138,27.881817330280708,38.60772999554956,-1.703075055656286,55.62630399534433,47.73697181100642,97.61366981252652};

    // std::tuple<int, int> EPAwNSADet(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample,double delta,int NSAiter,int iter,std::vector<double> alphaVec,std::vector<double> accuVec) 

    auto start = std::chrono::system_clock::now();
    auto [errorBits, errorFrames] = EPAwNSADet(TxAntNum, RxAntNum, ModType, SNRdB, sample,0.9,10,7,alphaVec,accuVec);
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "EPAwNSA: "
              << "BER " << errorBits << std::endl;
}