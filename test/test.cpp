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
    double SNRdB = 18;
    int sample = 100000;

    std::vector<double> alphaVec = {0.23478529776512855,2.7482475294009006,0.7876444963359318,0.3177415931724516,6.672742701679837,0.5158914335245294,0.5523547456356083,0.3560593542330146,0.40290725469817174,0.5569215000263725};
    std::vector<double> accuVec = {7.58849286951194,19.152180685104764,-0.21186167766993272,0.38862635614792684,24.006274316051897,9.970199793366808,7.388816330950275,0.13749695641325727,0.16301707018457567,38.373480083106784};

    // std::tuple<int, int> EPAwNSADet(int TxAntNum, int RxAntNum, int ModType, double SNRdB, int sample,double delta,int NSAiter,int iter,std::vector<double> alphaVec,std::vector<double> accuVec) 

    auto start = std::chrono::system_clock::now();
    auto [errorBits, errorFrames] = EPAwNSADet(TxAntNum, RxAntNum, ModType, SNRdB, sample,0.9,10,7,alphaVec,accuVec);
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "EPAwNSA: "
              << "BER " << errorBits << std::endl;
}