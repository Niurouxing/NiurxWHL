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
    int TxAntNum = 32;
    int RxAntNum = 128;
    int ModType = 6;
    // std::vector<double> SNRdBVec =  {6,8,10,12,14};
    std::vector<double> SNRdBVec =  {10,12,14,16,18};
    // std::vector<double> SNRdBVec =  {14,16,18,20,22};
    int sample = 100000;

    std::vector<double> para = {0.6065062532342977,0.42134845302385443,0.8614308486933289,1.2618686744944017,0.6309871139655454,3.9386826676319653,-44.31916665429276,22.79765821972129,0.06694902448475279,-21.606207554330886};
    // para的前半部分是alphaVec, 后半部分是accuVec
    std::vector<double> alphaVec =  std::vector<double>(para.begin(), para.begin() + para.size() / 2);
    std::vector<double> accuVec =  std::vector<double>(para.begin() + para.size() / 2, para.end());

    for (double SNRdB : SNRdBVec)
    {
        auto start = std::chrono::system_clock::now();
        auto [errorBits, errorFrames] = EPAwNSADet(TxAntNum, RxAntNum, ModType, SNRdB, sample,0.9,para.size()/2,4,alphaVec,accuVec);
        auto end = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "SNRdB: " << SNRdB << " EPAwNSA: "
                  << "BER " << errorBits/ (double)sample  / (double)TxAntNum  / (double)ModType << " FER " << errorFrames / (double)sample << " time " << duration.count() << "ms" << std::endl;
    }

    return 0;
}