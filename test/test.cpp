#include "MMSE.h"
#include "utils.h"
#include "Detection.h"
#include "EP.h"

#include <iostream>
#include <chrono>

 

int main(){

    int TxAntNum = 16;
    int RxAntNum = 32;
    int ModType = 4;
    double SNRdB = 10;
    int sample = 1000;

    Detection * det = new DetectionCD(TxAntNum, RxAntNum, ModType, SNRdB);
    DetectionAlgorithm * alg = new MMSECD();

    alg->bind(det);

    auto start = std::chrono::high_resolution_clock::now();

    for(int i=0;i<sample;i++){
        det->generate();
        alg->execute();
        alg->check();
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Error bits: " << alg->getErrorBits() << std::endl;
    std::cout << "Error frames: " << alg->getErrorFrames() << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;


    return 0;
}   