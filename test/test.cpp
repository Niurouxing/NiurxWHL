#include "MMSE.h"
#include "utils.h"
#include "Detection.h"
#include "EP.h"
#include "ExBsP.h"

#include <iostream>
#include <chrono>

#include "ExBsP_NB.h"
#include "MIMO.h"
#include "NBLDPC.h"

int main()
{

    openblas_set_num_threads(1);
    int TxAntNum = 4;
    int RxAntNum = 8;
    int ModType = 6;
    double SNRdB = 7;
    int sample = 10;

    // Detection * det = new DetectionRD(TxAntNum, RxAntNum, ModType, SNRdB);
    // DetectionAlgorithm * alg = new EP(5,0.9);

    // alg->bind(det);

    // auto start = std::chrono::high_resolution_clock::now();

    // for(int i=0;i<sample;i++){
    //     det->generate();
    //     alg->execute();
    //     alg->check();
    // }

    // auto end = std::chrono::high_resolution_clock::now();

    // std::cout << "Error bits: " << alg->getErrorBits() << std::endl;
    // std::cout << "Error frames: " << alg->getErrorFrames() << std::endl;
    // std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

    NBLDPC *nbldpc = new NBLDPC(96, 48, 64, 20, 20);

    // print code->mat, code->matValue
    for (int i = 0; i < nbldpc->code->M; i++)
    {
        for (int j = 0; j < nbldpc->code->rowDegree[i]; j++)
        {
            std::cout << nbldpc->code->mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    for (int i = 0; i < nbldpc->code->M; i++)
    {
        for (int j = 0; j < nbldpc->code->rowDegree[i]; j++)
        {
            std::cout << nbldpc->code->matValue[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print matUT
    for (int i = 0; i < nbldpc->code->M; i++)
    {
        for (int j = 0; j < nbldpc->code->N; j++)
        {
            std::cout << nbldpc->code->matUT[i][j] << " ";
        }
        std::cout << std::endl;
    }

    auto mimo = MIMO::getMIMO();
    mimo->addCode(nbldpc);
    mimo->addDetection(true, TxAntNum, RxAntNum, ModType, SNRdB);

    std::cout << "Code length: " << nbldpc->codeLength << std::endl;
    std::cout << "Block number: " << mimo->blockNum << std::endl;

    ExBsP_NB *exbsp = new ExBsP_NB(5, 1,1, 100, 0.3);

    std::cout << "bind success" << std::endl;

    mimo->generate();

    std::cout << "generate success" << std::endl;

    // codeWords
    std::cout << "codeWords: " << std::endl;
    for (int i = 0; i < exbsp->getCode()->code->N; i++)
    {
        std::cout << exbsp->getCode()->codeWords[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MIMO::getMIMO()->blockNum; i++)
    {
        std::cout << "block " << i << std::endl;
        for (int j = 0; j < mimo->detections[i]->TxAntNum; j++)
        {
            for (int k = 0; k < mimo->detections[i]->ModType; k++)
            {
                std::cout << mimo->detections[i]->TxBits[j * mimo->detections[i]->ModType + k] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::cout << "Tx Symbols: " << std::endl;
    for (int i = 0; i < MIMO::getMIMO()->blockNum; i++)
    {
        std::cout << "block " << i << std::endl;
        for (int j = 0; j < mimo->detections[i]->TxAntNum; j++)
        {

            std::cout << dynamic_cast<DetectionCD *>(mimo->detections[i])->TxSymbols[j] << " ";

            std::cout << std::endl;
        }
    }

    std::cout << "Rx Symbols: " << std::endl;
    for (int i = 0; i < MIMO::getMIMO()->blockNum; i++)
    {
        std::cout << "block " << i << std::endl;
        for (int j = 0; j < mimo->detections[i]->RxAntNum; j++)
        {

            std::cout << dynamic_cast<DetectionCD *>(mimo->detections[i])->RxSymbols[j] << " ";

            std::cout << std::endl;
        }
    }

    exbsp->execute();

    std::cout << "decide :" << std::endl;
    for (int i = 0; i < 16; i++)
    {
        std::cout << exbsp->getCode()->decide[i] << ' ';
    }

    std::cout << std::endl;

    std::cout << "error bits: " << exbsp->getCode()->getErrorBits() << std::endl;

    return 0;
}