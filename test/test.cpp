#include "MMSE.h"
#include "utils.h"
#include "Mimo.h"
#include "EP.h"

#include <iostream>
#include <chrono>

 

static std::mt19937 rng(123456);
static std::normal_distribution<double> distribution(0.0, 1.0);  
inline double randomGaussian() {
    return distribution(rng);
}

int main(){


    int row=5;
    int col=3;

    double ** matrix = new double*[row];
    for (int i = 0; i < row; i++) {
        matrix[i] = new double[col];
    }

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++){
            matrix[i][j] = randomGaussian();
        }
    }

    std::cout << "Matrix : " << std::endl;
    printMatrix(matrix, row, col);

    double ** matrixTMatrix = new double*[col];
    for (int i = 0; i < col; i++) {
        matrixTMatrix[i] = new double[col];
    }

    MatrixTransposeMultiplyMatrix(matrix, matrix, row, col, col, matrixTMatrix);

    std::cout << "MatrixTMatrix : " << std::endl;
    printMatrix(matrixTMatrix, col, col);

    CholeskyInv * choleskyInv = new CholeskyInv(col);

    double ** matrixTMatrixInv = new double*[col];
    for (int i = 0; i < col; i++) {
        matrixTMatrixInv[i] = new double[col];
    }

    choleskyInv->execute(matrixTMatrix, matrixTMatrixInv);

    std::cout << "MatrixTMatrixInv : " << std::endl;
    printMatrix(matrixTMatrixInv, col, col);




    Mimo::createMimo(8,16,4,10);

    std::cout << "Cons : " << std::endl;
    printVector(Mimo::getMimo()->Cons, Mimo::getMimo()->ConSize);

    std::cout << "bitCons : " << std::endl;
    printVector(Mimo::getMimo()->bitCons, Mimo::getMimo()->ConSize * Mimo::getMimo()->bitLength);

    Algorithm * alg = new EP(5,0.9);
    // Algorithm * alg = new MMSE();

    Mimo::getMimo()->reset();

    std::cout << "H : " << std::endl;
    printMatrix(Mimo::getMimo()->H, Mimo::getMimo()->RxAntNum2, Mimo::getMimo()->TxAntNum2);

    std::cout << "RxSymbols : " << std::endl;
    printVector(Mimo::getMimo()->RxSymbols, Mimo::getMimo()->RxAntNum2);

    std::cout << "TxSymbols : " << std::endl;
    printVector(Mimo::getMimo()->TxSymbols, Mimo::getMimo()->TxAntNum2);

    std::cout << "TxBits : " << std::endl;
    printMatrix(Mimo::getMimo()->TxBits, Mimo::getMimo()->TxAntNum2, Mimo::getMimo()->bitLength);

    alg->execute();

    std::cout << "TxBitsEst : " << std::endl;
    printMatrix(alg->getTxBitsEst(), Mimo::getMimo()->TxAntNum2, Mimo::getMimo()->bitLength);
    

    for (int i = 0; i < 10000; i++) {
        Mimo::getMimo()->reset();
        alg->execute();
        alg->check();
    }

    std::cout << "Error Bits : " << alg->getErrorBits() << std::endl;
    std::cout << "Error Frames : " << alg->getErrorFrames() << std::endl;

    return 0;
}