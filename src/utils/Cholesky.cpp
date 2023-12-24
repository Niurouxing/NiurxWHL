#include <cstring>
#include <cmath>
 

#include "CholeskyInv.h"



CholeskyInv::CholeskyInv(int size, bool isComplex) {
    this->size = size;
    this->isComplex = isComplex;

    if (isComplex) {
        lComplex = new std::complex<double>[size * size];
        yComplex = new std::complex<double>[size];
    } else {
        lReal = new double[size * size];
        yReal = new double[size];
    }
}


void CholeskyInv::execute(double * A, double * AInv) {
    std::memset(lReal, 0, size * size * sizeof(double));
    
    // Cholesky 分解
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += lReal[i * size + k] * lReal[j * size + k];
            }
            lReal[i * size + j] = (i == j) ? std::sqrt(A[i * size + i] - sum) : (1.0 / lReal[j * size + j] * (A[i * size + j] - sum));
        }
    }

    // 求逆
    for (int i = 0; i < size; ++i) {
        memset(yReal, 0, size * sizeof(double));

        // 前向替换求解 L * y = e_i
        for (int j = 0; j < size; ++j) {
            yReal[j] = (j == i) ? 1.0 : 0.0;
            for (int k = 0; k < j; ++k) {
                yReal[j] -= lReal[j * size + k] * yReal[k];
            }
            yReal[j] /= lReal[j * size + j];
        }

        // 后向替换求解 L^T * x = y
        for (int j = size - 1; j >= 0; --j) {
            AInv[j * size + i] = yReal[j];
            for (int k = j + 1; k < size; ++k) {
                AInv[j * size + i] -= lReal[k * size + j] * AInv[k * size + i];
            }
            AInv[j * size + i] /= lReal[j * size + j];
        }
    }

}

void CholeskyInv::execute(double ** A, double ** AInv) {
    std::memset(lReal, 0, size * size * sizeof(double));

    // Cholesky 分解
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += lReal[i * size + k] * lReal[j * size + k];
            }
            lReal[i * size + j] = (i == j) ? std::sqrt(A[i][j] - sum) : (1.0 / lReal[j * size + j] * (A[i][j] - sum));
        }
    }

    // 求逆
    for (int i = 0; i < size; ++i) {
        memset(yReal, 0, size * sizeof(double));

        // 前向替换求解 L * y = e_i
        for (int j = 0; j < size; ++j) {
            yReal[j] = (j == i) ? 1.0 : 0.0;
            for (int k = 0; k < j; ++k) {
                yReal[j] -= lReal[j * size + k] * yReal[k];
            }
            yReal[j] /= lReal[j * size + j];
        }

        // 后向替换求解 L^T * x = y
        for (int j = size - 1; j >= 0; --j) {
            AInv[j][i] = yReal[j];
            for (int k = j + 1; k < size; ++k) {
                AInv[j][i] -= lReal[k * size + j] * AInv[k][i];
            }
            AInv[j][i] /= lReal[j * size + j];
        }
    }
}

void CholeskyInv::execute(std::complex<double> ** a, std::complex<double> ** aInv) {
    std::memset(lComplex, 0, size * size * sizeof(std::complex<double>));

    // Cholesky 分解 - 复数版
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j <= i; ++j) {
            std::complex<double> sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += lComplex[i * size + k] * std::conj(lComplex[j * size + k]);
            }
            if (i == j) {
                lComplex[i * size + j] = std::sqrt(a[i][j] - sum);
            } else {
                lComplex[i * size + j] = (1.0 / lComplex[j * size + j] * (a[i][j] - sum));
            }
        }
    }

    // 求逆 - 复数版
    for (int i = 0; i < size; ++i) {
        std::memset(yComplex, 0, size * sizeof(std::complex<double>));

        // 前向替换求解 L * y = e_i
        for (int j = 0; j < size; ++j) {
            yComplex[j] = (j == i) ? 1.0 : 0.0;
            for (int k = 0; k < j; ++k) {
                yComplex[j] -= lComplex[j * size + k] * yComplex[k];
            }
            yComplex[j] /= lComplex[j * size + j];
        }

        // 后向替换求解 L^T * x = y
        for (int j = size - 1; j >= 0; --j) {
            aInv[j][i] = yComplex[j];
            for (int k = j + 1; k < size; ++k) {
                aInv[j][i] -= std::conj(lComplex[k * size + j]) * aInv[k][i];
            }
            aInv[j][i] /= std::conj(lComplex[j * size + j]);
        }
    }
}

CholeskyInv::~CholeskyInv() {
    if (isComplex) {
        delete[] lComplex;
        delete[] yComplex;
    } else {
        delete[] lReal;
        delete[] yReal;
    }
}