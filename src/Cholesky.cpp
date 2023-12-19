#include <cstring>
#include <cmath>

#include "CholeskyInv.h"



CholeskyInv::CholeskyInv(int size) {
    this->size = size;
    L = new double[size * size];
    y = new double[size];
}

CholeskyInv::~CholeskyInv() {
    delete[] L;
    delete[] y;
}

void CholeskyInv::execute(double * A, double * AInv) {
    std::memset(L, 0, size * size * sizeof(double));
    
    // Cholesky 分解
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += L[i * size + k] * L[j * size + k];
            }
            L[i * size + j] = (i == j) ? std::sqrt(A[i * size + i] - sum) : (1.0 / L[j * size + j] * (A[i * size + j] - sum));
        }
    }

    // 求逆
    for (int i = 0; i < size; ++i) {
        memset(y, 0, size * sizeof(double));

        // 前向替换求解 L * y = e_i
        for (int j = 0; j < size; ++j) {
            y[j] = (j == i) ? 1.0 : 0.0;
            for (int k = 0; k < j; ++k) {
                y[j] -= L[j * size + k] * y[k];
            }
            y[j] /= L[j * size + j];
        }

        // 后向替换求解 L^T * x = y
        for (int j = size - 1; j >= 0; --j) {
            AInv[j * size + i] = y[j];
            for (int k = j + 1; k < size; ++k) {
                AInv[j * size + i] -= L[k * size + j] * AInv[k * size + i];
            }
            AInv[j * size + i] /= L[j * size + j];
        }
    }

}
