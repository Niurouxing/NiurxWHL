#pragma once


#include <vector>
#include <iostream>
#include <complex>

#define sqrt2 1.4142135623730951
#define sqrt2_div 0.7071067811865475



inline void initializeMatrix(double** Mat, int row, int col) {
    
    Mat = new double*[row];
    for (int i = 0; i < row; i++) {
        Mat[i] = new double[col];
    }
}

inline void MatrixMultiplyVector(double** Mat, double* Vec, int row, int col, double* result) {
    for (int i = 0; i < row; i++) {
        result[i] = 0;
        for (int j = 0; j < col; j++) {
            result[i] += Mat[i][j] * Vec[j];
        }
    }
}

inline void MatrixMultiplyVector(double* Mat, double* Vec, int row, int col, double* result) {
    for (int i = 0; i < row; i++) {
        result[i] = 0;
        for (int j = 0; j < col; j++) {
            result[i] += Mat[i * col + j] * Vec[j];
        }
    }
}

inline void MatrixMultiplyVector(std::complex<double>** Mat, std::complex<double>* Vec, int row, int col, std::complex<double>* result) {
    for (int i = 0; i < row; i++) {
        result[i] = 0;
        for (int j = 0; j < col; j++) {
            result[i] += Mat[i][j] * Vec[j];
        }
    }
}



inline void MatrixTransposeMultiplyVector(double** Mat, double* Vec, int row, int col, double* result) {
    for (int i = 0; i < col; i++) {
        result[i] = 0;
        for (int j = 0; j < row; j++) {
            result[i] += Mat[j][i] * Vec[j];
        }
    }
}

inline void MatrixTransposeMultiplyVector(double* Mat, double* Vec, int row, int col, double* result) {
    for (int i = 0; i < col; i++) {
        result[i] = 0;
        for (int j = 0; j < row; j++) {
            result[i] += Mat[j * col + i] * Vec[j];
        }
    }
}

inline void MatrixTransposeMultiplyVector(std::complex<double>** Mat, std::complex<double>* Vec, int row, int col, std::complex<double>* result) {
    for (int i = 0; i < col; i++) {
        result[i] = 0;
        for (int j = 0; j < row; j++) {
            result[i] += std::conj(Mat[j][i]) * Vec[j];
        }
    }
}

inline void MatrixMultiplyMatrix(double** Mat1, double** Mat2, int row1, int col1, int col2, double** result) {
    for (int i = 0; i < row1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i][j] = 0;
            for (int k = 0; k < col1; k++) {
                result[i][j] += Mat1[i][k] * Mat2[k][j];
            }
        }
    }
}

inline void MatrixMultiplyMatrix(double* Mat1, double* Mat2, int row1, int col1, int col2, double* result) {
    for (int i = 0; i < row1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i * col2 + j] = 0;
            for (int k = 0; k < col1; k++) {
                result[i * col2 + j] += Mat1[i * col1 + k] * Mat2[k * col2 + j];
            }
        }
    }
}

inline void MatrixMultiplyMatrix(double** Mat1, double** Mat2, int row1, int col1, int col2, double* result) {
    for (int i = 0; i < row1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i * col2 + j] = 0;
            for (int k = 0; k < col1; k++) {
                result[i * col2 + j] += Mat1[i][k] * Mat2[k][j];
            }
        }
    }
}

inline void MatrixMultiplyMatrix(std::complex<double>** Mat1, std::complex<double>** Mat2, int row1, int col1, int col2, std::complex<double>** result) {
    for (int i = 0; i < row1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i][j] = std::complex<double>(0, 0);
            for (int k = 0; k < col1; k++) {
                result[i][j] += Mat1[i][k] * Mat2[k][j];
            }
        }
    }
}

inline void MatrixTransposeMultiplyMatrix(double** Mat1, double** Mat2, int row1, int col1, int col2, double** result) {
    for (int i = 0; i < col1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i][j] = 0;
            for (int k = 0; k < row1; k++) {
                result[i][j] += Mat1[k][i] * Mat2[k][j];
            }
        }
    }
}

inline void MatrixTransposeMultiplyMatrix(double* Mat1, double* Mat2, int row1, int col1, int col2, double* result) {
    for (int i = 0; i < col1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i * col2 + j] = 0;
            for (int k = 0; k < row1; k++) {
                result[i * col2 + j] += Mat1[k * col1 + i] * Mat2[k * col2 + j];
            }
        }
    }
}

inline void MatrixTransposeMultiplyMatrix(double** Mat1, double** Mat2, int row1, int col1, int col2, double* result) {
    for (int i = 0; i < col1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i * col2 + j] = 0;
            for (int k = 0; k < row1; k++) {
                result[i * col2 + j] += Mat1[k][i] * Mat2[k][j];
            }
        }
    }
}

inline void MatrixTransposeMultiplyMatrix(std::complex<double>** Mat1, std::complex<double>** Mat2, int row1, int col1, int col2, std::complex<double>** result) {
    for (int i = 0; i < col1; i++) {
        for (int j = 0; j < col2; j++) {
            result[i][j] = std::complex<double>(0, 0);
            for (int k = 0; k < row1; k++) {
                result[i][j] += std::conj(Mat1[k][i]) * Mat2[k][j];
            }
        }
    }
}
    
// only for test, never use it in real project
template <typename T>
void printVector(T* Vec, int length) {
    for (int i = 0; i < length; i++) {
        std::cout << Vec[i] << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void printMatrix(T** Mat, int row, int col) {
    for (int i = 0; i < row; i++) {
        printVector(Mat[i], col);
    }
    std::cout << std::endl;
}

template <typename T>
void printMatrix(T* Mat, int row, int col) {
    for (int i = 0; i < row; i++) {
        printVector(Mat + i * col, col);
    }
    std::cout << std::endl;
}


