#pragma once

#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>


#include "cblas.h"
#include <lapacke.h>


#define sqrt2 1.4142135623730951
#define sqrt2_div 0.7071067811865475

inline void sortIndice(double *arr, int *result, int n)
{
    // Initialize the indices
    for (int i = 0; i < n; i++)
    {
        result[i] = i;
    }

    // Use std::sort with a lambda function as the comparator
    std::sort(result, result + n, [&arr](const int &a, const int &b)
              { return arr[a] > arr[b]; });
}

// 堆基于排序，O(nlogk)，有个快排的O(n)算法，但是太难写
inline void mink(double *arr, int n, int *result, int k)
{
    // 初始化result数组
    std::fill_n(result, k, -1);

    for (int i = 0; i < n; ++i)
    {
        // 查找应插入的位置
        int j = 0;
        while (j < k && result[j] != -1 && arr[i] >= arr[result[j]])
        {
            ++j;
        }

        if (j < k)
        {
            // 向后移动元素以为新元素腾出空间
            for (int m = k - 1; m > j; --m)
            {
                result[m] = result[m - 1];
            }
            result[j] = i; // 插入新索引
        }
    }
}

inline void maxk(double *arr, int n, int *result, int k)
{
    // 初始化result数组
    std::fill_n(result, k, -1);

    for (int i = 0; i < n; ++i)
    {
        // 查找应插入的位置
        int j = 0;
        while (j < k && result[j] != -1 && arr[i] <= arr[result[j]])
        {
            ++j;
        }

        if (j < k)
        {
            // 向后移动元素以为新元素腾出空间
            for (int m = k - 1; m > j; --m)
            {
                result[m] = result[m - 1];
            }
            result[j] = i; // 插入新索引
        }
    }
}


inline void MatrixMultiplyVector(double *Mat, double *Vec, int row, int col, double *result)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                row, col,
                1.0, Mat, col,
                Vec, 1,
                0.0, result, 1);
}

inline void MatrixMultiplyVector(std::complex<double> *Mat, std::complex<double> *Vec, int row, int col, std::complex<double> *result)
{
    static const std::complex<double> alpha = 1.0;
    static const std::complex<double> beta = 0.0;
    cblas_zgemv(CblasRowMajor, CblasNoTrans,
                row, col,
                &alpha, Mat, col,
                Vec, 1,
                &beta, result, 1);
}


inline void MatrixTransposeMultiplyVector(double *Mat, double *Vec, int row, int col, double *result)
{
    cblas_dgemv(CblasRowMajor, CblasTrans,
                row, col,
                1.0, Mat, col,
                Vec, 1,
                0.0, result, 1);
}

inline void MatrixTransposeMultiplyVector(std::complex<double> *Mat, std::complex<double> *Vec, int row, int col, std::complex<double> *result)
{
    static const std::complex<double> alpha = 1.0;
    static const std::complex<double> beta = 0.0;
    cblas_zgemv(CblasRowMajor, CblasConjTrans,
                row, col,
                &alpha, Mat, col,
                Vec, 1,
                &beta, result, 1);
}



//  BLAS版本
inline void MatrixMultiplyMatrix(double *Mat1, double *Mat2, int row1, int col1, int col2, double *result)
{
    double alpha = 1.0;
    double beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                row1, col2, col1,
                alpha, Mat1, col1, Mat2, col2, beta, result, col2);
}


inline void MatrixTransposeMultiplyMatrix(double *Mat1, double *Mat2, int row1, int col1, int col2, double *result)
{
    for (int i = 0; i < col1; i++)
    {
        for (int j = 0; j < col2; j++)
        {
            result[i * col2 + j] = 0;
            for (int k = 0; k < row1; k++)
            {
                result[i * col2 + j] += Mat1[k * col1 + i] * Mat2[k * col2 + j];
            }
        }
    }
}

inline void MatrixTransposeMultiplyMatrix(std::complex<double> *Mat1, std::complex<double> *Mat2, int row1, int col1, int col2, std::complex<double> *result)
{
    for (int i = 0; i < col1; i++)
    {
        for (int j = 0; j < col2; j++)
        {
            result[i * col2 + j] = 0;
            for (int k = 0; k < row1; k++)
            {
                result[i * col2 + j] += std::conj(Mat1[k * col1 + i]) * Mat2[k * col2 + j];
            }
        }
    }
}

inline void MatrixTransposeMultiplySelf(double *Mat, int row, int col, double *result)
{
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                col, col, row,
                1.0, Mat, col,
                Mat, col,
                0.0, result, col);
}


inline void MatrixTransposeMultiplySelf(std::complex<double> *Mat, int row, int col, std::complex<double> *result)
{
    // 注意：BLAS中复数类型用 void* 传递。
    static const std::complex<double> alpha = 1.0;
    static const std::complex<double> beta = 0.0;

    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
                col, col, row,
                &alpha, Mat, col,
                Mat, col,
                &beta, result, col);
}

 // Function to solve A*X = B for a Hermitian positive definite matrix A, stored in row-major order
inline bool solveHermitianPositiveDefiniteSystem(
    double* A,  // Input matrix A, stored in row-major order, modified on output
    double* B,  // Input matrix B, stored in row-major order, solution X on output
    int n,      // The order of matrix A and number of rows in B
    int nrhs =1    // Number of columns in B (number of right-hand sides)
) {
    // Perform Cholesky decomposition A = L*L^H in row-major order
    int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', n, A, n);
    if (info != 0) {
        std::cerr << "Cholesky decomposition failed with info = " << info << ".\n";
        return false;
    }

    // Solve the equation A*X = B for X in row-major order
    info = LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'U', n, nrhs, A, n, B, 1);
    if (info != 0) {
        std::cerr << "Solving linear system failed with info = " << info << ".\n";
        return false;
    }

    return true;
}

// solve Ax = b, A is a Hermitian positive definite matrix like HtH
// Warning: this function will change the input matrix A, and the result will be stored in b
inline bool solveHermitianPositiveDefiniteSystem(
    std::complex<double>* A,  // Input matrix A, stored in row-major order, modified on output
    std::complex<double>* B,  // Input matrix B, stored in row-major order, solution X on output
    int n,       // The order of matrix A and number of rows in B
    int nrhs=1     // Number of columns in B (number of right-hand sides)
) {
    // Perform Cholesky decomposition A = U^H * U
    int info = LAPACKE_zpotrf(LAPACK_ROW_MAJOR, 'U', n, reinterpret_cast<lapack_complex_double*>(A), n);
    if (info != 0) {
        std::cerr << "Cholesky decomposition failed with info = " << info << ".\n";
        return false;
    }

    // Solve the equation A*X = B for X in row-major order
    info = LAPACKE_zpotrs(LAPACK_ROW_MAJOR, 'U', n, nrhs, reinterpret_cast<lapack_complex_double*>(A), n, reinterpret_cast<lapack_complex_double*>(B), 1); // https://github.com/OpenMathLib/OpenBLAS/issues/2188
    if (info != 0) {
        std::cerr << "Solving linear system failed with info = " << info << ".\n";
        return false;
    }

    return true;
}


// solve inv(A), A is a Hermitian positive definite matrix like HtH
// Warning: this function will change the input matrix A, and the result will be stored in A
// Warning: this function will only return the upper triangular part of the matrix A
inline bool solveHermitianPositiveDefiniteInv(
    double* A,  // Input matrix A, stored in row-major order, modified on output
    int n      // The order of matrix A and number of rows in B
) {
    // Perform Cholesky decomposition A = L*L^H in row-major order
    int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', n, A, n);
    // if (info != 0) {
    //     std::cerr << "Cholesky decomposition failed with info = " << info << ".\n";
    //     return false;
    // }
    
    info = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', n,  A, n);
    // if (info != 0) {
    //     std::cerr << "Solving linear system failed with info = " << info << ".\n";
    //     return false;
    // }

    return true;
}

// only for test, never use it in real project
template <typename T>
void printVector(T *Vec, int length, std::string name = "")
{
    if (name != "")
    {
        std::cout << name << ": " << std::endl;
    }
    for (int i = 0; i < length; i++)
    {
        std::cout << Vec[i] << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void printMatrix(T **Mat, int row, int col, std::string name = "")
{
    if (name != "")
    {
        std::cout << name << ": " << std::endl;
    }
    for (int i = 0; i < row; i++)
    {
        printVector(Mat[i], col);
    }
    std::cout << std::endl;
}

template <typename T>
void printMatrix(T *Mat, int row, int col, std::string name = "")
{
    if (name != "")
    {
        std::cout << name << ": " << std::endl;
    }
    for (int i = 0; i < row; i++)
    {
        printVector(Mat + i * col, col);
    }
    std::cout << std::endl;
}

template <typename T>
void printCube(T ***Mat, int row, int col, int height, std::string name = "")
{
    if (name != "")
    {
        std::cout << name << ": " << std::endl;
    }
    for (int i = 0; i < row; i++)
    {
        printMatrix(Mat[i], col, height, name + "[" + std::to_string(i) + "]");
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
