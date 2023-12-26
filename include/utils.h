#pragma once

#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>


#include "cblas.h"


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
    for (int i = 0; i < row; i++)
    {
        result[i] = 0;
        for (int j = 0; j < col; j++)
        {
            result[i] += Mat[i * col + j] * Vec[j];
        }
    }
}

inline void MatrixMultiplyVector(std::complex<double> *Mat, std::complex<double> *Vec, int row, int col, std::complex<double> *result)
{
    for (int i = 0; i < row; i++)
    {
        result[i] = 0;
        for (int j = 0; j < col; j++)
        {
            result[i] += Mat[i * col + j] * Vec[j];
        }
    }
}


inline void MatrixTransposeMultiplyVector(double *Mat, double *Vec, int row, int col, double *result)
{
    for (int i = 0; i < col; i++)
    {
        result[i] = 0;
        for (int j = 0; j < row; j++)
        {
            result[i] += Mat[j * col + i] * Vec[j];
        }
    }
}

inline void MatrixTransposeMultiplyVector(std::complex<double> *Mat, std::complex<double> *Vec, int row, int col, std::complex<double> *result)
{
    for (int i = 0; i < col; i++)
    {
        result[i] = 0;
        for (int j = 0; j < row; j++)
        {
            result[i] += std::conj(Mat[j * col + i]) * Vec[j];
        }
    }
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
