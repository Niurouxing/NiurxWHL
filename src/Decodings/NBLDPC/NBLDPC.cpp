#include "NBLDPC.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "codeData.h"
#include "BinGF.h"
#include <cstring>
#include "NBTools.h"

void NBLDPC::initialize(const int N, const int K, const int GF, const int n_cv, const int n_vc)
{
    code = new NBLDPCCode();
    table = new NBLDPCTable();
    decoder = new NBLDPCDecoder();

    auto it = codeDataMap.find(std::make_tuple(N, K, GF));
    if (it != codeDataMap.end())
    {
        // 找到了对应的CodeData实例
        CodeData &codeData = it->second;
        code->LoadCode(codeData);
    }
    else
    {
        // 未找到对应的CodeData实例
        std::cerr << "Error: No matching CodeData found for N=" << N << ", K=" << K << ", GF=" << GF << std::endl;
        exit(EXIT_FAILURE);
    }

    table->LoadTable(GF, code->logGF);
    decoder->AllocateDecoder(code, n_cv, n_vc);
    GaussianElimination(code, table);
}

std::map<std::tuple<int, int, int>, CodeData> NBLDPC::codeDataMap = {
    {{96, 48, 64}, N96_K48_GF64},
    {{96, 48, 256}, N96_K48_GF256},
    {{128, 64, 256}, N128_K64_GF256},
    {{512, 256, 256}, N512_K256_GF256},
    };

void NBLDPCCode::LoadCode(const CodeData &codeData)
{
    N = codeData.N;
    M = codeData.M;
    GF = codeData.GF;
    logGF = codeData.logGF;
    K = N - M;
    rate = static_cast<double>(K) / N;

    columnDegree = new int[N];
    std::copy(codeData.columnDegree.begin(), codeData.columnDegree.end(), columnDegree);

    // 直接使用codeData中的rowDegree数据
    rowDegree = new int[M];
    std::copy(codeData.rowDegree.begin(), codeData.rowDegree.end(), rowDegree);

    mat = new int *[M];
    matValue = new int *[M];

    for (int m = 0; m < M; ++m)
    {
        mat[m] = new int[rowDegree[m]];
        matValue[m] = new int[rowDegree[m]];
        // copy the data from codeData
        std::copy(codeData.mat[m].begin(), codeData.mat[m].end(), mat[m]);
        std::copy(codeData.matValue[m].begin(), codeData.matValue[m].end(), matValue[m]);
    }

    nbBranch = 0;
    for (int m = 0; m < M; ++m)
    {
        nbBranch += rowDegree[m];
    }

    // compute the max row degree
    maxRowDegree = 0;
    for (int m = 0; m < M; ++m)
    {
        if (rowDegree[m] > maxRowDegree)
        {
            maxRowDegree = rowDegree[m];
        }
    }
}

NBLDPCCode::~NBLDPCCode()
{
    delete[] columnDegree;
    delete[] rowDegree;
    if (mat != nullptr)
    {
        for (int i = 0; i < M; ++i)
        {
            delete[] mat[i];
        }
        delete[] mat;
    }
    if (matValue != nullptr)
    {
        for (int i = 0; i < M; ++i)
        {
            delete[] matValue[i];
        }
        delete[] matValue;
    }
}

void NBLDPCTable::LoadTable(int GF, int logGF)
{
    int nbRow, nbCol, g, k, l;
    int i, j;

    if (GF != 16 && GF != 32 && GF != 64 && GF != 256 && GF != 4096)
    {
        std::cout << "The binary image of GF(" << GF << ") is not available in this version of the program. Please try GF(64) or GF(256)" << std::endl;
        exit(EXIT_FAILURE);
    }

    nbRow = GF;
    nbCol = logGF;
    // Allocate and initialize BINGF
    BINGF = new int *[nbRow];
    BINGF[0] = new int[nbRow * nbCol];
    for (k = 1; k < nbRow; k++)
    {
        BINGF[k] = BINGF[0] + k * nbCol;
    }

    // Allocate and initialize ADDGF
    ADDGF = new int *[GF];
    ADDGF[0] = new int[GF * GF];
    for (k = 1; k < GF; k++)
    {
        ADDGF[k] = ADDGF[0] + k * GF;
    }

    // Allocate DECGF and GFDEC
    DECGF = new int[GF];
    GFDEC = new int[GF];

    // Allocate and initialize MULGF
    MULGF = new int *[GF];
    MULGF[0] = new int[GF * GF];
    for (k = 1; k < GF; k++)
    {
        MULGF[k] = MULGF[0] + k * GF;
    }

    // Allocate and initialize DIVGF
    DIVGF = new int *[GF];
    DIVGF[0] = new int[GF * GF];
    for (k = 1; k < GF; k++)
    {
        DIVGF[k] = DIVGF[0] + k * GF;
    }

    // Allocate and initialize MULDEC
    MULDEC = new int *[GF];
    MULDEC[0] = new int[GF * GF];
    for (k = 1; k < GF; k++)
    {
        MULDEC[k] = MULDEC[0] + k * GF;
    }

    // Allocate and initialize DIVDEC
    DIVDEC = new int *[GF];
    DIVDEC[0] = new int[GF * GF];
    for (k = 1; k < GF; k++)
    {
        DIVDEC[k] = DIVDEC[0] + k * GF;
    }

    switch (GF)
    {
    case 16:
        memcpy(BINGF[0], BinGF_16, nbRow * nbCol * sizeof(int));
        break;
    case 32:
        memcpy(BINGF[0], BinGF_32, nbRow * nbCol * sizeof(int));
        break;
    case 64:
        memcpy(BINGF[0], BinGF_64, nbRow * nbCol * sizeof(int));
        break;
    case 256:
        memcpy(BINGF[0], BinGF_256, nbRow * nbCol * sizeof(int));
        break;
    // case 4096:
    //     memcpy(BINGF[0], BinGF_4096, nbRow * nbCol * sizeof(int));
    //     break;
    default:
        std::cout << "The binary image of GF(" << GF << ") is not available in this version of the program. Please try GF(64) or GF(256)" << std::endl;
        exit(EXIT_FAILURE);
    }

    // bin2dec
    int sum;
    int tmp;
    for (j = 0; j < GF; j++)
    {
        sum = 0;
        for (i = 0; i < logGF; i++)
        {
            tmp = BINGF[j][i];
            // printf("%d",tmp);
            sum = sum + (tmp << i);
        }
        DECGF[j] = sum;
        GFDEC[sum] = j;
        // printf(" \n bin2dec of GF %d is %d \n",j,sum);
    }
    // getchar();

    // compute the multiplication table in GF(q)
    for (i = 0; i < GF; i++)
    {
        for (j = 0; j < GF; j++)
        {
            if (i == 0 || j == 0)
                MULGF[i][j] = 0;
            else if (i == 1)
                MULGF[i][j] = j;
            else if (j == 1)
                MULGF[i][j] = i;
            else
            {
                tmp = i + j - 2;
                if (tmp < GF - 1)
                    MULGF[i][j] = tmp + 1;
                else
                    MULGF[i][j] = (tmp % (GF - 1)) + 1;
            }
        }
    }

    // compute the division table in GF(q)
    int nb = GF - 1;
    for (i = 0; i < GF; i++)
    {
        for (j = 0; j < GF; j++)
        {
            if (j == 0)
            {
                DIVGF[i][j] = 0;
            }
            else if (i == 0)
            {
                DIVGF[i][j] = 0;
            }
            else if (j == 1)
            {
                DIVGF[i][j] = i;
            }
            else
            {
                DIVGF[i][j] = nb--;
            };
            if (nb < 1)
            {
                nb = GF - 1;
            }
        }
    }

    // compute the multiplication table in decimal
    for (i = 0; i < GF; i++)
    {
        for (j = 0; j < GF; j++)
        {
            MULDEC[i][j] = DECGF[MULGF[i][j]];
        }
    }

    // compute the division table in decimal
    for (i = 0; i < GF; i++)
    {
        for (j = 0; j < GF; j++)
        {
            DIVDEC[DECGF[i]][j] = DIVGF[i][j];
        }
    }
}

NBLDPCTable::~NBLDPCTable()
{
    if (BINGF != nullptr)
    {
        delete[] BINGF[0];
        delete[] BINGF;
    }
    if (ADDGF != nullptr)
    {
        delete[] ADDGF[0];
        delete[] ADDGF;
    }
    if (DECGF != nullptr)
    {
        delete[] DECGF;
    }
    if (GFDEC != nullptr)
    {
        delete[] GFDEC;
    }
    if (MULGF != nullptr)
    {
        delete[] MULGF[0];
        delete[] MULGF;
    }
    if (DIVGF != nullptr)
    {
        delete[] DIVGF[0];
        delete[] DIVGF;
    }
    if (MULDEC != nullptr)
    {
        delete[] MULDEC[0];
        delete[] MULDEC;
    }
    if (DIVDEC != nullptr)
    {
        delete[] DIVDEC[0];
        delete[] DIVDEC;
    }
}

void NBLDPCDecoder::AllocateDecoder(NBLDPCCode *code, int n_cv, int n_vc)
{
    const int N = code->N;
    int nbRow, nbCol, k, nbRow_arr;

    this->nbBranch = code->nbBranch;
    this->N = code->N;

    this->n_cv = n_cv;
    this->n_vc = n_vc;

    nbRow = code->nbBranch;
    nbCol = this->n_vc;
    nbRow_arr = code->rowDegree[0];

    /* VtoC [nbBranch][nbMax] */
    VtoC = new double *[nbRow];
    VtoC[0] = new double[nbRow * code->GF];
    for (k = 1; k < nbRow; k++)
        VtoC[k] = VtoC[0] + k * code->GF;

    /* M_CtoV_LLR [nbBranch][nbMax] */
    M_CtoV_LLR = new double *[nbRow_arr];
    M_CtoV_LLR[0] = new double[nbRow_arr * code->GF];
    for (k = 1; k < nbRow_arr; k++)
        M_CtoV_LLR[k] = M_CtoV_LLR[0] + k * code->GF;

    /* M_VtoC_LLR [nbBranch][nbMax] */
    M_VtoC_LLR = new double *[nbRow_arr];
    M_VtoC_LLR[0] = new double[nbRow_arr * nbCol];
    for (k = 1; k < nbRow_arr; k++)
        M_VtoC_LLR[k] = M_VtoC_LLR[0] + k * nbCol;

    /* M_CtoV_GF [rowDegree[0]][nbMax] */
    M_CtoV_GF = new int *[nbRow_arr];
    M_CtoV_GF[0] = new int[nbRow_arr * code->GF];
    for (k = 1; k < nbRow_arr; k++)
        M_CtoV_GF[k] = M_CtoV_GF[0] + k * code->GF;

    /* M_VtoC_GF [rowDegree[0]][nbMax] */
    M_VtoC_GF = new int *[nbRow_arr];
    M_VtoC_GF[0] = new int[nbRow_arr * nbCol];
    for (k = 1; k < nbRow_arr; k++)
        M_VtoC_GF[k] = M_VtoC_GF[0] + k * nbCol;

    nbRow = N;
    nbCol = code->GF;
    /* APP [N][GF] */
    APP = new double *[nbRow];
    APP[0] = new double[nbRow * nbCol];
    for (k = 1; k < nbRow; k++)
        APP[k] = APP[0] + k * nbCol;

    nbRow = N;
    nbCol = code->GF;
    /* intrinsic_LLR [N][GF] */ /* VN modified */
    intrinsic_LLR = new double *[nbRow];
    intrinsic_LLR[0] = new double[nbRow * nbCol];
    for (k = 1; k < nbRow; k++)
        intrinsic_LLR[k] = intrinsic_LLR[0] + k * nbCol;

    /* intrinsic_GF [N][GF] */ /* VN modified */
    intrinsic_GF = new int *[nbRow];
    intrinsic_GF[0] = new int[nbRow * nbCol];
    for (k = 1; k < nbRow; k++)
        intrinsic_GF[k] = intrinsic_GF[0] + k * nbCol;
}

NBLDPCDecoder::~NBLDPCDecoder()
{
    if (CtoV != nullptr)
    {
        delete[] CtoV[0];
        delete[] CtoV;
    }
    if (VtoC != nullptr)
    {
        delete[] VtoC[0];
        delete[] VtoC;
    }
    if (intrinsic_LLR != nullptr)
    {
        delete[] intrinsic_LLR[0];
        delete[] intrinsic_LLR;
    }
    if (intrinsic_GF != nullptr)
    {
        delete[] intrinsic_GF[0];
        delete[] intrinsic_GF;
    }
    if (APP != nullptr)
    {
        delete[] APP[0];
        delete[] APP;
    }
    if (M_VtoC_LLR != nullptr)
    {
        delete[] M_VtoC_LLR[0];
        delete[] M_VtoC_LLR;
    }
    if (M_VtoC_GF != nullptr)
    {
        delete[] M_VtoC_GF[0];
        delete[] M_VtoC_GF;
    }
    if (M_CtoV_LLR != nullptr)
    {
        delete[] M_CtoV_LLR[0];
        delete[] M_CtoV_LLR;
    }
    if (M_CtoV_GF != nullptr)
    {
        delete[] M_CtoV_GF[0];
        delete[] M_CtoV_GF;
    }
}