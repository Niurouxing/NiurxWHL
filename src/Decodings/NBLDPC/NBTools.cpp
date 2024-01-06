
#include "NBTools.h"
#include "NBLDPC.h"
#include <iostream>
#include <fstream>
#include "BinGF.h"
#include <random>
#include <vector>
#include <chrono>
#include <cstring>
#include <stdio.h>
#include <limits.h> // 对于INT_MAX
#include <string.h> // 对于memset
#include <algorithm> // 对于std::partial_sort
 

/**
 * \fn GaussianElimination
 * \brief Perform a Gaussian elimination on the parity-check matrix.
 * 		   The procedure stops when the
 * 		   The redundancy symbols are initialized to 0.
 * Inputs
 * 	- code->mat   : parity-check matrix
 * 	- table       : lookup tables for computations in GF(q)
 * Outputs       :
 *      - code->matUT : Upper triangular matrix used for encoding
 *      - code->Perm : Column permutation
 */
void GaussianElimination(NBLDPCCode *code, NBLDPCTable *table) {
    const int N = code->N;
    const int M = code->M;

    code->matUT = new int *[M];
    code->matUT[0] = new int[M * N];
    for (int k = 1; k < M; k++) {
        code->matUT[k] = code->matUT[0] + k * N;
    }

    code->Perm = new int[N];
    std::iota(code->Perm, code->Perm + N, 0); // Replace loop for initializing Perm

    for (int m = 0; m < M; m++) {
        for (int k = 0; k < code->rowDegree[m]; k++) {
            code->matUT[m][code->mat[m][k]] = code->matValue[m][k];
        }
    }

    for (int m = 0; m < M; m++) {
        int ind = m;
        while (ind < N && code->matUT[m][ind] == 0) {
            ++ind;
        }

        if (ind == N) {
            std::cout << "The matrix is not full rank (" << m << "," << ind << ")" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::swap(code->Perm[ind], code->Perm[m]);
        for (int m1 = 0; m1 < M; m1++) {
            std::swap(code->matUT[m1][m], code->matUT[m1][ind]);
        }

        for (int m1 = m + 1; m1 < M; m1++) {
            if (code->matUT[m1][m] != 0) {
                int buf = code->matUT[m1][m];
                for (int n = m; n < N; n++) {
                    if (code->matUT[m1][n] != 0) {
                        code->matUT[m1][n] = table->DIVGF[code->matUT[m1][n]][buf];
                    }
                }
                for (int n = m; n < N; n++) {
                    if (code->matUT[m1][n] != 0) {
                        code->matUT[m1][n] = table->MULGF[code->matUT[m1][n]][code->matUT[m][m]];
                    }
                }
                for (int n = m; n < N; n++) {
                    int temp[12];
                    for (int i = 0; i < code->logGF; i++) {
                        temp[i] = table->BINGF[code->matUT[m1][n]][i] ^ table->BINGF[code->matUT[m][n]][i];
                    }
                    code->matUT[m1][n] = Bin2GF(temp, code->GF, code->logGF, table->BINGF);
                }
            }
        }
    }
}

/*!
 * \fn Bin2GF
 * \brief compute a GF(q) symbol corresponding to a frame of log_2(GF) bits
 * Parameters    :
 * Inputs        :
 * 	- int *U    : array representing logGF bits
 * 	- int logGF : size of the array U. logGF = log2 (GF)
 * 	- int GF    : order of the field
 * 	- int ** BINGF: binary mapping table
 * Outputs       :
 *      - index of the non-binary symbol
 */

int Bin2GF(int *U, int GF, int logGF, int **BINGF)
{
    int sum = 0;
    for (int i = 0; i < logGF; i++)
    {
        sum = sum + (U[i] << (logGF - i - 1));
    }

    switch (GF)
    {
    case 16:
        return BinGF_16_aux[sum];
        break;
    case 32:
        return BinGF_32_aux[sum];
        break;
    case 64:
        return BinGF_64_aux[sum];
        break;
    case 256:
        return BinGF_256_aux[sum];
        break;
    // case 4096:
    //     return BinGF_4096_aux[sum];
    //     break;
    default:
        std::cout << "The binary image of GF(" << GF << ") is not available in this version of the program. Please try GF(64) or GF(256)" << std::endl;
        exit(EXIT_FAILURE);
    }
}

 

/**
 * \fn Encoding
 * \brief Encode the information bits into a codeword.
 * 		   matUT beeing upper triangular, the backsubstitution method is used.
 * 		   The M first symbols in NSYMB are redundancy symbols (before deinterleaving)
 * Inputs
 * 	- KSYMB  ( KSYMB are information symbols)
 * Outputs
 *      - Codeword
 *      - NBIN : binary copy of the codeword
 */

void Encoding(NBLDPCCode *code, NBLDPCTable *table, int *CodeWord, int **NBIN, int *KSYMB)
{
    const int N = code->N;
    const int M = code->M;
    const int logGF = code->logGF;
    int k, n, m, q, buf;
    int NSYMB[N];
    int temp[logGF];
    int mulGF_val, binGF_val;

    // 初始化符号
    memcpy(&NSYMB[M], KSYMB, sizeof(int) * (N - M));

    // 回代
    for (m = M - 1; m >= 0; m--)
    {
        buf = 0;
        for (n = m + 1; n < N; n++)
        {
            int matUT_val = code->matUT[m][n];
            if (matUT_val != 0)
            {
                mulGF_val = table->MULGF[matUT_val][NSYMB[n]];
                for (int i = 0; i < logGF; i++)
                {
                    binGF_val = table->BINGF[mulGF_val][i];
                    temp[i] = (table->BINGF[buf][i]) ^ binGF_val;
                }
                buf = Bin2GF(temp, code->GF, logGF);
            }
        }
        NSYMB[m] = table->DIVGF[buf][code->matUT[m][m]];
    }

    // 解交织
    for (n = 0; n < N; n++)
    {
        CodeWord[code->Perm[n]] = NSYMB[n];
    }

    // 二进制化码字
    for (n = 0; n < N; n++)
    {
        int *binGF_CodeWord_n = table->BINGF[CodeWord[n]];
        memcpy(NBIN[n], binGF_CodeWord_n, sizeof(int) * logGF);
    }
}
void EMS(NBLDPCCode *code, NBLDPCTable *table, NBLDPCDecoder *decoder, int NbOper, double offset, int NbIterMax, int *decide)
{
    double Mvc_temp[16][256];
    int synd = 0;
    int node = 0;
    int k, l, n, iter, i, g;

    // init Mvc with intrinsic
    for (int n = 0; n < code->N; n++) {
        std::copy(decoder->intrinsic_LLR[n], decoder->intrinsic_LLR[n] + code->GF, decoder->VtoC[n]);
    }

    // Decoding iterations
    for (iter = 0; iter < NbIterMax; iter++) // NbIterMax=5
    {

        for (node = 0; node < code->M; node++) // Loop for the M Check nodes
        {
            for (i = 0; i < code->rowDegree[node]; i++) // rowDegree[i] = the i^th check node degree 对每一个v_i
            {
                memcpy(Mvc_temp[i], decoder->VtoC[code->mat[node][i]], code->GF * sizeof(double));
            }

            // sorting Mvc values
            for (i = 0; i < code->rowDegree[node]; i++)
            {
                // Use nth_element to partition the first n_vc elements
                std::nth_element(&Mvc_temp[i][0], &Mvc_temp[i][decoder->n_vc], &Mvc_temp[i][code->GF]);

                std::sort(&Mvc_temp[i][0], &Mvc_temp[i][decoder->n_vc]);

                for (k = 0; k < decoder->n_vc; k++)
                {
                    decoder->M_VtoC_LLR[i][k] = Mvc_temp[i][k];
                    // Find the index of the element in the original array
                    decoder->M_VtoC_GF[i][k] = std::find(&decoder->VtoC[code->mat[node][i]][0], &decoder->VtoC[code->mat[node][i]][code->GF], Mvc_temp[i][k]) - &decoder->VtoC[code->mat[node][i]][0];
                }

                // Normalisation
                for (g = 1; g < decoder->n_vc; g++)
                {
                    decoder->M_VtoC_LLR[i][g] = decoder->M_VtoC_LLR[i][g] - decoder->M_VtoC_LLR[i][0];
                }
                decoder->M_VtoC_LLR[i][0] = 0.0;
            }

            CheckPassLogEMS(node, decoder, code, table, NbOper, offset);
            // CheckPassLogEMS_dc3(node, &decoder, &code, &table, NbOper, offset);

            // compute SO
            for (i = 0; i < code->rowDegree[node]; i++)
            {
                for (k = 0; k < code->GF; k++)
                {
                    decoder->APP[code->mat[node][i]][k] = decoder->M_CtoV_LLR[i][k] + decoder->VtoC[code->mat[node][i]][k];
                    decoder->VtoC[code->mat[node][i]][k] = decoder->M_CtoV_LLR[i][k] + decoder->intrinsic_LLR[code->mat[node][i]][k];
                }
            }

        } // End of the node update
        // Decision(decide, decoder->APP, code->N, code->GF);
        // synd = Syndrom(code, decide, table);
        // if (synd == 0)
        //     break;
    }
    Decision(decide, decoder->APP, code->N, code->GF);
    return;
}







void CheckPassLogEMS(int node, NBLDPCDecoder *decoder, NBLDPCCode *code, NBLDPCTable *table, int NbOper, double offset)
{
    const int nbMax = decoder->n_vc;
    const int rowDegreeNode = code->rowDegree[node];
    const int S = 2 * (rowDegreeNode - 2);
    int Stp = 0;

    double **MatriceInter = new double*[S];
    int **MatriceInterIndice = new int*[S];
    double *OutForward = new double[nbMax];
    double *OutBackward = new double[nbMax];
    double *OutForward1 = new double[nbMax];
    double *OutBackward1 = new double[nbMax];
    int *OutForwardIndice = new int[nbMax];
    int *OutBackwardIndice = new int[nbMax];
    int *OutForwardIndice1 = new int[nbMax];
    int *OutBackwardIndice1 = new int[nbMax];
    double LLR_tmp[code->GF];

    // Initialization
    std::fill_n(OutForward, nbMax, 1e5);
    std::fill_n(OutBackward, nbMax, 1e5);
    std::fill_n(OutForwardIndice, nbMax, -1);
    std::fill_n(OutBackwardIndice, nbMax, -1);
    
    for (int i = 0; i < S; ++i) {
        MatriceInter[i] = new double[nbMax];
        std::fill_n(MatriceInter[i], nbMax, 1e5);
        MatriceInterIndice[i] = new int[nbMax];
        std::fill_n(MatriceInterIndice[i], nbMax, -1);
    }

    // The following nested loops are removed and replaced with std::fill_n above

    // Rotate non-zero values in matrix for VtoC
    for (int t = 0; t < rowDegreeNode; ++t) {
        for (int k = 0; k < nbMax; ++k) {
            int c = decoder->M_VtoC_GF[t][k];
            if (c != -1) {
                decoder->M_VtoC_GF[t][k] = table->MULDEC[c][code->matValue[node][t]];
            }
        }
    }

    // Initialization of LLR values from M_VtoC_LLR
    for (int k = 0; k < nbMax; k++) {
        OutForward[k] = decoder->M_VtoC_LLR[0][k];
        OutBackward[k] = decoder->M_VtoC_LLR[rowDegreeNode - 1][k];
        OutForwardIndice[k] = decoder->M_VtoC_GF[0][k];
        OutBackwardIndice[k] = decoder->M_VtoC_GF[rowDegreeNode - 1][k];
    }

    // Main processing loop with ElementaryStep calls
    // Omitted for brevity - no changes here


    // Start of the recusivity开始递归
    for (int kk = 1; kk < (code->rowDegree[node] - 1); kk++) // the middle dc-2 ones
    {

        // for VN No.kk
        // storage M_VtoC_LLR[kk][] to Outforward1(tmp)
        // storage M_VtoC_LLR[dc-kk-1][] to Backforward1(tmp)
        for (int k = 0; k < nbMax; k++)
        {
            // forward step
            OutForward1[k] = decoder->M_VtoC_LLR[kk][k];
            OutForwardIndice1[k] = decoder->M_VtoC_GF[kk][k];
            // backward step
            OutBackward1[k] = decoder->M_VtoC_LLR[code->rowDegree[node] - kk - 1][k];
            OutBackwardIndice1[k] = decoder->M_VtoC_GF[code->rowDegree[node] - kk - 1][k];
        }

        // MatriceInter
        for (int k = 0; k < nbMax; k++)
        {
            // Storage of the intermediate vectors
            MatriceInter[kk - 1][k] = OutForward[k];
            MatriceInterIndice[kk - 1][k] = OutForwardIndice[k];
            MatriceInter[2 * (code->rowDegree[node] - 2) - kk][k] = OutBackward[k];
            MatriceInterIndice[2 * (code->rowDegree[node] - 2) - kk][k] = OutBackwardIndice[k];
        }

        // forward step
        // printf(" \n forward step %d \n",kk);
        if (kk < (code->rowDegree[node] - 1))
            ElementaryStep(OutForward, OutForward1, OutForwardIndice, OutForwardIndice1, OutForward, OutForwardIndice, table->ADDGF, code->GF, nbMax, NbOper);

        // backward step
        // printf(" \n backward step %d \n",kk);
        if (kk < (code->rowDegree[node] - 1))
            ElementaryStep(OutBackward, OutBackward1, OutBackwardIndice, OutBackwardIndice1, OutBackward, OutBackwardIndice, table->ADDGF, code->GF, nbMax, NbOper);
    }

    // Update of vectors CtoV (first and last)
    for (int k = 0; k < nbMax; k++)
    {
        // last vector M_CtoV_LLR
        decoder->M_CtoV_LLR[code->rowDegree[node] - 1][k] = OutForward[k];
        decoder->M_CtoV_GF[code->rowDegree[node] - 1][k] = OutForwardIndice[k];
        // first vector M_CtoV_LLR
        decoder->M_CtoV_LLR[0][k] = OutBackward[k];
        decoder->M_CtoV_GF[0][k] = OutBackwardIndice[k];
    }

    // printf("\n merge \n");
    // Update of the others vectors CtoV
    for (int k = 0; k < (code->rowDegree[node] - 2); k++)
    {

        ElementaryStep(MatriceInter[k], MatriceInter[(code->rowDegree[node] - 2) + k], MatriceInterIndice[k], MatriceInterIndice[(code->rowDegree[node] - 2) + k], OutForward, OutForwardIndice, table->ADDGF, code->GF, nbMax, NbOper);
        for (int g = 0; g < nbMax; g++)
        {
            decoder->M_CtoV_LLR[k + 1][g] = OutForward[g];
            decoder->M_CtoV_GF[k + 1][g] = OutForwardIndice[g];
        }
    }

    // Rotation par la valeur non nulle dans la matrice pour CtoV 矩阵中非零值的旋转
    for (int t = 0; t < code->rowDegree[node]; t++)
    {
        for (int k = 0; k < nbMax; k++)
        {

            if (decoder->M_CtoV_GF[t][k] == -1)
            {
                Stp = k;
                break;
            }
            else
                Stp = nbMax;
        }

        // back to GF symbol and devise
        for (int k = 0; k < Stp; k++)
        {
            // decoder->M_CtoV_GF[t][k]=table->DIVGF[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
            decoder->M_CtoV_GF[t][k] = table->DIVDEC[decoder->M_CtoV_GF[t][k]][code->matValue[node][t]];
        }

        /// reorder in GF order and add offset
        // printf("sat=%0.3f \n", decoder->M_CtoV_LLR[t][Stp-1] + offset );
        // getchar();

        for (int k = 0; k < code->GF; k++)
        {
            LLR_tmp[k] = decoder->M_CtoV_LLR[t][Stp - 1] + offset;
        }

        for (int k = 0; k < Stp; k++)
        {
            LLR_tmp[decoder->M_CtoV_GF[t][k]] = decoder->M_CtoV_LLR[t][k];
        }

        for (int k = 0; k < code->GF; k++)
        {
            decoder->M_CtoV_GF[t][k] = k;
            decoder->M_CtoV_LLR[t][k] = LLR_tmp[k];
        }
    }
    for (int i = 0; i < S; ++i) {
        delete[] MatriceInter[i];
        delete[] MatriceInterIndice[i];
    }
    delete[] MatriceInter;
    delete[] MatriceInterIndice;
    delete[] OutForward;
    delete[] OutBackward;
    delete[] OutForward1;
    delete[] OutBackward1;
    delete[] OutForwardIndice;
    delete[] OutBackwardIndice;
    delete[] OutForwardIndice1;
    delete[] OutBackwardIndice1;
}


void ElementaryStep(double *Input1, double *Input2, int *IndiceInput1, int *IndiceInput2, double *Output, int *IndiceOut, 
                    int **ADDGF, int GF, int nbMax, int nbOper)
{
    // 由于数组始终是256大小，可以在函数外定义为宏以避免硬编码
    #define MAX_SIZE 256

    double loc_Output[MAX_SIZE];
    int loc_IndiceOut[MAX_SIZE];

    int i, j, s;
    int nb_bubble = 8;
    int pos;

    double tab_aux[MAX_SIZE][MAX_SIZE];
    memset(tab_aux, 0, sizeof(tab_aux)); // 只需一次memset来代替两个循环

    int GFvalues_already_in_Out[MAX_SIZE];
    memset(GFvalues_already_in_Out, -1, sizeof(GFvalues_already_in_Out)); // 使用memset替代循环

    // 初始条件现在仅通过一次赋值完成
    for (i = 0; i < nbMax; i++)
    {
        loc_Output[i] = 1e5;
        loc_IndiceOut[i] = -1;
    }

    // pre-compute the addition of bubble check
    int half_nb_bubble = nb_bubble >> 1;
    for (i = 0; i < nbMax; i++)
    {
        for (j = 0; j < half_nb_bubble; j++)
        {
            tab_aux[i][j] = Input1[i] + Input2[j];
        }
        for (j = half_nb_bubble; j < nbMax; j++)
        {
            tab_aux[i][j] = Input1[i] + Input2[j];
        }
    }

    double tab_comp[3][MAX_SIZE]; // Assuming that 100 can be replaced with MAX_SIZE safely

    // 用结构体替换表，便于排序和管理
    struct {
        double value;
        int idx_i;
        int idx_j;
    } competitors[nb_bubble];

    // 初始化竞争者表
    for (j = 0; j < half_nb_bubble; j++)
    {
        competitors[j].value = tab_aux[j][0];
        competitors[j].idx_i = j;
        competitors[j].idx_j = 0;
    }
    for (j = 0; j < half_nb_bubble; j++)
    {
        competitors[j + half_nb_bubble].value = tab_aux[half_nb_bubble][j];
        competitors[j + half_nb_bubble].idx_i = half_nb_bubble;
        competitors[j + half_nb_bubble].idx_j = j;
    }

    // 实现一个简单的选择排序，这样可以在后续操作中减少minimum函数的调用次数
    for (s = 0; s < nbOper && s < nbMax; s++)
    {
        pos = 0;
        for (j = 1; j < nb_bubble; j++)
        {
            if (competitors[j].value < competitors[pos].value)
            {
                pos = j;
            }
        }

        int indice_aux = IndiceInput1[competitors[pos].idx_i] ^ IndiceInput2[competitors[pos].idx_j];
        if (GFvalues_already_in_Out[indice_aux] == -1)
        {
            loc_Output[s] = competitors[pos].value;
            loc_IndiceOut[s] = indice_aux;
            GFvalues_already_in_Out[indice_aux] = s; // 存储位置信息，而不仅仅是标记
        }

        // 更新竞争者表
        if (competitors[pos].idx_i < nbMax - 1 && competitors[pos].idx_j < nbMax - 1)
        {
            if (pos >= half_nb_bubble)
            {
                competitors[pos].idx_i++;
            }
            else
            {
                competitors[pos].idx_j++;
            }
            competitors[pos].value = tab_aux[competitors[pos].idx_i][competitors[pos].idx_j];
        }
        else
        {
            // 如果到达边界，则用极大值替换，确保不会被再次选择
            competitors[pos].value = 1e5;
        }
    }

    memcpy(Output, loc_Output, sizeof(double) * nbMax); // 使用memcpy替换循环
    memcpy(IndiceOut, loc_IndiceOut, sizeof(int) * nbMax); // 使用memcpy替换循环

    #undef MAX_SIZE
}

void Decision(int *decision, double **APP, int N, int GF)
{
    for (int n = 0; n < N; n++)
    {
        auto minIt = std::min_element(APP[n], APP[n] + GF);
 
        decision[n] = std::distance(APP[n], minIt);
    }
}

int Syndrom(NBLDPCCode *code, int *decide, NBLDPCTable *tableGF)
{
    int k, l, i;
    int synd = 0;
    int temp[code->logGF]; // 假设code->logGF最大为12，否则需要动态分配
    
    for (k = 0; k < code->M && synd == 0; k++) // 避免synd != 0后的无效循环
    {
        for (l = 0; l < code->rowDegree[k]; l++)
        {
            int *currentBINGF = tableGF->BINGF[synd]; // 提前获取索引为synd的BINGF
            int *mulGF = tableGF->MULGF[code->matValue[k][l]];
            int *decidedBINGF = tableGF->BINGF[mulGF[decide[code->mat[k][l]]]];

            for (i = 0; i < code->logGF; i++)
            {
                temp[i] = currentBINGF[i] ^ decidedBINGF[i];
            }
            synd = Bin2GF(temp, code->GF, code->logGF, tableGF->BINGF); // 评估Bin2GF是否可能优化
        }
    }

    return synd;
}