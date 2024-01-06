#pragma once
#include <map>

#include "baseCode.h"
struct CodeData;

class NBLDPCCode
{
public:
    int N;             /* number of columns in H */
    int M;             /* number of rows in H */
    int K;             /* number of information symbols : K = N-M */
    int GF;            /* Field order (eg. GF=256) */
    int logGF;         /* logGF = log2(GF)  ( given GF = 2^q => logGF = q) : logGF is the number of bits forming GF symbols */
    int **mat;         /* the Parity-check matrix : tableau bidimensionnel contenant les VN (i.e les colonnes) qui participe dans chaque contrainte
                       de parit�  (i.e les lignes) */
    int **matValue;    /* Parity-check matrix non-binary coefficients : contient les coefficient GF(q) pour chaque ligne */
    int *rowDegree;    /* rowDegree[i] = the i^th check node degree */
    int *columnDegree; /* columnDegree[j] = the j^th variable node degree */
    int nbBranch;      /* number of edges in the Tanner graph */
    double rate;       /* Code rate */
    int **matUT;       /* Upper Triangular form of the parity-check matrix after Gaussian elimination. matUT is used for encoding*/
    int *Perm;         /* Permutation used for Gaussian Elimination  */

    int maxRowDegree; /* maximum row degree */

    void LoadCode(const CodeData &codeData);
    ~NBLDPCCode();
};

class NBLDPCTable
{
public:
    int **BINGF;  /* Mapping symbol GFq -> ensemble de symboles binaires */
    int **ADDGF;  /* Addition table in GFq */
    int **MULGF;  /* Multiplication table in GFq */
    int **DIVGF;  /* Division table in GFq */
    int *DECGF;   /*Mapping symbol GFq -> binary converted to decimal */
    int *GFDEC;   /*Mapping decimal to GF symbol */
    int **MULDEC; /* Multiplication in decimal */
    int **DIVDEC; /* Division in decimal */

    void LoadTable(int GF, int logGF);
    ~NBLDPCTable();
};

class NBLDPCDecoder
{
public:
    int N;
    int n_cv;               /* only the n_cv most reliable GF are transmitted from Check Node to Variable Node  */
    int n_vc;               /* only the n_vc most reliable GF are transmitted from Check Node to Variable Node  */
    int nbBranch;           /* nombre de branches dans le graphe de Tanner */
    double **CtoV;          /* An array to store the nbMax largest Check_to_Variable reliabilities for all edges on the graph CtoV Array size = (nbBranch x nbMax) */
    double **VtoC;          /* An array to store the nbMax largest Variable to Check reliabilities for all edges on the graph VtoC Array size = (nbBranch x nbMax) */
    double **intrinsic_LLR; /* An array to store intrinsic Log intrinsic_LLR Ratios received from channel. Each variable node has GF intrinsic LLRs. So, the size of intrinsic_LLR is (N x GF) where N is the number of VNs and GF is the order of the Galois field. The values are sorted in decreasing way */
    int **intrinsic_GF;     /*  Galois field symbols associated to the LLRs values in 'intrinsic_LLR' Same size as for 'intrinsic_LLR' */
    double **APP;           /* Array to store A Posteriori Probability used to make hard decision on symbol value Same size as for 'intrinsic_LLR'*/
    double **M_VtoC_LLR;    // LLR inputs to one Check node processor.
    int **M_VtoC_GF;        // GF index corresponding to M_VtoC_LLR.
    double **M_CtoV_LLR;    // LLR output from one Check node processor.
    int **M_CtoV_GF;        // GF index corresponding to M_VtoC_LLR.

    // used in CheckPassLogEMS()
    double **MatriceInter;
    int **MatriceInterIndice;
    double *OutForward;
    double *OutBackward;
    double *OutForward1;
    double *OutBackward1;
    int *OutForwardIndice;
    int *OutBackwardIndice;
    int *OutForwardIndice1;
    int *OutBackwardIndice1;
    double *LLR_tmp;

    void AllocateDecoder(NBLDPCCode *code, int n_cv = 20, int n_vc = 20);
    ~NBLDPCDecoder();
};

class NBLDPC : public BaseCode
{
public:
    NBLDPCCode *code;
    NBLDPCTable *table;
    NBLDPCDecoder *decoder;

    int **NBIN;
    int **KBIN;
    // int * NSYM;
    int *KSYM;

    int *codeWords;
    int *decide;

    NBLDPC(const int N = 96, const int K = 48, const int GF = 256, const int n_cv = 20, const int n_vc = 20);
    ~NBLDPC();

    static std::map<std::tuple<int, int, int>, CodeData> codeDataMap;

    void encode();
    void decode(int iter, int NbOper, double offset);
};
