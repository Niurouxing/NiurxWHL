#pragma once
 

class NBLDPCCode;
class NBLDPCTable;
class NBLDPCDecoder;

void GaussianElimination (NBLDPCCode * code, NBLDPCTable* table);

int Bin2GF(int *U,int GF,int logGF, int **BINGF=nullptr);

 

void Encoding(NBLDPCCode *code, NBLDPCTable *table, int *CodeWord, int **NBIN, int *KSYMB);

void EMS(NBLDPCCode *code, NBLDPCTable *table, NBLDPCDecoder *decoder, int NbOper, double offset, int NbIterMax, int *decide);
void CheckPassLogEMS (int node,NBLDPCDecoder *decoder, NBLDPCCode *code, NBLDPCTable *table,int NbOper, double offset);
void ElementaryStep(double *Input1,double *Input2,int *IndiceInput1,int *IndiceInput2,double *Output,int *IndiceOut,int **ADDGF,int GF,int nbMax,int nbOper);

 
void Decision(int* decision, double** APP, int N, int GF);
int Syndrom(NBLDPCCode* code, int* decide, NBLDPCTable* tableGF);