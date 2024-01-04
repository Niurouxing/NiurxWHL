#include "baseCode.h"
#include <cstdlib> 

BaseCode::~BaseCode()
{
    delete[] codedBits;
    delete[] infoBits;
}

void plainCode::encode()
{
    for(int i=0;i<infoLength;i++){
        codedBits[i] = rand()%2;  
    }
}