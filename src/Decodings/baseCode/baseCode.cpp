#include "baseCode.h"
#include <cstdlib> 

BaseCode::~BaseCode()
{
    delete[] codedBits;
    delete[] infoBits;
    delete[] decodedBits;
}


BaseCode::BaseCode()
{
    codedBits = nullptr;
    infoBits = nullptr;
    decodedBits = nullptr;
}

void BaseCode::check()
{
    newErrorBits = 0;
    newErrorFrames = 0;
    for (int i = 0; i < codeLength; i++)
    {
        if (codedBits[i] != decodedBits[i])
        {
            newErrorBits++;
        }
    }
    if (newErrorBits > 0)
    {
        newErrorFrames = 1;
    }
    errorBitsAll += newErrorBits;
    errorFramesAll += newErrorFrames;
}