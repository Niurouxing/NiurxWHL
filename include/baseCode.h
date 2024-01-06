#pragma once

class BaseCode
{
protected:
    int errorBitsAll = 0;
    int errorFramesAll = 0;

    int newErrorBits = 0;
    int newErrorFrames = 0;

public:
    int codeLength;
    int infoLength;

    int *codedBits;
    int *infoBits;

    int *decodedBits;

    virtual void encode() = 0;
 
    int getErrorBits() const { return errorBitsAll; }
    int getErrorFrames() const { return errorFramesAll; }

    int getNewErrorBits() const { return newErrorBits; }
    int getNewErrorFrames() const { return newErrorFrames; }

    void resetError();

    virtual void check();

    BaseCode();
    ~BaseCode();
};
