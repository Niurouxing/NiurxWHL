#pragma once

class BaseCode {
public:
    int codeLength;
    int infoLength;

    int * codedBits;
    int * infoBits;
    virtual void encode() = 0;
    ~BaseCode();
};


//  used in uncoded system
class plainCode : public BaseCode {
public:
    void encode() override;
};


 