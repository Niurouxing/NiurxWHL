#pragma once
#include "Mimo.h"


// class Algorithm, a virtual base class for all Mimo detection algorithms
class Algorithm {
    protected:
        Mimo * mimo;
        int errorBits;
        int errorFrames;
    
    public:
        Algorithm();
        virtual void execute() = 0;
        virtual ~Algorithm() = default;
        void check();
        int getErrorBits() const { return errorBits; }
        int getErrorFrames() const { return errorFrames; }
};
        


