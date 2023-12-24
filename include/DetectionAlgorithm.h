#pragma once

#include <complex>
class DetectionRD;
class DetectionCD;
class Detection;

// class DetectionAlgorithm, a virtual base class for all Mimo detection algorithms
class DetectionAlgorithm {
    protected:
        Detection * detection;
        int errorBits;
        int errorFrames;

        double  Nv;
        double  NvInv;
        double  sqrtNvDiv2;
        int TxAntNum, RxAntNum, ConSize, bitLength;

        int * TxBitsEst;
    
    public:
        DetectionAlgorithm();
        virtual void bind(Detection* detection)=0;
        virtual void execute() = 0;
        virtual ~DetectionAlgorithm()=default;
        virtual void check()=0;

        int * getTxBitsEst() const { return TxBitsEst; }
        int getErrorBits() const { return errorBits; }
        int getErrorFrames() const { return errorFrames; }
};

class DetectionAlgorithmRD : public DetectionAlgorithm {
    protected:
        DetectionRD * detectionRD;
        int TxAntNum2, RxAntNum2;
        double ** H;
        double * RxSymbols;
        const double * Cons;
        const int * bitCons;

    public:
        DetectionAlgorithmRD();
        void bind(Detection* detection) override;
        void check() override;
        void symbolsToBits(double * TxSymbolsEst);
        ~DetectionAlgorithmRD() override;
};

class DetectionAlgorithmCD : public DetectionAlgorithm {
    protected:
        DetectionCD * detectionCD;
        std::complex<double> ** H;
        std::complex<double> * RxSymbols;

        const std::complex<double> * ConsComplex;
        const double * ConsReal;

        const int * bitConsComplex;
        const int * bitConsReal;

    public:
        DetectionAlgorithmCD();
        void bind(Detection* detection) override;
        void check() override;
        void symbolsToBits(std::complex<double> * TxSymbolsEst);
        ~DetectionAlgorithmCD() override;
};


        


