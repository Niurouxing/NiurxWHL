#pragma once


class DetectionRD;
class Detection;

// class DetectionAlgorithm, a virtual base class for all Mimo detection algorithms
class DetectionAlgorithm {
    protected:
        Detection * detection;
        int errorBits;
        int errorFrames;

        double  Nv;
        int TxAntNum, RxAntNum, TxAntNum2, RxAntNum2, ConSize, bitLength;

        int * TxBitsEst;
    
    public:
        DetectionAlgorithm();
        virtual void bind(Detection* detection)=0;
        virtual void execute() = 0;
        virtual ~DetectionAlgorithm()=default;
        virtual void check()=0;
        virtual void symbolsToBits(double * TxSymbolsEst)=0;

        int * getTxBitsEst() const { return TxBitsEst; }
        int getErrorBits() const { return errorBits; }
        int getErrorFrames() const { return errorFrames; }
};

class DetectionAlgorithmRD : public DetectionAlgorithm {
    protected:
        DetectionRD * detectionRD;
        double ** H;
        double * RxSymbols;
        const double * Cons;
        const int * bitCons;

    public:
        DetectionAlgorithmRD();
        void bind(Detection* detection) override;
        void check() override;
        void symbolsToBits(double * TxSymbolsEst) override;
        ~DetectionAlgorithmRD() override = default;

};


        


