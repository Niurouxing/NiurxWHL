#pragma once


#include <complex>

class Detection{
    public:

        Detection(int TxAntNum, int RxAntNum, int ModType, double SNRdB);

        // used in both real and complex domain
        int TxAntNum, RxAntNum, ModType, ConSize, bitLength;
        int TxAntNum2, RxAntNum2;

        int * TxIndice;
        int * TxBits;
        double SNRdB;
        double Nv;
        double sqrtNvDiv2;
        double NvInv;
        std::uniform_int_distribution<int> randomInt;

        virtual void generateChannel()=0;
        virtual void generateTxSignals()=0;
        virtual void generateTxSignals(int * bits)=0;
        virtual void generateRxSignalsWithNoise()=0;
        virtual void generate();
        virtual void generate(int * bits);
        virtual ~Detection();
};


class DetectionRD : public Detection{
    public:
        DetectionRD(int TxAntNum, int RxAntNum, int ModType, double SNRdB);

        // used in real domain
        double * TxSymbols;
        double * H;
        double * RxSymbols;

        const double * Cons;
        const double * Cons2;
        const int * bitCons;

        void generateChannel() override;
        void generateTxSignals() override;
        void generateTxSignals(int * bits) override;
        void generateRxSignalsWithNoise() override;

        ~DetectionRD() override;
};

class DetectionCD : public Detection{
    public:
        DetectionCD(int TxAntNum, int RxAntNum, int ModType, double SNRdB);

        // used in complex domain
        std::complex<double> * TxSymbols;
        std::complex<double> * H;
        std::complex<double> * RxSymbols;

        const double * ConsReal;
        const std::complex<double> * ConsComplex;

        const int * bitConsReal;
        const int * bitConsComplex;

        void generateChannel() override;
        void generateTxSignals() override;
        void generateTxSignals(int * bits) override;
        void generateRxSignalsWithNoise() override;

        ~DetectionCD() override;
};