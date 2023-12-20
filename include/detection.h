#pragma once

class Detection{
    private:
        static Detection * detection;
        Detection(int TxAntNum, int RxAntNum, int ModType, double SNRdB);

        void generateChannel();
        void generateTxSignals();
        void generateRxSignalsWithNoise();
    public:
        Detection(const Detection&) = delete;
        Detection& operator=(const Detection&) = delete;

        static void createDetection(int TxAntNum, int RxAntNum, int ModType, double SNRdB);
        static Detection * getDetection();


        int TxAntNum, RxAntNum, ModType, ConSize, bitLength;
        int TxAntNum2, RxAntNum2;
        double * TxSymbols;
        double * RxSymbols;
        const double * Cons;
        const double * Cons2;
        const int * bitCons;

        int * TxBits;
        double ** H;
        double SNRdB;
        double Nv;
        double sqrtNvDiv2;
        double NvInv;

        void reset();

        ~Detection()=default;
};