#pragma once

class Detection{
    private:
        static Detection * detection;
        Detection(int TxAntNum, int RxAntNum, int ModType, double SNRdB);


    public:
        Detection(const Detection&) = delete;
        Detection& operator=(const Detection&) = delete;

        static void createDetection(int TxAntNum, int RxAntNum, int ModType, double SNRdB);
        static Detection * getDetection();

        // used in both real and complex domain
        int TxAntNum, RxAntNum, ModType, ConSize, bitLength;
        int TxAntNum2, RxAntNum2;

        int * TxBits;
        double SNRdB;
        double Nv;
        double sqrtNvDiv2;
        double NvInv;


        // used in real domain
        double * TxSymbols;
        double ** H;
        double * RxSymbols;
        const double * Cons;
        const double * Cons2;
        const int * bitCons;





        void generateChannel();
        void generateTxSignals();
        void generateRxSignalsWithNoise();
        void generate();
        ~Detection()=default;
};