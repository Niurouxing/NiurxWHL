#pragma once

class Mimo{
    private:
        static Mimo * mimo;
        Mimo(int TxAntNum, int RxAntNum, int ModType, double SNRdB);

        void generateChannel();
        void generateTxSignals();
        void generateRxSignalsWithNoise();
    public:
        Mimo(const Mimo&) = delete;
        Mimo& operator=(const Mimo&) = delete;

        static void createMimo(int TxAntNum, int RxAntNum, int ModType, double SNRdB);
        static Mimo * getMimo();


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

        ~Mimo()=default;
};