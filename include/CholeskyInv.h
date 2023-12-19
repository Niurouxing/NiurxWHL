


class CholeskyInv
{
    private:
        double * L;
        double * y;
        int size;

    public:
        CholeskyInv(int size);
        void execute(double * A, double * AInv);
        ~CholeskyInv();
};